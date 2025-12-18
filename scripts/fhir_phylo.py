import json
import argparse
import os
import csv
import sys
import re

def parse_fhir(file_path):
    with open(file_path, 'r') as f:
        data = json.load(f)
    
    sample_id = os.path.basename(file_path).replace('.merged.fhir.json', '').replace('.fhir.json', '')
    cgmlst_dict = {}
    
    metadata = {
        "sample_id": sample_id,
        "patient_id": "NA",
        "latitude": "NA",
        "longitude": "NA",
        "conclusion": "NA",
        "mlst_st": "Unknown"
    }
    
    if 'entry' in data:
        for entry in data['entry']:
            res = entry.get('resource', {})
            r_type = res.get('resourceType')

            if r_type == 'Patient':
                metadata["patient_id"] = res.get('id', 'NA')
                
                for addr in res.get('address', []):
                    for ext in addr.get('extension', []):
                        if ext.get('url') == 'http://hl7.org/fhir/StructureDefinition/geolocation':
                            for sub_ext in ext.get('extension', []):
                                if sub_ext.get('url') == 'latitude':
                                    metadata["latitude"] = sub_ext.get('valueDecimal')
                                elif sub_ext.get('url') == 'longitude':
                                    metadata["longitude"] = sub_ext.get('valueDecimal')

            elif r_type == 'DiagnosticReport':
                conclusions = []
                if res.get('conclusion'):
                    conclusions.append(res.get('conclusion'))
                
                for cc in res.get('conclusionCode', []):
                    if cc.get('text'):
                        conclusions.append(cc.get('text'))
                
                if conclusions:
                    metadata["conclusion"] = "; ".join(list(set(conclusions)))

            elif r_type == 'Observation':
                code_obj = res.get('code', {})
                
                is_mlst = False
                for coding in code_obj.get('coding', []):
                    if coding.get('code') == '612-2':
                        is_mlst = True
                        break
                
                if is_mlst:
                    st_val = res.get('valueCodeableConcept', {}).get('text')
                    if not st_val:
                        codings = res.get('valueCodeableConcept', {}).get('coding', [])
                        if codings:
                            st_val = codings[0].get('display')
                    
                    if st_val:
                        match = re.search(r"(ST\d+)", st_val)
                        if match:
                            metadata['mlst_st'] = match.group(1)
                        else:
                            metadata['mlst_st'] = st_val
                    continue 

                if metadata['mlst_st'] == "Unknown":
                    div_text = res.get('text', {}).get('div', '')
                    match = re.search(r"\b(ST\d+)\b", div_text)
                    if match:
                        metadata['mlst_st'] = match.group(1)
                    else:
                        match = re.search(r"\bST[:\s]+(\d+)\b", div_text, re.IGNORECASE)
                        if match:
                            metadata['mlst_st'] = f"ST{match.group(1)}"

                is_cgmlst = False
                for coding in code_obj.get('coding', []):
                    if coding.get('code') == '95643-1':
                        is_cgmlst = True
                        break
                
                if not is_cgmlst and code_obj.get('text') == 'Core Genome MLST Profile':
                    is_cgmlst = True

                if is_cgmlst:
                    for i, comp in enumerate(res.get('component', [])):
                        if 'valueQuantity' in comp: continue

                        locus_name = f"Locus_{i+1}"
                        if 'code' in comp:
                            if 'coding' in comp['code'] and len(comp['code']['coding']) > 0:
                                locus_name = comp['code']['coding'][0].get('code', locus_name)
                            elif 'text' in comp['code']:
                                locus_name = comp['code']['text']
                        
                        val = comp.get('valueString', '-')
                        if val in ['NA', '', 'LNF', 'NIPH', 'ALM', 'NIPHEM']:
                            val = '-'
                        
                        cgmlst_dict[locus_name] = val
                        
    return sample_id, cgmlst_dict, metadata

def calculate_distance(p1_dict, p2_dict, all_loci):
    dist = 0
    missing_values = {'-'}
    for locus in all_loci:
        a1 = p1_dict.get(locus, '-')
        a2 = p2_dict.get(locus, '-')
        if a1 in missing_values or a2 in missing_values: continue
        if a1 != a2: dist += 1
    return dist

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', nargs='+', required=True, help="List of FHIR JSON files")
    args = parser.parse_args()

    sample_data = {}
    all_metadata = []
    samples = []
    all_loci = set()

    for f in args.inputs:
        sid, profile_dict, meta = parse_fhir(f)
        if not profile_dict:
            print(f"Warning: No cgMLST found for {sid}", file=sys.stderr)
            continue
        sample_data[sid] = profile_dict
        all_metadata.append(meta)
        samples.append(sid)
        all_loci.update(profile_dict.keys())

    if not samples:
        sys.exit(1)

    sorted_loci = sorted(list(all_loci))

    with open("metadata.tsv", "w", newline='') as f:
        fieldnames = ["sample_id", "patient_id", "latitude", "longitude", "conclusion", "mlst_st"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(all_metadata)

    with open("cgmlst_profile.tsv", "w", newline='') as f:
        header = ["#Names"] + sorted_loci
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(header)
        for sid in samples:
            row = [sid]
            for locus in sorted_loci:
                row.append(sample_data[sid].get(locus, '-'))
            writer.writerow(row)

    n = len(samples)
    matrix_data = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            s1 = samples[i]
            s2 = samples[j]
            dist = calculate_distance(sample_data[s1], sample_data[s2], sorted_loci)
            matrix_data[i][j] = dist
            matrix_data[j][i] = dist
            
    with open("distance_matrix.tsv", "w") as f:
        f.write("cgmlst-dists\t" + "\t".join(samples) + "\n")
        for i, row in enumerate(matrix_data):
            f.write(samples[i] + "\t" + "\t".join(map(str, row)) + "\n")

if __name__ == "__main__":
    main()