# KP FHIR Phylogeny Analysis Pipeline

This pipeline processes FHIR bundle JSON files containing Klebsiella pneumoniae genomics data to generate cgMLST matrices and phylogenetic trees. Please refer to our full documentation: https://kp-pipeline-docs.readthedocs.io/en/latest/


## Features

- **Direct FHIR Genomics JSON input**
- **Phylogenetic tree**
- **Transmission network visualization:** Graph and statistical plots (histogram, heatmap, violin plot).
- **Clinical metadata integration:** Extracted from Bundle Genomics FHIR.

## Usage

### Requirements

- [Nextflow](https://www.nextflow.io/)
- Python 3.8+
- Python packages: `biopython`, `pandas`, `pyvis`, `matplotlib`, `seaborn`

Install Python dependencies:
```bash
pip install biopython pandas networkx pyvis matplotlib seaborn
```

### Run the Pipeline

```bash
nextflow run main.nf
```

### Input

- Place FHIR Bundle Genomics JSON files in `data/JSON/`
