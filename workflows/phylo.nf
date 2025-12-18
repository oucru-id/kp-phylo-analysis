nextflow.enable.dsl = 2

process FHIR_ANALYSIS {
    publishDir "${params.results_dir}/phylo", mode: 'copy'
    
    input:
    path fhir_files

    output:
    path "distance_matrix.tsv", emit: matrix
    path "cgmlst_profile.tsv",  emit: profile
    path "metadata.tsv",        emit: metadata

    script:
    """
    python3 $baseDir/scripts/fhir_phylo.py \\
        --inputs ${fhir_files}
    """
}

process GRAPETREE_MST {
    publishDir "${params.results_dir}/phylo", mode: 'copy'

    input:
    path profile

    output:
    path "grapetree.nwk", emit: tree

    script:
    """
    grapetree \
        --profile ${profile} \
        --method MSTreeV2 \
        > grapetree.nwk
    """
}

workflow PHYLO_ANALYSIS {
    take:
    fhir_files

    main:
    FHIR_ANALYSIS(fhir_files.collect())
    
    GRAPETREE_MST(
        FHIR_ANALYSIS.out.profile
    )

    emit:
    matrix   = FHIR_ANALYSIS.out.matrix
    profile  = FHIR_ANALYSIS.out.profile
    tree     = GRAPETREE_MST.out.tree
    metadata = FHIR_ANALYSIS.out.metadata
}