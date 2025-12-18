#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info """
    Klebsiella pneumoniae (cgMLST) Phylogeny & Visualization Pipeline (v${params.version})
    Created by: SPHERES Lab Team
"""

include { PHYLO_ANALYSIS } from './workflows/phylo.nf'
include { VISUALIZATION }  from './workflows/visualization.nf'
include { VERSIONS }  from './workflows/utils.nf'

workflow {
    fhir_ch = Channel.fromPath("${params.fhir_dir}/*.json", checkIfExists: true)
    
    PHYLO_ANALYSIS(fhir_ch)
    
    VISUALIZATION(
        PHYLO_ANALYSIS.out.matrix, 
        PHYLO_ANALYSIS.out.metadata, 
        PHYLO_ANALYSIS.out.tree
    )
    
    VERSIONS()
}