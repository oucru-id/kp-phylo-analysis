nextflow.enable.dsl = 2

process VISUALIZE_REPORT {
    publishDir "${params.results_dir}/visualization", mode: 'copy'
    
    input:
    path matrix
    path metadata
    path tree

    output:
    path "stats_histogram.png",       emit: hist
    path "stats_heatmap.png",         emit: heatmap
    path "stats_violin.png",          emit: violin
    path "grapetree_viz_rectangular.png", emit: rect_tree
    path "grapetree_viz_radial.png",      emit: radial_tree

    script:
    """
    python3 $baseDir/scripts/visualize_results.py \\
        --matrix ${matrix} \\
        --metadata ${metadata} \\
        --tree ${tree}
    """
}

workflow VISUALIZATION {
    take:
    matrix
    metadata
    tree

    main:
    VISUALIZE_REPORT(matrix, metadata, tree)
}