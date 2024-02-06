nextflow.enable.dsl=2

workflow {
    comparison_conditions_channel = Channel.fromList(params.comparison_condtitions)
    comparison_conditions.map { folder_name -> conditions_compared = folder_name.conditions_compared }
    other_keys_channel = Channel.from{params_file.other_keys}

    output_folder()

    write_params()

    building_dds_obj(params.count_data, params.metadata, comparison_conditions_channel)

    if(!perform_batch_correction){
        batch_correction()
    }

    if(shrink_LFC_using || MA){
        shrink_LFC()
    }

    if(list_all_genes || list_DEG){
        detection_DEG()
    }

    if(mean_vs_sd_before_stabilisation || plot_sample_clustering){
        perform_variance_stabilisation()
    }

    if(PCA){
        plot_PCA()
    }

    if(counts){
        plot_counts()
    }
}



process output_folder{
    input:
    val conditions_compared from comparison_conditions
    
    output:
    val conditions_compared into output_folder
    
    script:
    """
    mkdir $conditions_compared
    mv $PWD/output_folder $conditions_compared
    """
}

process write_params{
    output:
    path "*.txt"
    
    shell:
    '''
    suffix=$(date "+%Y_%m_%d_%H_%M")
    echo "Paramenters of the job on" $suffix >> params_used_${suffix}.txt
    echo "!{params}" >> params_used_${suffix}.txt
    '''
}


process building_dds_obj{
    publishDir params.outdir, mode: 'copy', overwrite: true
    input:
    path count_data
    path metadata
    tuple path(count_data), path(metadata), val(comparison_key)
    output:
    path "*.pdf"
    path "dds_obj.rds", emit: dds_obj

    shell:
    '''
    building_deseq2_object.r --count_data=!{count_data} --metadata=!{metadata} --comparison_key=!{comparison_key} --other_keys=!{params.other_keys} --minCounts=!{params.minCounts} --smallestGroupSize=!{params.minSamples} --detect_sample_outliers=!{params.detect_sample_outliers} --RLE_threshold_max=!{params.RLE_threshold_max}
    '''
}


process group_comparison{
    publishDir params.outdir, mode: 'copy', overwrite: true
    input:
    tuple path("dds_obj.rds"), val(conditions_compared), val(pValue)

    output:
    file 

    script:
    """
    Rscript group_comparison.r --conditions_compared=${conditions_compared} --pValue=${pValue}
    """
}

process batch_correction{
    publishDir params.outdir, mode: 'copy', overwrite: true
    input:
    tuple path("dds_obj.rds"), path("dds_design.rds"), path("dds_results.rds"), val(comparison_key), val(RUV_threshold_not_sig), val(conditions_compared)

    output:
    file pdf_"*"_plot

    script:
    """
    Rscript batch_correction_deseq2.r --comparison_key=${comparison_key} --RUV_threshold_not_sig=${RUV_threshold_not_sig} --conditions_compared=${conditions_compared}
    """
}

process shrink_LFC{
    input:
    tuple path("dds_obj.rds"), path("dds_results.rds"), val(conditions_compared), val(MA), val(shrink_LFC_using)

    output:
    file pdf_MA
    
    script:
    """
    Rscript LFCshrink_deseq2.r --conditions_compared=${conditions_compared} --MA=${MA} --shrink_LFC_using=${shrink_LFC_using}
    """
}

process detection_DEG{
    publishDir params.outdir, mode: 'copy', overwrite: true
    input:
    tuple path("dds_results.rds"), val(condition_compared), val(pValue), val(log2fc_threshold_DEG), val(list_all_genes), val(list_DEG)

    output:
    file xlsx_allGenes
    file xlsx_DEG


    script:
    """
    Rscript DEG_detection_deseq2.r --condition_compared=${condition_compared} --pValue=${pValue} --log2fc_threshold_DEG=${log2fc_threshold_DEG} --list_all_genes=${list_all_genes} --list_DEG=${list_DEG}
    """
}

process perform_variance_stabilisation{ //needs optimizing if other_keys are given
    input:
    tuple path("dds_results.rds"), val(perform_variance_stabilisation), val(do_blind_stabilization), val(mean_vs_sd_before_stabilisation), val(mean_vs_sd_after_stabilisation), val(plot_sample_clustering), val(other_keys)

    output:
    file pdf_transform
    file pdf_heatmap

    script:
    """
    Rscript variance_stabilisation_deseq2.r --perform_variance_stabilisation=${perform_variance_stabilisation} --do_blind_stabilization=${do_blind_stabilization} --mean_vs_sd_before_stabilisation=${mean_vs_sd_before_stabilisation} --mean_vs_sd_after_stabilisation=${mean_vs_sd_after_stabilisation} --plot_sample_clustering=${plot_sample_clustering} --other_keys=${other_keys}
    """
}

process plot_PCA{ //needs optimizing if other_keys are given 
//I don't know why, but it needs another round of pre-filtering
    input:
    tuple path("dds_results.rds"), val(other_keys), val(PCA), val(minCounts), val(smallestGroupSize)

    output:
    file pdf_gpca

    script:
    """
    Rscript PCAplot_deseq2.r --other_keys=${other_keys} --PCA=${PCA} --minCounts=${minCounts} --smallestGroupSize=${smallestGroupSize}
    """
}

//needs optimisation
process plot_counts{
    input:
    tuple path("dds_results.rds"), val(comparison_key), val(genes_of_interest), val(counts)

    output:
    file pdf_geneCount

    script:
    """
    Rscript counts_plot_deseq2.r --comparison_key${comparison_key} --genes_of_interest${genes_of_interest} --counts=${counts}
    """
}
