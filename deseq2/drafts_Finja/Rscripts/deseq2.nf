workflow {
    comparison_conditions_channel = Channel.fromList(params.comparison_condtitions)
    comparison_conditions.map { folder_name -> conditions_compared = folder_name.conditions_compared }
    other_keys_channel = Channel.from{params_file.other_keys}

    output_folder()

    write_params()

    building_object()

    if(!perform_batch_correction){
        batch_correction()
    }

    plotting()
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
    input:

    output:

    script:
}


//COLORING !!!!!
process building_DDS_obj {
    publishDir params.outdir, mode: 'copy', overwrite: true
    input:
    file(params_file) from params.params_file

    tuple OTHER_KEYS, path(count_data), path(metadata), val(comparison_key), val(reference_condition), val(alternate_condition), val(conditions_compared)
    val(other_keys) from other_keys_channel

    output:
    file pdf_dispersion_plot
    file pdf_pValue_plot

    script:
    """
    Rscript building_object_deseq2.r --params_file=${params_file} --count_data=${count_data} --metadata=${metadata} --comparison_key=${comparison_key} --reference_condition=${reference_condition} --alternate_condition=${alternate_condition} --conditions_compared=${conditions_compared}
    """
}

process batch_correction{
    publishDir params.outdir, mode: 'copy', overwrite: true
    input:
    

    output:
    file pdf_ruv_plot

    script:
    """
    Rscript batch_correction_deseq2.r 
    """
}


process plotting{
    publishDir params.outdir, mode: 'copy', overwrite: true
    input:
    

    output:
    file xlsx_allGenes
    file xlsx_DEG


    script:
    """
    Rscript plots_deseq2.r
    """
}

