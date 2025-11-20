#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define the parameters that will be passed to Nextflow

params.output = 'results'
params.variant_mapping = ''

variant_mapping_ch = Channel.fromPath(params.variant_mapping).collect()
sqtl_path_ch = Channel.fromPath(params.sqtl_path)

process HARMONIZE_SQTL {
        //scratch true // Needs to be set to true!
        scratch '$TMPDIR'
        // scratch false
        publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: 'eqtlgen_parquet_dataset/phenotype_id=*'
        tag "$id"
        beforeScript "ml R"

    input:
        path gwas
        path variant_reference

    output:
        tuple val(id), path("eqtlgen_parquet_dataset/phenotype_id=*")
        path("*.log.txt"), emit: log

    shell:
        '''
        harmonize_sqtl.R \
        --gwas-file !{gwas.join(" ")} \
        --out-dir . \
        --eQTLGen-index !{variant_reference}
        '''
}


workflow {
    // harmonized_ch = HARMONIZE_GWAS(combined_ch, format_ch)
    harmonization_output = HARMONIZE_SQTL(sqtl_path_ch, variant_mapping_ch)
    harmonization_output.log.collectFile(name: 'harmonization_log.txt', keepHeader: true, skip: 1, storeDir: params.output)
}
