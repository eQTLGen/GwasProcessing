#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

params.output = 'results'
params.variant_mapping = ''

variant_mapping_ch = Channel.fromPath(params.variant_mapping).collect()

// pQTL arguments
pqtl_files_ch = Channel.fromPath("${params.pqtl_files}/*.tar").
    map { file ->
          def fileSplit = file.name.toString().split('_')
          def assay = fileSplit[0]
          def uniprot = fileSplit[1]
          def olinkId = fileSplit[2]

          def key = "${assay}_${uniprot}_${olinkId}"

          return tuple(key, file) }


process HARMONIZE_UKB_PPP {
    //scratch true // Needs to be set to true!
    scratch '$TMPDIR'
    // scratch false
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: 'eqtlgen_parquet_dataset/gwas_id=*'
    tag "$id"
    beforeScript "ml R"

    input:
        tuple val(id), path(pqtl)
        path variant_reference

    output:
        tuple val(id), path("eqtlgen_parquet_dataset/gwas_id=*")
        path("*.log.txt")

    shell:
        '''
        mkdir untar
        tar -C untar -xvf !{pqtl}
        pqtl_path=!{pqtl}
        pqtl_path=$(echo "${pqtl_path%.*}")

        mkdir output

        harmonize_pqtl.R \
        --pqtl-folder untar/${pqtl_path} \
        --pqtl-id !{id} \
        --eQTLGen-index !{variant_reference} \
        --out-dir .

        rm -r untar
        '''
}

workflow {
    HARMONIZE_UKB_PPP(pqtl_files_ch, variant_mapping_ch)
}
