#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define the parameters that will be passed to Nextflow

params.output = 'results'
params.variant_mapping = ''

variant_mapping_ch = Channel.fromPath(params.variant_mapping).collect()
chain_file_ch = Channel.fromPath(params.chain_file).collect()
format_ch = Channel.fromPath(params.format_file).collect()

trait_file = params.trait_file

process HARMONIZE_GWAS {
        //scratch true // Needs to be set to true!
        scratch '$TMPDIR'
        // scratch false
        publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: 'eqtlgen_parquet_dataset/gwas_id=*'
        tag "$id"
        beforeScript "ml R"

    input:
        tuple val(id), val(format_id), path(gwas), val(sample_size)
        path variant_reference
        path chain_file
        path format

    output:
        tuple val(id), path("eqtlgen_parquet_dataset/gwas_id=*")
        path("*.log.txt"), emit: log

    shell:
        '''
        harmonize_sumstats.R \
        --gwas-file !{gwas.join(" ")} \
        --gwas-id !{id} \
        --format-file !{format} \
        --format-id !{format_id} \
        --out-dir . \
        --liftover-chain !{chain_file} \
        --eQTLGen-index !{variant_reference} \
        --sample-size !{sample_size}
        '''
}


workflow {
    gwas_ch = Channel.fromPath(trait_file)
        .splitCsv(header: true, sep: '\t')
        .map { gwas ->
            def full_paths = file("${gwas.gwas_dir}/${gwas.gwas_file}") // Join and convert to file
            return [trait_identifier: gwas.trait_identifier, format_name: gwas.format_name, gwas_path: full_paths, sample_size: gwas.sample_size]
        }
        .branch { gwas ->
            existing: file(gwas.gwas_path[0]).exists()
            other: true
        }

    // harmonized_ch = HARMONIZE_GWAS(combined_ch, format_ch)
    harmonization_output = HARMONIZE_GWAS(gwas_ch.existing, variant_mapping_ch, chain_file_ch, format_ch)
    harmonization_output.log.collectFile(name: 'harmonization_log.txt', keepHeader: true, skip: 1, storeDir: params.output)

    // Write list of missing files
    gwas_ch.other
        .map { gwas -> gwas.gwas_path.toString() }
        .collectFile(name: "missing_files.txt", newLine: true, storeDir: params.output)
}
