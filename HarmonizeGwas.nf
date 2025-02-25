#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define the parameters that will be passed to Nextflow

Channel.fromPath(params.gwas_tab_file).collect().set { gwas_file_ch }


process JOIN_GWAS_TO_FORMAT {
    cpus 1
    memory '4GB'

    input:
    path gwas
    path format

    output:
    path "merged.txt"

    script:
    """
    ml purge
    ml RPlus/4.2.1-foss-2022a-v22.12.1
    export R_LIBS_USER=/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/users/umcg-rwarmerdam/Rlibs/x86_64-pc-linux-gnu-library/4.2

    merge_gwas_to_format.R --trait-file "gwas" --format-file "format"
    """
}

process HARMONIZE_GWAS {
    tag "$gwas"
    cpus 1
    memory '4GB'

    input:
    path gwas
    path format

    output:
    path "output"

    script:
    """
    ml purge
    ml RPlus/4.2.1-foss-2022a-v22.12.1
    export R_LIBS_USER=/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/users/umcg-rwarmerdam/Rlibs/x86_64-pc-linux-gnu-library/4.2

    mkdir output

    /groups/umcg-fg/tmp01/projects/eqtlgen-phase2/users/umcg-rwarmerdam/eqtlgen-phase2_scripts/harmonize_sumstats.R \
    --gwas-file ${gwas} \
    --format-file ${format} \
    --out-dir "output" \
    --liftover-chain /groups/umcg-fg/tmp01/projects/eqtlgen-phase2/public_data/liftOver/hg19ToHg38.over.chain \
    --eQTLGen-index /groups/umcg-fg/tmp01/projects/eqtlgen-phase2/processed_data/variants/1000G-30x_index.parquet \
    --z
    """
}

process HARMONIZE_UKB_PPP {
    tag "$gwas"
    cpus 2
    memory '6GB'

    input:
        path pqtl
        tuple

    output:
        tuple val(gwas), path("${gwas}_processed")

    shell:
        '''
        mkdir untar
        tar -C untar -xvf !{pqtl}
        pqtl_path=!{pqtl}
        pqtl_path=$(echo "${pqtl_path%.*}")

        harmonize_pqtl.R \
        --pqtl_folder untar/${pqtl_path} \
        --eQTLGen-index /groups/umcg-fg/tmp01/projects/eqtlgen-phase2/processed_data/variants/1000G-30x_index.parquet
        --trait-id !{id}

        rm -r untar
        '''
}


workflow {
    gwas_ch = Channel.fromPath(trait_file).collect().splitCsv(header: true).map { gwas ->
      return tuple(trait_name, file)
      }
    // format_ch = Channel.fromPath(format_file).collect()

    gwas_ch.branch { gwas ->
        ukb_ppp: gwas.format_name == "ukb_ppp"
        standard: gwas.format_name != "ukb_ppp"
    }.set { ukb_ppp_ch, standard_ch }

    // harmonized_ch = HARMONIZE_GWAS(combined_ch, format_ch)
    processed_ukb_ch = HARMONIZE_UKB_PPP(ukb_ppp_ch)

}