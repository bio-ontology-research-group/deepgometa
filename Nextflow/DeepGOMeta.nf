#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BORG/DeepGOMeta
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/bio-ontology-research-group/deepgometa
----------------------------------------------------------------------------------------
*/

params.OTU_table = 'otutab_relative_withtaxa_merged.tsv' // OTU table of relative abundance and assigned taxa using RDP
params.protein_fasta = '' // Prodigal output of amino acid sequences in Fasta format
params.out_dir = 'results' // Default output directory
params.amplicon = false
params.wgs = false
params.python_paths = "/data_and_scripts/"
params.pkl_dir = "/home/tawfiqre/standardize_methods/DeepGOMeta/data_and_scripts/rdp_genomes_preds/"

/*
    ================================================================================
                                    16S Data
    ================================================================================
*/

process ExtractAccessions {
    label "amplicon"
    publishDir "${params.out_dir}"

    input:
    path otu_table
    path rdp_genomes_clean

    output:
    path "extracted_accessions.tsv"

    script:
    """
    python3 ${params.python_paths}extract_accessions.py \
        --otu_table $otu_table \
        --rdp_genomes_clean $rdp_genomes_clean \
        --output extracted_accessions.tsv
    """
}

process FunctionAbundance {
    label "amplicon"
    publishDir "${params.out_dir}"

    input:
    path extractedAccessions
    path otu_table

    output:
    path "GO_terms.csv"
    path "OTU_function_abundance.tsv"

    script:
    """
    python3 ${params.python_paths}16s_function_abundance.py \
        --pkl_dir ${params.pkl_dir} \
        --accessions $extractedAccessions \
        --otu_table $otu_table
    """
}
/*
    ================================================================================
                                    WGS Data
    ================================================================================
*/
process RunDeepGOMetaModel {
    label "WGS"
    input:
    path fasta_file

    output:
    path "predictions.*"

    script:
    """
    python predict.py -if ${fasta_file}
    """
}



workflow {
    rdp_genomes_clean = Channel.fromPath("${params.pkl_dir}/rdp-genomes_clean.txt")

    if (params.amplicon) {
    otu_table = Channel.fromPath(params.OTU_table)
    extractedAccessions = ExtractAccessions(otu_table, rdp_genomes_clean)
    funAbundOut = FunctionAbundance(extractedAccessions, otu_table)
    }

    if (params.wgs) {
    
    }
}
