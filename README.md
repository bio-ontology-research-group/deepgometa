# deepgometa
This repository contains the scripts and datafiles used in the DeepGOmeta manuscript.

## Paired Datasets
1. **Data and metadata**: download from SRA and MG-RAST using [sample accessions](PairedDatasets/Sample_data.csv)
2. **Processing reads**:
   * 16S reads - generate OTU tables using the Nextflow [16SProcessing workflow](https://github.com/bio-ontology-research-group/16SProcessing)
   * WGS reads - obtain protein sequences using the [assembly pipeline](PairedDatasets/WGSPipeline.py)
3. **Functional annotation**:
   * OTU tables - [generate](PairedDatasets/16S_function_abundance.py) a weighted functional profile for each OTU table using DeepGOmeta predictions
   * Protein fasta - [run]() DeepGOmeta on Prodigal output from metagenome assemblies, and [generate](PairedDatasets/WGS_function_binary.py) a binary functional profile for each dataset
4. **Clustering and Purity**: 
