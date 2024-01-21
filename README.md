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
4. **Clustering and Purity**: use a metadata file and the functional profile to apply PCA, k-means clustering, calculating purity, and generating plots for [16S datasets](PairedDatasets/16S_pheno_PCA.py) and [WGS datasets](PairedDatasets/WGS_pheno_PCA.py)
5. **Information Content Calculation**: create a .txt file for each sample containing the 16S predicted functions and WGS predicted functions on separate lines (e.g. 16Ssample'\t'GO1'\t'GO2'\n'WGSsample'\t'GO2'\t'GO3), and get [IC](PairedDatasets/ICVectorSim.groovy) for each function, then run a [t-test](PairedDatasets/t-test_IC.py)
