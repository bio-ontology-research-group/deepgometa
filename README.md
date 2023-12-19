# deepgometa
This repository contains the scripts and datafiles used in the DeepGOmeta manuscript.

## Paired Datasets
1. **Data and metadata**: download from SRA and MG-RAST using [sample accessions](PairedDatasets/Sample_data.csv)
2. **Processing reads**:
   * 16S reads - generate OTU tables using the Nextflow [16SProcessing workflow](https://github.com/bio-ontology-research-group/16SProcessing)
   * WGS reads - obtain protein sequences using the [assembly pipeline](PairedDatasets/WGSPipeline.py)
