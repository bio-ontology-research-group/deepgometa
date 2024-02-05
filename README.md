# DeepGOMeta
This repository contains the scripts and datafiles used in the DeepGOmeta manuscript.

# Dependencies
* The code was developed and tested using python 3.10.
* Clone the repository: `git clone https://github.com/bio-ontology-research-group/deepgometa.git`
* Create virtual environment with Conda or python3-venv module. 
* Install PyTorch: `pip install torch==2.0.1 torchvision==0.15.2 torchaudio==2.0.2`
* Install DGL: `pip install dgl==1.1.2+cu117 -f https://data.dgl.ai/wheels/cu117/repo.html`
* Install other requirements:
  `pip install -r requirements.txt`


# Running DeepGOMeta model
Follow these instructions to obtain predictions for your proteins. You'll need
around 30Gb storage and a GPU with >16Gb memory (or you can use CPU)
* Download the [data.tar.gz](https://deepgo.cbrc.kaust.edu.sa/data/deepgometa/data.tar.gz)
* Extract `tar xvzf data.tar.gz`
* Run the model `python predict.py -if data/example.fa`


# Docker container
We also provide a docker container with all dependencies installed:
`docker pull coolmaksat/deepgometa` \
This repository is installed at /deepgometa directory. To run the scripts you'll
need to mount the data directory. Example: \
`docker run --gpus all -v $(pwd)/data:/workspace/deepgometa/data coolmaksat/deepgometa python predict.py -if data/example.fa`

# Nextflow
DeepGOMeta can be run as a Nextflow workflow using the docker image for easier execution.

Requirements:
* For amplicon data: OTU table of relative abundance, where OTUs are classified using the RDP database
* For WGS data: Protein sequences in FASTA format 

1. After cloning the repository, navigate to the Nextflow directory: `cd Nextflow`
2. Update the runOptions paths in [nextflow.config](Nextflow/nextflow.config)
3. Navigate to the data directory `cd data` and download the [genome annotations](https://bio2vec.cbrc.kaust.edu.sa/data/deepgometa/rdp_genomes.tar.gz)
4. Run workflow. Example: `nextflow run DeepGOMeta.nf -profile docker/singularity --amplicon true --OTU_table otu_relative_abd.tsv --pkl_dir /PATH/TO/PKL/DIR/`

# Paired Datasets
1. **Data and metadata**: download from SRA and MG-RAST using [sample accessions](PairedDatasets/Sample_data.csv)
2. **Processing reads**:
   * 16S reads - generate OTU tables using the Nextflow [16SProcessing workflow](https://github.com/bio-ontology-research-group/16SProcessing)
   * WGS reads - obtain protein sequences using the [assembly pipeline](PairedDatasets/WGSPipeline.py)
3. **Functional annotation**:
   * OTU tables - [generate](PairedDatasets/16S_function_abundance.py) a weighted functional profile for each OTU table using DeepGOmeta predictions
   * Protein fasta - [run]() DeepGOmeta on Prodigal output from metagenome assemblies, and [generate](PairedDatasets/WGS_function_binary.py) a binary functional profile for each dataset
4. **Clustering and Purity**: use a metadata file and the functional profile to apply PCA, k-means clustering, calculating purity, and generating plots for [16S datasets](PairedDatasets/16S_pheno_PCA.py) and [WGS datasets](PairedDatasets/WGS_pheno_PCA.py)
5. **Information Content Calculation**: create a .txt file for each sample containing the 16S predicted functions and WGS predicted functions on separate lines (e.g. 16Ssample'\t'GO1'\t'GO2'\n'WGSsample'\t'GO2'\t'GO3), and get [IC](PairedDatasets/ICVectorSim.groovy) for each function, then run a [t-test](PairedDatasets/t-test_IC.py)
