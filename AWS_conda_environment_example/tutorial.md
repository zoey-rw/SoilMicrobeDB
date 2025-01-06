## Using SoilMicrobeDB to classify metagenomic soil samples

SoilMicrobeDB is a Kraken2 database with extensive representation of high-quality genomes of soil organisms, including uncultured and fungal species. Key applications of SoilMicrobeDB include advancing metagenomic research in soil ecology, identifying microbial taxa in diverse soil environments, and supporting agricultural and environmental monitoring. For details on how this database was designed, and potential applications, please see our manuscript (Werbin et al., 2025).

The tutorial below is for locally downloading this database to classify samples or update the database with new genomes. See section 6 ("Cloud Computing") for an alternative approach using an Amazon Machine Image. This tutorial walks through the process for evaluating the relative abundance of species within a metagenomic soil sample sequenced by the National Ecological Observatory Network.  

Downloading the Kraken2 database (Step 1) requires 289 GB of storage, and a similarly large amount of RAM for loading the database into memory when classifying samples (Step 2). Steps 3 and 4 are not RAM-intensive, but require installing software other than Kraken2. Step 4 requires downloading additional data files (approximately 400 MB of storage).

The simplest approach to setup is to leverage the Struo2 pipeline for database development (Youngblut et al. 2022), with full installation details available here. This includes the software for steps 2 and 4, but the quality-filtering software (Architeuthis) must be installed separately. 

### 0. Install software, define paths and environment variables

#### 0.1 Installation
To set up the Struo2 environment, run the following commands:

```
git clone --recurse-submodules https://github.com/leylabmpi/Struo2.git 
cd Struo2
conda env create --name struo2_conda -f conda_env.yaml
conda activate struo2_conda
mamba install kraken2
```

Install quality-filter software (Architeuthis)
```
conda install -c conda-forge -c bioconda architeuthis
```

Install helper for quality-filter outputs
```
wget https://github.com/daydream-boost/kraken2/archive/refs/heads/master.zip
unzip master.zip
cd ./kraken2-master/src
g++ -O3 -std=c++11 mmap_file.cc reports.cc taxonomy.cc kraken2-report.cpp -o kraken2-report
```

Install Bracken
Additional details on Bracken software available from the Bracken repository. 

```
wget https://github.com/jenniferlu717/Bracken/archive/refs/heads/master.zip
unzip master.zip -d ./Bracken
bash ./Bracken/Bracken-master/install_bracken.sh
```

#### 0.2 Specify input and output data paths

To ensure setup is completed smoothly, we provide a single "test" sample in the project repository, which contains forward (R1) and reverse (R2) metagenomic reads. Download the test samples and assign them to variables using the following code:

```
# Placeholder - replace with test files once AWS URL is available
wget -P ./test_sample ./test_sample_R1.fastq.gz ./test_sample_R2.fastq.gz

samp_dir=./test_sample
samp_name="$(basename -- $samp_dir)"
path_samp_prefix=${samp_dir}/${samp_name}
READ1=${path_samp_prefix}_R1.fastq.gz
READ2=${path_samp_prefix}_R2.fastq.gz

# Create output directories
mkdir -p ./test_output/{01_kraken,02_architeuthis,03_bracken}

input_dir_path=./test_sample
kraken_dir_path=./test_output/01_kraken
architeuthis_dir_path=./test_output/02_architeuthis
bracken_dir_path=./test_output/03_bracken

```
### 1. Downloading the SoilMicrobeDB:

Download using AWS command line toolkit:

```
# Placeholder - replace with SoilGenomeDB once AWS URL is available
aws s3 ls --no-sign-request s3://kraken2-ncbi-refseq-complete-v205/
DBNAME=SoilMicrobeDB
```

After download, ensure each of the following paths points to a location on your system:
```
# Placeholder - replace with SoilGenomeDB file paths once AWS URL is available
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxonomy/
DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxo.k2d
```

### 2. Classifying a sample

To classify a sample, Kraken2 must be installed (Step 0). Modify the value of "NTHREADS" to reflect how many threads should run in parallel on your system, and specify variables to reflect the input and output locations. 

```
NTHREADS=28 
KRAKEN_REPORT=${kraken_dir_path}/${samp_name}_${DBNAME}_kraken.kreport
KRAKEN_OUTPUT=${kraken_dir_path}/${samp_name}_${DBNAME}_kraken.output
```
Then, classify a sample with the following command:
```
kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --threads $NTHREADS --report-zero-counts --memory-mapping --report-minimizer-data --paired $READ1 $READ2
```

Note: the "--memory-mapping" checks for the database in memory rather than loading it from disk space. However, it is not suitable for all workflows. If you are using a high-performance cluster, check with system administrators to prevent excessive use of shared memory.

### 3. Quality-filtering 

To run the quality-filter step, we recommend the software Architeuthis (Diener et al. 2024). Install using the instructions provided by Architeuthis. The output files provide per-read results, but must be converted to the Kraken2 report format before passing to Bracken (Step 4). Currently, the only way to do this is by compiling a C++ script provided here, but this tutorial will be updated if an easier solution is found. (create a pull request or open an issue!)

```
ARCHITEUTHIS_FILTERED=${architeuthis_dir_path}/${samp_name}_filtered.output
ARCHITEUTHIS_FILTERED_REPORT=${architeuthis_dir_path}/${samp_name}_filtered_kraken.kreport
```

Run Architeuthis:

```
architeuthis mapping filter $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_FILTERED 

# Convert to Kraken report format
./kraken2-master/src/kraken2-report $DB_taxo $ARCHITEUTHIS_FILTERED $ARCHITEUTHIS_FILTERED_REPORT	
```

### 4. Estimating relative abundance

To estimate relative abundance with Bracken, you may use the default Bracken database provided with the SoilMicrobeDB Kraken2 database. It assumes your reads are 150bp; for different read lengths, a new Bracken database will need to be generated using developer instructions.

Install Bracken software using instructions in Step 0. Then run as follows:

```
BRACKEN_OUTPUT=${architeuthis_dir_path}/${samp_name}_filtered.b2

./Bracken/Bracken-master/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT
```

The above Bracken command can be modified to provide information at specific taxonomic ranks, such as species, genus, phylum, or domain. The output values represent the relative abundance of taxa among the classified reads from each sample.

### 5. Custom database updating

To create the database from a genome list, or update the existing SoilGenomeDB with a new set of genomes, we recommend using the Struo2 pipeline (Youngblut et al. 2021). This involves creating a Struo2 environment, and customizing input and output filepaths within the configuration file (YAML format). The genome input file is a TSV with columns labeling the genome, its taxonomy ID, and filepath. The file used to generate SoilMicrobeDB is available here:

[path]

To add new genomes, create a second TSV with the same format (genome name, taxonomy ID, and filepath). Genome names must be unique, and taxonomy IDs should be consistent with NCBI taxonomy structure. Then, specify the path to this new genome TSV file within the configuration file.

 Configuration settings (such as paths for the input TSV and database output directories) are specified in a YAML file that is provided to Snakemake. Additional genomes can be added to the database using a secondary TSV file and the "update" mode in the configuration YAML.

### 6. Cloud computing:

This database can also be integrated into cloud computing workflows using an Amazon Machine Image (AMI), offering benefits such as scalability for handling large datasets and cost-efficiency by optimizing resource usage during computational tasks. 

A detailed tutorial, provided by Robyn Wright, walks through the process of setting up an AMI for use with Kraken2. Follow the instructions, replacing the database download portion ("Get the Kraken2 database") with the following command: 

```
# Placeholder - replace with SoilGenomeDB once AWS URL is available
aws s3 cp --recursive --no-sign-request s3://kraken2-ncbi-refseq-complete-v205/Kraken2_RefSeqCompleteV205
```
Then, proceed with steps 0 and 2-4 of this tutorial to conduct analyses.
