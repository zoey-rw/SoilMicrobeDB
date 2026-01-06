## User Guide: SoilMicrobeDB AWS Tutorial

SoilMicrobeDB is a Kraken2 database with extensive representation of high-quality genomes of soil organisms, including uncultured and fungal species. Key applications of SoilMicrobeDB include advancing metagenomic research in soil ecology, identifying microbial taxa in diverse soil environments, and supporting agricultural and environmental monitoring. For details on how this database was designed, and potential applications, please refer to our manuscript (Werbin et al., 2025).

The tutorial below is for locally downloading this database to classify samples or update the database with new genomes. See section 6 ("Cloud Computing") for an alternative approach using an Amazon Machine Image. This tutorial walks through the process of evaluating the relative abundance of species within a metagenomic soil sample sequenced by the National Ecological Observatory Network. This tutorial is divided into the following steps:

0. Install software, define paths and environment variables  
1. Download the SoilMicrobeDB  
2. Classify sample reads  
3. Filter classifications for false-positives  
4. Estimate relative abundances of taxa

**Resource requirements:** 
- **Storage:** Downloading the complete SoilMicrobeDB (Step 1) requires approximately **289 GB** of disk space. This includes all Kraken2 classification files (~289 GB) and Bracken abundance estimation files (~400 MB). The Bracken files are included in the full download and are small relative to the main database files.
- **RAM:** Classification (Step 2) benefits from substantial RAM (~300 GB for optimal performance), though the `--memory-mapping` flag allows operation with lower RAM by accessing the database from disk. Steps 3 and 4 are not RAM-intensive.
- **Note:** If you only need classification without abundance estimation, you can skip the Bracken files, but they represent a small fraction of the total size (~400 MB out of 289 GB).

The recommended approach to using SoilMicrobeDB is to leverage the Struo2 pipeline for database development (Youngblut et al. 2022), with full installation details available here. Struo2 allows for customizing the database through additional genomes, and the pipeline includes the Kraken and Bracken software used in Steps 2 and 4. Step 3's false-positive filtering software (Architeuthis) must be installed separately (outlined in below in section 0.1).

### Quick start

```bash
# 0) Create/activate environment (see §0.1 for full install details)
git clone --recurse-submodules https://github.com/leylabmpi/Struo2.git
cd Struo2
conda env create --name struo2_conda -f conda_env.yaml
conda activate struo2_conda
mamba install -c bioconda kraken2
conda install -c conda-forge -c bioconda architeuthis

# 0.2) Inputs/outputs
mkdir -p ./test_sample ./test_output/{01_kraken,02_architeuthis,03_bracken}

# Download test sample files from S3 (or use your own files)
aws s3 cp --no-sign-request s3://kraken2-soil-microbe-database/test_sample_R1.fastq.gz ./test_sample/
aws s3 cp --no-sign-request s3://kraken2-soil-microbe-database/test_sample_R2.fastq.gz ./test_sample/

samp_dir=./test_sample
samp_name="$(basename -- "$samp_dir")"
path_samp_prefix="${samp_dir}/${samp_name}"
READ1="${path_samp_prefix}_R1.fastq.gz"
READ2="${path_samp_prefix}_R2.fastq.gz"

kraken_dir_path=./test_output/01_kraken
architeuthis_dir_path=./test_output/02_architeuthis
bracken_dir_path=./test_output/03_bracken

# 1) Download database (see §1 for full details)
# Create local database directory
mkdir -p /databases/SoilMicrobeDB/kraken2

# Download Kraken2 database files from S3
aws s3 sync --no-sign-request s3://kraken2-soil-microbe-database/kraken/ /databases/SoilMicrobeDB/kraken2/

# Set database variables
DBNAME=SoilMicrobeDB
DBDIR=/databases/SoilMicrobeDB/kraken2
DB_taxonomy_dir=/databases/SoilMicrobeDB/kraken2/taxonomy/
DB_taxo=$DBDIR/taxo.k2d

# 2) Kraken2 classification
NTHREADS=28
KRAKEN_REPORT="${kraken_dir_path}/${samp_name}_${DBNAME}_kraken.kreport"
KRAKEN_OUTPUT="${kraken_dir_path}/${samp_name}_${DBNAME}_kraken.output"

kraken2 --db "$DBDIR" --report "$KRAKEN_REPORT" --output "$KRAKEN_OUTPUT" \
  --threads "$NTHREADS" --report-zero-counts --memory-mapping --report-minimizer-data \
  --paired "$READ1" "$READ2"

# 3) False-positive filtering (Architeuthis)
ARCHITEUTHIS_FILTERED="${architeuthis_dir_path}/${samp_name}_filtered.output"
ARCHITEUTHIS_FILTERED_REPORT="${architeuthis_dir_path}/${samp_name}_filtered_kraken.kreport"

architeuthis mapping filter "$KRAKEN_OUTPUT" --db "$DBDIR" --data-dir "$DB_taxonomy_dir" \
  --out "$ARCHITEUTHIS_FILTERED"

./kraken2-master/src/kraken2-report "$DB_taxo" "$ARCHITEUTHIS_FILTERED" "$ARCHITEUTHIS_FILTERED_REPORT"

# 4) Bracken abundance estimates
BRACKEN_OUTPUT="${bracken_dir_path}/${samp_name}_abundance.tsv"
./Bracken/Bracken-master/bracken -r 150 -d "$DBDIR" -i "$ARCHITEUTHIS_FILTERED_REPORT" -o "$BRACKEN_OUTPUT"
```


### 0. Install software, define paths and environment variables

#### 0.1 Installation

To set up the Struo2 environment, run the following commands:

```bash
git clone --recurse-submodules https://github.com/leylabmpi/Struo2.git
cd Struo2
conda env create --name struo2_conda -f conda_env.yaml
conda activate struo2_conda
mamba install -c bioconda kraken2
```

Install quality-filter software (Architeuthis):

```bash
conda activate struo2_conda
conda install -c conda-forge -c bioconda architeuthis
```

(Optional) build `kraken2-report` from the Kraken2 source (needed in Step 3 if `kraken2-report` is not already available on your PATH):

```bash
wget https://github.com/DerrickWood/kraken2/archive/refs/heads/master.zip
unzip master.zip
cd kraken2-master/src
g++ -O3 -std=c++11 mmap_file.cc reports.cc taxonomy.cc kraken2-report.cpp -o kraken2-report
cd ../../
```

Install Bracken  
Additional details on Bracken software are available from the Bracken repository.

```bash
wget https://github.com/jenniferlu717/Bracken/archive/refs/heads/master.zip
unzip master.zip -d ./Bracken

# Sanity check what the zip unpacked into (usually ./Bracken/Bracken-master/)
ls ./Bracken

# Typical install path (adjust if your directory name differs)
bash ./Bracken/Bracken-master/install_bracken.sh

# If the zip unpacks directly into ./Bracken (no subfolder), use instead:
# bash ./Bracken/install_bracken.sh
```

#### 0.2 Specify input and output data paths

To ensure setup is completed smoothly, we provide a single "test" sample (forward/R1 and reverse/R2 reads) in the S3 bucket. Download the test samples and assign them to variables:

```bash
# Create test sample directory
mkdir -p ./test_sample

# Download test sample files from S3
aws s3 cp --no-sign-request s3://kraken2-soil-microbe-database/test_sample_R1.fastq.gz ./test_sample/
aws s3 cp --no-sign-request s3://kraken2-soil-microbe-database/test_sample_R2.fastq.gz ./test_sample/

# Confirm the expected files exist before proceeding
ls -lh ./test_sample

samp_dir=./test_sample
samp_name="$(basename -- "$samp_dir")"
path_samp_prefix="${samp_dir}/${samp_name}"
READ1="${path_samp_prefix}_R1.fastq.gz"
READ2="${path_samp_prefix}_R2.fastq.gz"

# Create output directories
mkdir -p ./test_output/{01_kraken,02_architeuthis,03_bracken}

kraken_dir_path=./test_output/01_kraken
architeuthis_dir_path=./test_output/02_architeuthis
bracken_dir_path=./test_output/03_bracken
```

### 1. Download the SoilMicrobeDB

Download the SoilMicrobeDB from the AWS S3 bucket using the AWS command line toolkit. The database files are organized in the `kraken/` subdirectory of the S3 bucket.

First, list the contents of the S3 bucket to verify access:

```bash
aws s3 ls --no-sign-request s3://kraken2-soil-microbe-database/
```

Create a local directory for the database and download the Kraken2 files:

```bash
# Create local database directory (adjust path as needed for your system)
mkdir -p /databases/SoilMicrobeDB/kraken2

# Download all Kraken2 database files from the kraken/ subdirectory
aws s3 sync --no-sign-request s3://kraken2-soil-microbe-database/kraken/ /databases/SoilMicrobeDB/kraken2/
```

**Using an external drive:** The database (~289 GB) can be stored on an external drive (e.g., a 1TB USB drive). This is useful if you don't have sufficient space on your main drive. To use an external drive:

1. Mount or connect your external drive and note its mount point (e.g., `/Volumes/MyDrive` on macOS or `/media/username/drive` on Linux)
2. Create the database directory on the external drive:
   ```bash
   # Example for macOS external drive
   mkdir -p /Volumes/MyDrive/SoilMicrobeDB/kraken2
   
   # Example for Linux external drive
   mkdir -p /media/username/drive/SoilMicrobeDB/kraken2
   ```
3. Download to the external drive location:
   ```bash
   aws s3 sync --no-sign-request s3://kraken2-soil-microbe-database/kraken/ /Volumes/MyDrive/SoilMicrobeDB/kraken2/
   ```
4. Update the `DBDIR` variable to point to the external drive location
5. **Important:** Use the `--memory-mapping` flag when running Kraken2 (already included in the tutorial commands) to access the database from disk rather than loading it entirely into RAM. This is especially important for external drives.

**Performance considerations:** External drives (especially USB 3.0+ or Thunderbolt) work well with Kraken2's memory-mapping mode. Classification may be slower than with an internal SSD, but the workflow will function correctly. Ensure the drive is properly mounted and has sufficient free space (~289 GB).

#### 1.1 Transferring database from remote server to external drive

If the database is already on a remote server and you want to copy it to an external drive for local testing, here are the fastest methods:

**Method 1: rsync (Recommended - resumable and efficient)**

`rsync` is the best choice for large file transfers because it:
- Can resume interrupted transfers
- Only transfers changed files (useful for updates)
- Shows progress
- Verifies file integrity

```bash
# Mount your external drive first, then:
# Replace with your actual paths
REMOTE_SERVER="user@server.edu"
REMOTE_DB_PATH="/path/to/remote/database/kraken2"
EXTERNAL_DRIVE="/Volumes/MyDrive/SoilMicrobeDB/kraken2"  # macOS
# EXTERNAL_DRIVE="/media/username/drive/SoilMicrobeDB/kraken2"  # Linux

# Create destination directory
mkdir -p "$EXTERNAL_DRIVE"

# Transfer with rsync (shows progress, preserves permissions)
rsync -avz --progress "${REMOTE_SERVER}:${REMOTE_DB_PATH}/" "$EXTERNAL_DRIVE/"

# Verify transfer completed successfully
rsync -avz --dry-run "${REMOTE_SERVER}:${REMOTE_DB_PATH}/" "$EXTERNAL_DRIVE/"
```

**Method 2: Direct copy from SSHFS mount (if already mounted)**

If you've already mounted the remote server using SSHFS (see `mount_remote.sh` in the repository):

```bash
# Assuming remote is mounted at ~/remote_soil_microbe_db
REMOTE_DB="/Users/yourname/remote_soil_microbe_db/databases/soil_microbe_db/kraken2"
EXTERNAL_DB="/Volumes/MyDrive/SoilMicrobeDB/kraken2"

# Use rsync for local-to-local copy (faster than cp for large files)
rsync -av --progress "$REMOTE_DB/" "$EXTERNAL_DB/"
```

**Method 3: scp (simpler but less efficient)**

For a one-time transfer, `scp` works but is less efficient than `rsync`:

```bash
# Transfer entire directory
scp -r user@server:/path/to/database/kraken2 /Volumes/MyDrive/SoilMicrobeDB/
```

**Transfer speed tips:**
- Use a wired network connection if possible (faster than WiFi)
- Transfer during off-peak hours if on a shared network
- Consider using `screen` or `tmux` to keep the transfer running if you disconnect
- Monitor progress: `rsync` shows transfer speed and ETA
- Verify after transfer: Check that all essential files are present (see file list below)

**Estimated transfer time:**
- Gigabit Ethernet (1 Gbps): ~40-60 minutes for 289 GB
- 100 Mbps connection: ~6-8 hours
- WiFi (varies): 8-24+ hours depending on signal strength

After download or transfer, set the following variables to point to your database location:

```bash
DBNAME=SoilMicrobeDB
DBDIR=/databases/SoilMicrobeDB/kraken2
DB_taxonomy_dir=/databases/SoilMicrobeDB/kraken2/taxonomy/
DB_taxo=$DBDIR/taxo.k2d
```

Verify the database structure. The `$DBDIR` directory should contain:

**Essential files for classification (Step 2):**
- `taxo.k2d` - taxonomy index file
- `hash.k2d` - hash table file (largest file, ~280+ GB)
- `opts.k2d` - options file
- `seqid2taxid.map` - sequence ID to taxonomy ID mapping
- `taxonomy/` - directory containing taxonomy files (names.dmp, nodes.dmp, etc.)

**Additional files for abundance estimation (Step 4):**
- `database150mers.kraken` - Bracken database file for 150bp reads
- `database150mers.kmer_distrib` - Bracken k-mer distribution file
- `database.kraken` - Bracken database file

**Outputs**
- A Kraken2-formatted database directory at `$DBDIR` containing the `.k2d` files, `taxonomy/` subdirectory, and Bracken database files.
- You should be able to run `kraken2 --db "$DBDIR" --help` without errors related to missing DB files.


### 2. Classify sample reads

To classify a sample, Kraken2 must be installed (Step 0). Modify the value of `NTHREADS` to reflect how many threads should run in parallel on your system, and specify variables to reflect the input and output locations.

#### 2.1 Evaluating RAM requirements

Before running classification, evaluate whether your system has sufficient RAM:

**Check your current RAM:**
```bash
# On macOS
sysctl hw.memsize
# Or use Activity Monitor (GUI)

# On Linux
free -h
# Or
cat /proc/meminfo | grep MemTotal
```

**Understanding RAM requirements:**
- **With `--memory-mapping` (recommended for laptops):** Kraken2 accesses the database from disk rather than loading it entirely into RAM. This allows classification on systems with **8-16 GB RAM or more**, though performance will be slower than with more RAM. The `--memory-mapping` flag is already included in the tutorial commands.
- **Without `--memory-mapping`:** The database would need to be loaded into RAM, requiring ~300 GB RAM for optimal performance. This is typically only available on high-performance computing clusters.

**Testing RAM usage with a small sample:**
1. Start with the test sample provided in the tutorial (smaller than full datasets)
2. Monitor RAM usage during classification:
   ```bash
   # On macOS - monitor in another terminal
   top -l 1 | grep "PhysMem"
   
   # On Linux - monitor in another terminal
   watch -n 1 free -h
   ```
3. If classification completes successfully without running out of memory, your system should handle larger samples

**If you encounter memory issues:**
- Ensure `--memory-mapping` is included in your command (it's in the tutorial)
- Reduce the number of threads (`NTHREADS`) to lower memory usage
- Process smaller batches of reads if working with very large datasets
- Close other memory-intensive applications during classification

```
NTHREADS=28
KRAKEN_REPORT=${kraken_dir_path}/${samp_name}_${DBNAME}_kraken.kreport
KRAKEN_OUTPUT=${kraken_dir_path}/${samp_name}_${DBNAME}_kraken.output
```

Then, classify a sample with the following command:

```bash
kraken2 --db $DBDIR --report $KRAKEN_REPORT --output $KRAKEN_OUTPUT --threads $NTHREADS \
  --report-zero-counts --memory-mapping --report-minimizer-data --paired $READ1 $READ2
```

**Note on `--memory-mapping`:** This flag uses memory-mapped file I/O, which allows Kraken2 to access the database from disk without loading it entirely into RAM. This is essential for systems with limited RAM (like most laptops). With `--memory-mapping`, classification can run on systems with 8-16 GB RAM, though it will be slower than systems with more RAM. If your system has very large amounts of RAM (300GB+), you can remove `--memory-mapping` to load the database directly into memory for faster performance, but this is typically only available on high-performance computing clusters.

**Outputs**
- `$KRAKEN_OUTPUT`: per-read classifications (Kraken2 output format)
- `$KRAKEN_REPORT`: summary report (Kraken report / `.kreport`)


### 3. Filter classifications for false-positives

To run the false-positive filtering step, we recommend the software Architeuthis (Diener et al. 2024). Install using the instructions provided by Architeuthis. The output files provide per-read results, but must be converted to the Kraken2 report format before passing to Bracken (Step 4). Currently, the only way to do this is by compiling a C++ script provided here, but this tutorial will be updated if an easier solution is found. (Create a pull request or open an issue!)

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

**Outputs**
- `$ARCHITEUTHIS_FILTERED`: filtered classifications (Kraken2 output format)
- `$ARCHITEUTHIS_FILTERED_REPORT`: filtered report in `.kreport` format


### 4. Estimate relative abundances of taxa

To estimate relative abundance with Bracken, you may use the default Bracken database provided with the SoilMicrobeDB Kraken2 database. It assumes your reads are 150bp; for different read lengths, a new Bracken database will need to be generated using developer instructions.

Install Bracken software using instructions in Step 0. Then run as follows:

```
BRACKEN_OUTPUT=${bracken_dir_path}/${samp_name}_abundance.tsv

./Bracken/Bracken-master/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT
```

> If your Bracken install produced a different folder name or binary location, update the path above accordingly (e.g., `./Bracken/bracken`).


The above Bracken command can be modified to provide information at specific taxonomic ranks, such as species, genus, phylum, or domain. The output values represent the relative abundance of taxa among the classified reads from each sample.

**Outputs**
- `$BRACKEN_OUTPUT`: abundance estimates in tab-delimited format (TSV)


### 5. Custom database updating

To create the database from a genome list, or update the existing SoilMicrobeDB with a new set of genomes, we recommend using the Struo2 pipeline (Youngblut et al. 2021). This involves creating a Struo2 environment, and customizing input and output filepaths within the configuration file (YAML format). The genome input file is a TSV with columns labeling the genome, its taxonomy ID, and filepath. The file used to generate SoilMicrobeDB is available in the S3 bucket:

```bash
# Download the genome table file
aws s3 cp --no-sign-request s3://kraken2-soil-microbe-database/soil_microbe_db_genome_table.csv ./
```

To add new genomes, create a second TSV with the same format (genome name, taxonomy ID, and filepath). Genome names must be unique, and taxonomy IDs should be consistent with NCBI taxonomy structure. Then, specify the path to this new genome TSV file within the configuration file.

Configuration settings (such as paths for the input TSV and database output directories) are specified in a YAML file that is provided to Snakemake. Additional genomes can be added to the database using a secondary TSV file and the "update" mode in the configuration YAML.

### 6. Cloud computing

This database can also be integrated into cloud computing workflows using an Amazon Machine Image (AMI), offering benefits such as scalability for handling large datasets and cost-efficiency by optimizing resource usage during computational tasks.

A detailed tutorial, provided by Robyn Wright, walks through the process of setting up an AMI for use with Kraken2. Follow the instructions, replacing the database download portion ("Get the Kraken2 database") with the following command:

```bash
# Download the Kraken2 database files from S3
aws s3 sync --no-sign-request s3://kraken2-soil-microbe-database/kraken/ /path/to/local/kraken2/
```

Then, proceed with steps 0 and 2-4 of this tutorial to conduct analyses.

### 7. Extracting unclassified reads for secondary classification

In soil metagenomics, a portion of reads often remains unclassified due to extreme diversity or novel taxa. To maximize data recovery, you can extract these unclassified reads and pass them to a secondary, more generalized database to identify non-soil-specific microbes or contaminants.

#### 7.1 Extracting unclassified reads

To capture reads that did not match SoilMicrobeDB, you must re-run or modify your Step 2 Kraken2 command to include the `--unclassified-out` flag.

```
# Define output path for unclassified reads
UNCLASS_OUT=${kraken_dir_path}/${samp_name}_unclassified#.fastq

# Run Kraken2 with extraction flag
kraken2 --db $DBDIR   --report $KRAKEN_REPORT   --output $KRAKEN_OUTPUT   --threads $NTHREADS   --memory-mapping   --unclassified-out $UNCLASS_OUT   --paired $READ1 $READ2
```

Note: The `#` character in the filename is a wildcard. Kraken2 will automatically generate `_1.fastq` and `_2.fastq` for your paired-end data.

#### 7.2 Secondary database selection

Choosing the right secondary database depends on your research goals and available computational resources.

| Database | Contents | Disk/RAM (2026) |
|---|---|---|
| Standard | Archaea, Bacteria, Viral, Human, UniVec | ~75–80 GB |
| PlusPF | Standard + Protozoa and Fungi | ~140–160 GB |
| PlusPFP | Standard + Protozoa, Fungi, and Plants | ~190–210 GB |
| Core NT | GenBank, RefSeq, TPA, and PDB (Inclusive) | ~450–500 GB |

> **Note:** These disk/RAM figures are rough planning numbers and can change as Kraken2/Bracken databases evolve. Always check the published size of the specific DB build you are using before provisioning compute.


#### 7.3 Classification against a secondary database

Once you have extracted the unclassified reads, you can "subtract" the soil-specific knowledge and query the broader biological space.

```
# Define paths for secondary database
SEC_DB_PATH=/path/to/secondary/database
SEC_REPORT=${kraken_dir_path}/${samp_name}_secondary_refseq.kreport

# Run classification on the unclassified subset
kraken2 --db $SEC_DB_PATH   --threads $NTHREADS   --report $SEC_REPORT   --paired ${kraken_dir_path}/${samp_name}_unclassified_1.fastq          ${kraken_dir_path}/${samp_name}_unclassified_2.fastq
```

#### 7.4 Merging and interpreting results

After secondary classification, you will have two reports:

- **SoilMicrobeDB report:** high-resolution soil microbial taxa.
- **Secondary report:** broad-spectrum taxa (e.g., human/host contamination or rare environmental microbes).

To create a final abundance table that combines your SoilMicrobeDB results with your Secondary Database results, use `combine_kreports.py` from KrakenTools. This script merges multiple `.kreport` files into a single matrix.

##### 7.4.1 Install KrakenTools

If you haven't already, clone the repository:

```
git clone https://github.com/jenniferlu717/KrakenTools.git
```

##### 7.4.2 Run the merge command

This command will take your soil-specific report and your secondary "unclassified-subset" report and align them by taxonomic ID.

```
# Define the reports to merge
SOIL_REPORT=${kraken_dir_path}/${samp_name}_SoilMicrobeDB_kraken.kreport
SEC_REPORT=${kraken_dir_path}/${samp_name}_secondary_refseq.kreport
FINAL_COMBINED=${kraken_dir_path}/${samp_name}_combined_taxonomy.txt

# Run combine_kreports.py
python3 ./KrakenTools/combine_kreports.py   --report-files $SOIL_REPORT $SEC_REPORT   --output $FINAL_COMBINED   --display-headers
```

##### 7.4.3 Understanding the output

The combined file will contain columns for each report:

- Column 1–n: percentage and read counts for each input file.
- Final columns: aggregated sums across all samples/databases.

If you notice duplicate taxa appearing in both reports (e.g., a common bacterium in both SoilMicrobeDB and RefSeq), the combined report will show two distinct count columns for that TaxID. This is useful for identifying "leakage" where a specific soil organism was missed by the first database due to k-mer specificity but caught by the broader one.

**Outputs (Section 7)**
- Unclassified reads FASTQs (from `--unclassified-out`, typically `*_R1.fastq` and `*_R2.fastq` or gzipped equivalents)
- A secondary `.kreport` and `.output` generated by classifying the unclassified subset against the secondary DB
- A merged `.kreport` produced by `combine_kreports.py` (primary + secondary)

