#!/bin/bash -l

# Set SCC project
#$ -P dietzelab

# 512 GB node
#$ -pe omp 28
#$ -l mem_per_core=18G

#$ -j y # Merge stderr into the stdout file, to reduce clutter
#$ -o /projectnb/frpmars/soil_microbe_db/log_files/run_kraken_mars.log

# load kraken2

module load miniconda
conda activate struo2
module load kraken2

cd /projectnb/frpmars/soil_microbe_db/
	
NTHREADS=28

in_dir_path=/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_assemblies/metaGEM/workflow/dataset/
out_dir_path=/projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_classification/01_kraken_output
	

DBNAME=soil_microbe_db
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2

time_with_seconds=$(date +%T)
echo "Beginning Kraken2 loop at: $time_with_seconds"

for samp_dir in /projectnb/frpmars/soil_microbe_db/data/NEON_metagenome_assemblies/metaGEM/workflow/dataset/*; do
	samp_name="$(basename -- $samp_dir)"
	path_samp_prefix=${samp_dir}/${samp_name}
	READ1=${path_samp_prefix}_R1.fastq.gz
	READ2=${path_samp_prefix}_R2.fastq.gz
	echo "$samp_name kraken loop"

	REPORT_PATH=${out_dir_path}/${samp_name}_${DBNAME}_kraken.kreport
	OUTPUT_PATH=${out_dir_path}/${samp_name}_${DBNAME}_kraken.output

	if test -e $REPORT_PATH; then
		echo "$REPORT_PATH exists; skipping this run."
	else
		kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --threads $NTHREADS --report-zero-counts --memory-mapping --report-minimizer-data --paired $READ1 $READ2
	fi
done

time_with_seconds=$(date +%T)
echo "Finished Kraken2 loop at: $time_with_seconds"
