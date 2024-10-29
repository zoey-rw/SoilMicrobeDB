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
cd /projectnb/frpmars/soil_microbe_db/
architeuthis_dir_path=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output
kraken_dir_path=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/01_kraken_output


# Change to database before running loop
DBNAME=soil_microbe_db
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxonomy/
DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxo.k2d


DBNAME=pluspf
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/pluspf
DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/pluspf/taxo.k2d
DB_taxonomy_dir=/projectnb/microbiome/ref_db/NCBI-taxdump


DBNAME=gtdb_207
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxonomy/
DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxo.k2d


time_with_seconds=$(date +%T)
echo "Beginning Architeuthis/Bracken loop at: $time_with_seconds"

for samp_file in /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/01_kraken_output/*.kreport; do

samp_file_basename="$(basename -- $samp_file)"

samp_name=${samp_file_basename::-15}
echo "$samp_name architeuthis loop"

# Only run if file has the correct dbname
if [[ $samp_name != *$DBNAME* ]]; then
continue
else
echo "cont"

# not for existing soilmicrobedb output
#KRAKEN_OUTPUT=${kraken_dir_path}/${samp_name}_${DBNAME}_kraken.output
KRAKEN_OUTPUT=${kraken_dir_path}/${samp_name}_kraken.output


ARCHITEUTHIS_FILTERED=${architeuthis_dir_path}/${samp_name}_filtered.output
ARCHITEUTHIS_SUMMARY=${architeuthis_dir_path}/${samp_name}_summary.output
ARCHITEUTHIS_SCORES=${architeuthis_dir_path}/${samp_name}_scores.output
ARCHITEUTHIS_FILTERED_REPORT=${architeuthis_dir_path}/${samp_name}_filtered_kraken.kreport

BRACKEN_OUTPUT=${architeuthis_dir_path}/${samp_name}_filtered.b2
BRACKEN_OUTPUT_GENUS=${architeuthis_dir_path}/${samp_name}_genus_filtered.b2
BRACKEN_OUTPUT_DOMAIN=${architeuthis_dir_path}/${samp_name}_domain_filtered.b2
	
	# only run if filtering output does not exist
if test -e $ARCHITEUTHIS_FILTERED_REPORT; then
echo "$ARCHITEUTHIS_FILTERED_REPORT exists; skipping this run."

else

# Run filter
architeuthis mapping filter $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_FILTERED 
#architeuthis mapping summary $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_SUMMARY
#architeuthis mapping score $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_SCORES

# Convert back to Kraken output
cd /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/misc_scripts/mock_community/kraken2_report/daydream-boost-kraken2-e533c50/src
g++ -O3 -std=c++11 mmap_file.cc reports.cc taxonomy.cc kraken2-report.cpp -o kraken2-report
./kraken2-report $DB_taxo $ARCHITEUTHIS_FILTERED $ARCHITEUTHIS_FILTERED_REPORT	
fi

	
# only run if braken output does not exist
if test -e $BRACKEN_OUTPUT; then
echo "$BRACKEN_OUTPUT exists; skipping this run."
else
/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT
/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_GENUS -l G
/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_DOMAIN -l D
fi

fi
done

time_with_seconds=$(date +%T)
echo "Finished Architeuthis/Bracken loop at: $time_with_seconds"




# Everything below is miscellaneous/test code not integrated into workflow

# Merge all domain files to get % filtered 
architeuthis merge -o /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/soil_microbe_db_filtered_domain_merged.csv /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output/*_soil_microbe_db_domain_filtered.b2
architeuthis merge -o /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/pluspf_filtered_domain_merged.csv /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output/*_pluspf_domain_filtered.b2
architeuthis merge -o /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/gtdb_207_filtered_domain_merged.csv /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output/*_gtdb_207_domain_filtered.b2

# Merge all genus files to get % abundance 
architeuthis merge -o /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/soil_microbe_db_filtered_genus_merged.csv /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output/*_soil_microbe_db_genus_filtered.b2
architeuthis merge -o /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/pluspf_filtered_genus_merged.csv /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output/*_pluspf_genus_filtered.b2
architeuthis merge -o /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/gtdb_207_filtered_genus_merged.csv /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output/*_gtdb_207_genus_filtered.b2

# Compressing some of the largest Kraken output files to save space
module load pigz
cd /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/01_kraken_output

tar --remove-files --use-compress-program="pigz --best --recursive" -cf /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/01_kraken_output/CLBJ_soil_microbe_db_kraken_compressed.tar.gz /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/01_kraken_output/CLBJ*soil_microbe_db_kraken.output

# Checking for incompletely written files

# function from https://stackoverflow.com/questions/38746/how-to-detect-file-ends-in-newline
function file_ends_with_newline() {
    [[ $(tail -c1 "$1" | wc -l) -gt 0 ]]
}


for samp_file in /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output/*_filtered.output; do
if ! file_ends_with_newline $samp_file
then
    echo "$samp_file likely incomplete"
fi
done


for samp_file in /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output/*_scores.output; do
if ! file_ends_with_newline $samp_file
then
    echo "$samp_file likely incomplete"
fi
done


