
module load miniconda
conda activate struo2
cd /projectnb/frpmars/soil_microbe_db/

architeuthis_dir_path=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output
kraken_dir_path=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/01_kraken_output
bracken_outdir_path=/projectnb/dietzelab/zrwerbin/soil_genome_db_output/bracken_phylum

# Change to database before running loop
DBNAME=soil_microbe_db
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxonomy/
DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxo.k2d

#samp_file=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/01_kraken_output/ABBY_001-M-20170607-COMP_soil_microbe_db_kraken.kreport
time_with_seconds=$(date +%T)
echo "Beginning Architeuthis/Bracken loop at: $time_with_seconds"

for samp_file in /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/01_kraken_output/*.kreport; do

samp_file_basename="$(basename -- $samp_file)"
# remove ".kreport" (last 8 characters)
#samp_name=${samp_file_basename::-8}
samp_name=${samp_file_basename::-15}
echo "$samp_name bracken loop"

# Only run if file has the correct dbname
if [[ $samp_name != *$DBNAME* ]]; then
continue
else
    echo "cont"

ARCHITEUTHIS_FILTERED_REPORT=${architeuthis_dir_path}/${samp_name}_filtered_kraken.kreport
BRACKEN_OUTPUT_PHYLUM=${bracken_outdir_path}/${samp_name}_phylum_filtered.b2


# only run if braken output does not exist
if test -e $BRACKEN_OUTPUT_PHYLUM; then
echo "$BRACKEN_OUTPUT exists; skipping this run."
else
/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_PHYLUM -l P
fi

fi
done

time_with_seconds=$(date +%T)
echo "Finished Architeuthis/Bracken loop at: $time_with_seconds"






merged_file=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_phylum_merged.csv

lineage_merged_file=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_phylum_merged_lineage.csv

architeuthis merge -o $merged_file /projectnb/dietzelab/zrwerbin/soil_genome_db_output/bracken_phylum/*_soil_microbe_db_phylum_filtered.b2

architeuthis lineage $merged_file  --db $DBDIR --data-dir $DB_taxonomy_dir --out $lineage_merged_file

