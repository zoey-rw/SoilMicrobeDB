# After running Kraken, Architeuthis, and Bracken
# Merge Bracken outputs with "architeuthis merge"
# Then, lineage info can be added to the merged file with "architeuthis lineage"

module load miniconda
conda activate struo2
cd /projectnb/frpmars/soil_microbe_db/
    
    
DBNAME=soil_microbe_db
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxonomy/


# Merge at species-level for Soil Microbe DB
merged_file=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_species_merged.csv
lineage_file=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_species_merged_lineage.csv
architeuthis merge -o $merged_file /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output/*soil_microbe_db_filtered.b2
architeuthis lineage $merged_file  --db $DBDIR --data-dir $DB_taxonomy_dir --out $lineage_file

# genus level?
#architeuthis merge -o /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_genus_merged_lineage.csv /projectnb/dietzelab/zrwerbin/soil_genome_db_output/bracken_with_lineage/*_soil_microbe_db_genus_lineage.b2
#merged_file=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/soil_microbe_db_filtered_genus_merged.csv
#architeuthis lineage $merged_file  --db $DBDIR --data-dir $DB_taxonomy_dir --out /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/soil_microbe_db_filtered_genus_merged_lineage.csv






DBNAME=gtdb_207
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxonomy/

# Now merge at species-level for gtdb
merged_file=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/gtdb_species_merged.csv
lineage_file=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/gtdb_species_merged_lineage.csv
architeuthis merge -o $merged_file /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output/*gtdb_207_filtered.b2
architeuthis lineage $merged_file  --db $DBDIR --data-dir $DB_taxonomy_dir --out $lineage_file



DBNAME=pluspf
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/pluspf
DB_taxonomy_dir=/projectnb/microbiome/ref_db/NCBI-taxdump

# Now merge at species-level for pluspf
merged_file=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/pluspf_species_merged.csv
lineage_file=/projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/summary_files/pluspf_species_merged_lineage.csv
architeuthis merge -o $merged_file /projectnb/frpmars/soil_microbe_db/NEON_metagenome_classification/02_bracken_output/*pluspf_filtered.b2
architeuthis lineage $merged_file  --db $DBDIR --data-dir $DB_taxonomy_dir --out $lineage_file




