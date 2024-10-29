

cd /projectnb/frpmars/soil_microbe_db/mock_community_analysis/
module load miniconda
conda activate struo2

##### RUN FILTER ON KRAKEN2 OUTPUT #####
kraken_dir_path=/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/03_kraken_output/
architeuthis_dir_path=/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/04_architeuthis_output/
architeuthis_kraken_dir_path=/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/05_architeuthis_kraken_output/
	
	
samp_name=fungprop_5_readdepth_5m
samp_name=fungprop_5_readdepth_1m
samp_name=fungprop_5_readdepth_.1m
samp_name=fungprop_5_readdepth_10m

samp_name=fungprop_10_readdepth_.1m
samp_name=fungprop_10_readdepth_1m
samp_name=fungprop_10_readdepth_5m

samp_name=fungprop_20_readdepth_5m
samp_name=fungprop_20_readdepth_.1m
samp_name=fungprop_20_readdepth_1m
#samp_name=fungprop_20_readdepth_5m

DBNAME=gtdb_207
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxonomy/
DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxo.k2d


DBNAME=gtdb_207_unfiltered
#DBDIR=/projectnb/microbiome/ref_db/GTDB_207_kraken2
#DB_taxonomy_dir=/projectnb/microbiome/ref_db/GTDB_207_kraken2/taxonomy
DBDIR=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2
DB_taxonomy_dir=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/taxonomy
DB_taxo=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2/taxo.k2d



DBNAME=soil_microbe_db
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxonomy/
DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxo.k2d


DBNAME=pluspf
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/pluspf
DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/pluspf/taxo.k2d
#DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/pluspf/taxonomy
DB_taxonomy_dir=/projectnb/microbiome/ref_db/NCBI-taxdump


DBNAME=pluspfp8
DBDIR=/projectnb/microbiome/ref_db/taxonomy_db/pluspfp8_kraken2_db
#DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxonomy/
DB_taxonomy_dir=/projectnb/microbiome/ref_db/NCBI-taxdump
DB_taxo=/projectnb/microbiome/ref_db/taxonomy_db/pluspfp8_kraken2_db/taxo.k2d


# DB_taxo=/projectnb/frpmars/soil_microbe_db/databases/pluspf_kraken2_db/taxo.k2d
# DB_taxonomy_dir=/projectnb/microbiome/ref_db/taxonomy_db/pluspf_kraken2_db


samp_name=fungprop_5_readdepth_1m

for samp_name in fungprop_5_readdepth_1m fungprop_20_readdepth_.1m fungprop_20_readdepth_1m fungprop_20_readdepth_5m fungprop_5_readdepth_.1m fungprop_10_readdepth_.1m fungprop_10_readdepth_1m fungprop_10_readdepth_5m fungprop_5_readdepth_5m fungprop_5_readdepth_10m;

do
echo "$samp_name"


KRAKEN_OUTPUT=${kraken_dir_path}${samp_name}_${DBNAME}_kraken.output
ARCHITEUTHIS_FILTERED=${architeuthis_dir_path}${samp_name}_${DBNAME}_filtered.output
ARCHITEUTHIS_SUMMARY=${architeuthis_dir_path}${samp_name}_${DBNAME}_summary.output
ARCHITEUTHIS_FILTERED_REPORT=${architeuthis_kraken_dir_path}${samp_name}_${DBNAME}_filtered_kraken.kreport
ARCHITEUTHIS_SCORES=${architeuthis_dir_path}${samp_name}_${DBNAME}_scores.output

# Run filter
architeuthis mapping filter $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_FILTERED 
architeuthis mapping summary $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_SUMMARY
architeuthis mapping score $KRAKEN_OUTPUT --db $DBDIR --data-dir $DB_taxonomy_dir --out $ARCHITEUTHIS_SCORES

# Convert back to Kraken output
cd /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/misc_scripts/mock_community/kraken2_report/daydream-boost-kraken2-e533c50/src
#g++ -O3 -std=c++11 mmap_file.cc reports.cc taxonomy.cc kraken2-report.cpp -o kraken2-report
./kraken2-report $DB_taxo $ARCHITEUTHIS_FILTERED $ARCHITEUTHIS_FILTERED_REPORT

done

cd /projectnb/microbiome/ref_db/NCBI-taxdump/.snapshots


cd /projectnb/microbiome/ref_db/taxonomy_db/pluspf_kraken2_db/.snapshots

cd /projectnb/microbiome/ref_db/taxonomy_db/pluspf_kraken2_db/.snapshots

cd /projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/.snapshots

cd /projectnb/microbiome/ref_db/GTDB_207_kraken2/.snapshots
