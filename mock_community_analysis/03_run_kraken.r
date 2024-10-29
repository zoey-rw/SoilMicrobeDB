



# Step 3: Run Kraken2 on mock community
qrsh -pe omp 28 -l mem_per_core=16G

cd /projectnb/frpmars/soil_microbe_db/mock_community_analysis/
module load miniconda
conda activate struo2
module load kraken2/2.1.2

NTHREADS=28

	
out_dir_path=/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/03_kraken_output/
in_dir_path=/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/02_simulated_reads/

samp_name=fungprop_5_readdepth_1m
samp_name=fungprop_5_readdepth_5m
samp_name=fungprop_5_readdepth_10m

samp_name=fungprop_20_readdepth_1m
samp_name=fungprop_20_readdepth_.1m
samp_name=fungprop_20_readdepth_5m


samp_name=fungprop_10_readdepth_1m
samp_name=fungprop_10_readdepth_.1m
samp_name=fungprop_10_readdepth_5m

samp_name=fungprop_5_readdepth_.1m


READ1=${in_dir_path}${samp_name}_R1.fastq
READ2=${in_dir_path}${samp_name}_R2.fastq

DBNAME=gtdb_207
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2/taxonomy/

REPORT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.kreport
OUTPUT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.output
	
kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --threads $NTHREADS --report-zero-counts --memory-mapping --report-minimizer-data --paired $READ1 $READ2
	

DBNAME=soil_microbe_db
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2
DB_taxonomy_dir=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2/taxonomy/

REPORT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.kreport
OUTPUT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.output

kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --threads $NTHREADS --report-zero-counts --memory-mapping --report-minimizer-data --paired $READ1 $READ2


DBNAME=pluspf
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/pluspf
	
REPORT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.kreport
OUTPUT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.output

kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --threads $NTHREADS --report-zero-counts --memory-mapping --report-minimizer-data --paired $READ1 $READ2

	

DBNAME=pluspfp8
DBDIR=/projectnb/microbiome/ref_db/taxonomy_db/pluspfp8_kraken2_db


	
DBNAME=gtdb_207_unfiltered
#DBDIR=/projectnb/microbiome/ref_db/GTDB_207_kraken2
#DB_taxonomy_dir=/projectnb/microbiome/ref_db/GTDB_207_kraken2/taxonomy
DBDIR=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2
#DB_taxonomy_dir=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/taxonomy

fungprop_5_readdepth_1m


for samp_name in fungprop_20_readdepth_.1m fungprop_20_readdepth_1m fungprop_20_readdepth_5m;

for samp_name in fungprop_10_readdepth_.1m fungprop_10_readdepth_1m fungprop_10_readdepth_5m;

for samp_name in fungprop_20_readdepth_.1m fungprop_20_readdepth_1m fungprop_20_readdepth_5m fungprop_5_readdepth_.1m fungprop_10_readdepth_.1m fungprop_10_readdepth_1m fungprop_10_readdepth_5m fungprop_5_readdepth_5m fungprop_5_readdepth_10m;

do
echo "$samp_name"

READ1=${in_dir_path}${samp_name}_R1.fastq
READ2=${in_dir_path}${samp_name}_R2.fastq


REPORT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.kreport
OUTPUT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.output


kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --threads $NTHREADS --report-zero-counts --memory-mapping --report-minimizer-data --paired $READ1 $READ2

done
	
	
DB_taxonomy_dir=/projectnb/microbiome/ref_db/NCBI-taxdump
