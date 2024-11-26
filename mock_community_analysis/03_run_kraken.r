
# Step 3: Run Kraken2 on mock community
qrsh -pe omp 28 -l mem_per_core=16G

cd /projectnb/frpmars/soil_microbe_db/mock_community_analysis/
module load miniconda
conda activate struo2
module load kraken2/2.1.2

NTHREADS=28
	
out_dir_path=/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/03_kraken_output/
in_dir_path=/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/02_simulated_reads/

# Run for each database and each fung_prop

DBNAME=soil_microbe_db
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2

DBNAME=gtdb_207
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2

DBNAME=pluspf
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/pluspf

DBNAME=gtdb_207_unfiltered
DBDIR=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2

DBNAME=pluspfp8
DBDIR=/projectnb/microbiome/ref_db/taxonomy_db/pluspfp8_kraken2_db

proportion=fungprop_5
proportion=fungprop_10
proportion=fungprop_15




readdepth=5M
readdepth=1M
readdepth=.1M
readdepth=.5M
for readdepth in .1M .5M 1M 5M;

do

samp_name=${proportion}_readdepth_${readdepth}

echo $samp_name

READ1=${in_dir_path}${samp_name}_R1.fastq
READ2=${in_dir_path}${samp_name}_R2.fastq

REPORT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.kreport
OUTPUT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.output
	
kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --threads $NTHREADS --report-zero-counts --memory-mapping --report-minimizer-data --paired $READ1 $READ2
	
done






# for samp_name in fungprop_20_readdepth_.1m fungprop_20_readdepth_1m fungprop_20_readdepth_5m fungprop_5_readdepth_.1m fungprop_10_readdepth_.1m fungprop_10_readdepth_1m fungprop_10_readdepth_5m fungprop_5_readdepth_5m fungprop_5_readdepth_10m;

do
echo "$samp_name"

READ1=${in_dir_path}${samp_name}_R1.fastq
READ2=${in_dir_path}${samp_name}_R2.fastq

REPORT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.kreport
OUTPUT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.output

kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --threads $NTHREADS --report-zero-counts --memory-mapping --report-minimizer-data --paired $READ1 $READ2

done
	
	
#DB_taxonomy_dir=/projectnb/microbiome/ref_db/NCBI-taxdump
