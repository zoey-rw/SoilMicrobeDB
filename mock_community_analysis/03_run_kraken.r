
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

proportion=fungprop_1
proportion=fungprop_3
proportion=fungprop_20



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

BRACKEN_OUTPUT=${out_dir_path}${samp_name}_${DBNAME}_filtered.b2
BRACKEN_OUTPUT_GENUS=${out_dir_path}${samp_name}_${DBNAME}_genus_filtered.b2
BRACKEN_OUTPUT_DOMAIN=${out_dir_path}${samp_name}_${DBNAME}_domain_filtered.b2
BRACKEN_OUTPUT_PHYLUM=${out_dir_path}${samp_name}_${DBNAME}_phylum_filtered.b2
BRACKEN_OUTPUT_CLASS=${out_dir_path}${samp_name}_${DBNAME}_class_filtered.b2
BRACKEN_OUTPUT_ORDER=${out_dir_path}${samp_name}_${DBNAME}_order_filtered.b2
BRACKEN_OUTPUT_FAMILY=${out_dir_path}${samp_name}_${DBNAME}_family_filtered.b2

if test -e $BRACKEN_OUTPUT; then
    echo "$BRACKEN_OUTPUT exists; skipping species-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT
fi

if test -e $BRACKEN_OUTPUT_GENUS; then
    echo "$BRACKEN_OUTPUT_GENUS exists; skipping genus-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT_GENUS -l G
fi

if test -e $BRACKEN_OUTPUT_DOMAIN; then
    echo "$BRACKEN_OUTPUT_DOMAIN exists; skipping domain-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT_DOMAIN -l D
fi

if test -e $BRACKEN_OUTPUT_PHYLUM; then
    echo "$BRACKEN_OUTPUT_PHYLUM exists; skipping phylum-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT_PHYLUM -l P
fi

if test -e $BRACKEN_OUTPUT_CLASS; then
    echo "$BRACKEN_OUTPUT_CLASS exists; skipping class-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT_CLASS -l C
fi

if test -e $BRACKEN_OUTPUT_ORDER; then
    echo "$BRACKEN_OUTPUT_ORDER exists; skipping order-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT_ORDER -l O
fi

if test -e $BRACKEN_OUTPUT_FAMILY; then
    echo "$BRACKEN_OUTPUT_FAMILY exists; skipping family-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT_FAMILY -l F
fi
	
done






# for samp_name in fungprop_20_readdepth_.1m fungprop_20_readdepth_1m fungprop_20_readdepth_5m fungprop_5_readdepth_.1m fungprop_10_readdepth_.1m fungprop_10_readdepth_1m fungprop_10_readdepth_5m fungprop_5_readdepth_5m fungprop_5_readdepth_10m;

do
echo "$samp_name"

READ1=${in_dir_path}${samp_name}_R1.fastq
READ2=${in_dir_path}${samp_name}_R2.fastq

REPORT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.kreport
OUTPUT_PATH=${out_dir_path}${samp_name}_${DBNAME}_kraken.output

kraken2 --db $DBDIR --report $REPORT_PATH --output $OUTPUT_PATH --threads $NTHREADS --report-zero-counts --memory-mapping --report-minimizer-data --paired $READ1 $READ2

BRACKEN_OUTPUT=${out_dir_path}${samp_name}_${DBNAME}_filtered.b2
BRACKEN_OUTPUT_GENUS=${out_dir_path}${samp_name}_${DBNAME}_genus_filtered.b2
BRACKEN_OUTPUT_DOMAIN=${out_dir_path}${samp_name}_${DBNAME}_domain_filtered.b2
BRACKEN_OUTPUT_PHYLUM=${out_dir_path}${samp_name}_${DBNAME}_phylum_filtered.b2
BRACKEN_OUTPUT_CLASS=${out_dir_path}${samp_name}_${DBNAME}_class_filtered.b2
BRACKEN_OUTPUT_ORDER=${out_dir_path}${samp_name}_${DBNAME}_order_filtered.b2
BRACKEN_OUTPUT_FAMILY=${out_dir_path}${samp_name}_${DBNAME}_family_filtered.b2

if test -e $BRACKEN_OUTPUT; then
    echo "$BRACKEN_OUTPUT exists; skipping species-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT
fi

if test -e $BRACKEN_OUTPUT_GENUS; then
    echo "$BRACKEN_OUTPUT_GENUS exists; skipping genus-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT_GENUS -l G
fi

if test -e $BRACKEN_OUTPUT_DOMAIN; then
    echo "$BRACKEN_OUTPUT_DOMAIN exists; skipping domain-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT_DOMAIN -l D
fi

if test -e $BRACKEN_OUTPUT_PHYLUM; then
    echo "$BRACKEN_OUTPUT_PHYLUM exists; skipping phylum-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT_PHYLUM -l P
fi

if test -e $BRACKEN_OUTPUT_CLASS; then
    echo "$BRACKEN_OUTPUT_CLASS exists; skipping class-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT_CLASS -l C
fi

if test -e $BRACKEN_OUTPUT_ORDER; then
    echo "$BRACKEN_OUTPUT_ORDER exists; skipping order-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT_ORDER -l O
fi

if test -e $BRACKEN_OUTPUT_FAMILY; then
    echo "$BRACKEN_OUTPUT_FAMILY exists; skipping family-level Bracken."
else
    /projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $REPORT_PATH -o $BRACKEN_OUTPUT_FAMILY -l F
fi

done
	
	
#DB_taxonomy_dir=/projectnb/microbiome/ref_db/NCBI-taxdump
