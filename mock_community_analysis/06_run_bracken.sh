



cd /projectnb/frpmars/soil_microbe_db/mock_community_analysis/
module load miniconda
conda activate struo2

##### RUN BRACKEN ON FILTERED KRAKEN2 REPORT #####

bracken_dir_path=/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/06_bracken_output/
architeuthis_kraken_dir_path=/projectnb/frpmars/soil_microbe_db/mock_community_analysis/data/05_architeuthis_kraken_output/


DBNAME=gtdb_207
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/gtdb_207_filtered/kraken2


DBNAME=soil_microbe_db
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/soil_microbe_db/kraken2


declare -a samp_list=("fungprop_5_readdepth_.1m" 
                "fungprop_5_readdepth_1m" "fungprop_5_readdepth_5m"
                "fungprop_5_readdepth_10m"
                )


DBNAME=gtdb_207_unfiltered
DBDIR=/projectnb/microbiome/dgolden/Struo2/custom_dbs/GTDB_release207/kraken2

DBNAME=pluspf
DBDIR=/projectnb/frpmars/soil_microbe_db/databases/pluspf

DBNAME=pluspfp8
DBDIR=/projectnb/microbiome/ref_db/taxonomy_db/pluspfp8_kraken2_db

#for samp_name in fungprop_5_readdepth_.1m fungprop_5_readdepth_1m fungprop_5_readdepth_5m fungprop_5_readdepth_10m fungprop_10_readdepth_.1m fungprop_10_readdepth_1m fungprop_10_readdepth_5m fungprop_20_readdepth_.1m fungprop_20_readdepth_1m fungprop_20_readdepth_5m;


proportion=fungprop_5
proportion=fungprop_10
proportion=fungprop_15
proportion=fungprop_20
proportion=fungprop_1
proportion=fungprop_3


for readdepth in .1M 1M .5M 5M;

do

samp_name=${proportion}_readdepth_${readdepth}
echo "$samp_name"


ARCHITEUTHIS_FILTERED_REPORT=${architeuthis_kraken_dir_path}${samp_name}_${DBNAME}_filtered_kraken.kreport
BRACKEN_OUTPUT=${bracken_dir_path}/${samp_name}_${DBNAME}_filtered.b2
BRACKEN_OUTPUT_GENUS=${bracken_dir_path}/${samp_name}_${DBNAME}_genus_filtered.b2
BRACKEN_OUTPUT_DOMAIN=${bracken_dir_path}/${samp_name}_${DBNAME}_domain_filtered.b2

/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT

/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_GENUS -l G

/projectnb2/talbot-lab-data/zrwerbin/soil_genome_db/Struo2/Bracken/bracken -r 150 -d $DBDIR -i $ARCHITEUTHIS_FILTERED_REPORT -o $BRACKEN_OUTPUT_DOMAIN -l D

done

