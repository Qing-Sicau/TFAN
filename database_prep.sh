# Choose a location with ample storage space

export DB_DIR="/home/chenq/fareannotate/annocript_v2/raw_database"
if [[ ! -e "$DB_DIR" ]]; then
    mkdir -p "$DB_DIR"
fi

# --- 1 swiss-prot download
# Navigate to your database directory

cd $DB_DIR

# --- Swiss-Prot (High-Quality, Curated) ---
# aria2c -c -d ./ "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
# diamond makedb --in uniprot_sprot.fasta.gz --db uniprot_sprot.dmnd

# --- TrEMBL (Comprehensive, Unreviewed) - NEW ---
# aria2c -c -d ./ "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz"
# diamond makedb --in uniprot_trembl.fasta.gz --db uniprot_trembl.dmnd

# --- UniRef90 (Clustered at 90% Identity) - NEW ---
# aria2c -c -d ./ "https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz"
# diamond makedb --in uniref90.fasta.gz --db uniref90.dmnd


# --- 2 interproscan data

cd $DB_DIR
mkdir interproscan && cd interproscan
# Download the complete offline data package (this is a very large file)
# The exact URL/filename changes with versions. Find the latest on the InterPro website.
# aria2c -c -d ./ "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.75-106.0/interproscan-5.75-106.0-64-bit.tar.gz"

# Unpack the data. This will create a 'data' subdirectory.
#tar -pxzf interproscan-*.tar.gz
#cd interproscan-*
#python3 setup.py -f interproscan.properties
#add the interproscan.sh to the env PATH



# InterProScan requires a properties file to be configured. The pipeline will handle this,
# but the 'data' directory must be present.

# --- 3 eggNOG

cd $DB_DIR
# You must have eggnog-mapper installed to download its data
# conda activate nextflow; conda install -c bioconda eggnog-mapper

# Download the data (e.g., for eukaryotes). This creates a 'eggnog_db' directory.
# This step requires an internet connection and is done only once.
# mkdir eggnog_db
# download_eggnog_data.py -y --data_dir eggnog_db
# gunzip eggnog_proteins.dmnd.gz
# gunzip -k eggnog.db.gz


# --- 4 RFAM

cd $DB_DIR

# Download the covariance model file
#aria2c -c -d ./ "https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz"

# Uncompress and then create the press index
# gunzip Rfam.cm.gz
# cmpress Rfam.cm
