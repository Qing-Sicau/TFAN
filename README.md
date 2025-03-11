# TFAN
Transcriptome sequence functional annotation

# Description
This script was used to annotate transcript sequences in fasta format, either after genome annotation or assembly from RNA-seq reads

It will iteratively blastx the transcripts against swiss_prot, swiss_trembl, refseq90, pfam domains, and retrive the corresponding descriptions as well as th go terms corresponding to the best hit.

## How to use

- first, download corresponding data into a directory using either your favrate tool;
- put them into a directory;
- install necessary tools using conda, or your way;
- run `` python annotator.py `` once, this will generate a default config file in yaml;
- run again using parameters, -h might help with which could be modified;

## Logic

  it will use diamond to first blastx the query to the swiss-prot file, get the hit; left fasta transcript will subject to blastx to swiss-trembl, get the hit; and then left ones to uniref90; then all transcripts will be translated deduced by transdecoder, put the proteins into pfam search to scan the domains.

  Finally, all hits were parsed, using the AC number to retrieve the GO id from idmapping file of Uniprot, extract the Description of the protein from the corresponding dat file from uniprot, and pfam hit to the corresponding GO terms. All GO terms were combined, deplicates removed. 
  
  The left ones were blast to ncRNA database, to note as whether nc or c.
  
  These results were piped into a duckdb database for final outputs avoiding mysql issues.

  All large files were analysed in a way by seperated small trunks in case of memory problems.

## Database link

swissprotDB: [ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz) 

swissprotDBDAT: [ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz) 

tremblDB: [ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz)

tremblDBDAT: [ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz)

GODB: [ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz) 

ncRNADB: [ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz](ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/sequences/rnacentral_active.fasta.gz ) 

unirefDB: [ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.fasta.gz) 

unirefDBDAT: [ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.dat.gz](ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/uniref/uniref90/uniref90.dat.gz) 

pfam2GO: [http://geneontology.org/external2go/pfam2go](http://geneontology.org/external2go/pfam2go)

## Software

python=3.9
diamond
blast
transdecoder
fastparquet
duckdb
pyarrow
dask
hmmer
