#!/usr/bin/env python
import pandas as pd
import argparse
import logging
import os
import sys
import numpy as np

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] - %(message)s")

def file_is_valid(path):
    """Check if a file path is valid and not empty."""
    if path is None or path == 'null_file' or not os.path.exists(path):
        logging.warning(f"File {path} is invalid or does not exist.")
        return False
    if os.path.getsize(path) == 0:
        logging.warning(f"File {path} is empty.")
        return False
    return True

def add_transcript_id(df, column='qseqid_raw'):
    """Creates a 'transcript_id' column by cleaning the protein ID."""
    if column in df.columns and not df.empty:
        df['transcript_id'] = df[column].str.replace(r'\.p\d+$', '', regex=True)
    else:
        # This case is expected if a file is empty, so no warning is needed.
        pass
    return df

def parse_rfam_tblout(filepath):
    """Robustly parses the Rfam tblout format, keeping only the best hit."""
    logging.info(f"Robustly parsing Rfam hits from {filepath}")
    data = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#'): continue
                parts = line.split()
                if len(parts) < 18: continue
                transcript_id = parts[2]
                description = " ".join(parts[17:])
                data.append({'transcript_id': transcript_id, 'rfam_description': description})
        df = pd.DataFrame(data)
        if not df.empty:
            df.drop_duplicates(subset='transcript_id', keep='first', inplace=True)
        else:
            logging.warning(f"No valid data parsed from {filepath}.")
        return df
    except Exception as e:
        logging.error(f"Failed to parse Rfam file {filepath}: {e}")
        return pd.DataFrame(columns=['transcript_id', 'rfam_description'])

def apply_prioritization(row):
    """
    Applies the expert hierarchical logic to a single row to determine the best annotation.
    """
    # Helper function to check if a value is informative
    def is_informative(value):
        return pd.notna(value) and value != '-' and 'hypothetical' not in str(value).lower() and 'uncharacterized' not in str(value).lower()

    # Priority 1: Swiss-Prot
    if 'sprot_desc' in row and is_informative(row['sprot_desc']):
        return 'Protein (Swiss-Prot Hit)', row['sprot_desc']
    # Priority 2: TrEMBL
    if 'trembl_desc' in row and is_informative(row['trembl_desc']):
        return 'Protein (TrEMBL Hit)', row['trembl_desc']
    # Priority 3: eggNOG
    if 'EggNOG_Description' in row and is_informative(row['EggNOG_Description']):
        return 'Protein (eggNOG Hit)', row['EggNOG_Description']
    # Priority 4: InterProScan
    if 'signature_desc' in row and is_informative(row['signature_desc']):
        return 'Protein (InterProScan Domain)', row['signature_desc']
    # Priority 5: UniRef90 (last resort for protein)
    if 'uniref90_desc' in row and is_informative(row['uniref90_desc']):
        return 'Protein (UniRef90 Hit)', row['uniref90_desc']
    # Priority 6: Rfam
    if 'rfam_description' in row and pd.notna(row['rfam_description']):
        return 'Non-coding RNA (Rfam)', row['rfam_description']
    # Default
    return 'Unknown', 'No description found'

def main(args):
    logging.info("Starting definitive annotation integration using a robust merge-based strategy.")
    
    try:
        final_df = pd.read_csv(args.query_ids, header=None, names=['Transcript_ID']).set_index('Transcript_ID')
        total_queries = len(final_df)
        logging.info(f"Total sequences to annotate: {total_queries}")
    except Exception as e:
        logging.error(f"Could not read the query IDs file: {args.query_ids}. Error: {e}")
        sys.exit(1)

    # --- Load all data sources and prepare them with a 'transcript_id' index ---
    if file_is_valid(args.sprot):
        try:
            sprot_df = pd.read_csv(args.sprot, sep='\t', header=None, names=['qseqid_raw', 'sseqid', 'pident', 'evalue', 'bitscore', 'sprot_desc'])
            sprot_df = add_transcript_id(sprot_df).set_index('transcript_id')
            final_df = final_df.join(sprot_df[['sprot_desc']])
        except Exception as e:
            logging.warning(f"Failed to process Swiss-Prot file {args.sprot}: {e}")

    if file_is_valid(args.trembl):
        try:
            trembl_df = pd.read_csv(args.trembl, sep='\t', header=None, names=['qseqid_raw', 'sseqid', 'pident', 'evalue', 'bitscore', 'trembl_desc'])
            trembl_df = add_transcript_id(trembl_df).set_index('transcript_id')
            final_df = final_df.join(trembl_df[['trembl_desc']])
        except Exception as e:
            logging.warning(f"Failed to process TrEMBL file {args.trembl}: {e}")

    if file_is_valid(args.uniref90):
        try:
            uniref90_df = pd.read_csv(args.uniref90, sep='\t', header=None, names=['qseqid_raw', 'sseqid', 'pident', 'evalue', 'bitscore', 'uniref90_desc'])
            uniref90_df = add_transcript_id(uniref90_df).set_index('transcript_id')
            final_df = final_df.join(uniref90_df[['uniref90_desc']])
        except Exception as e:
            logging.warning(f"Failed to process UniRef90 file {args.uniref90}: {e}")

    if file_is_valid(args.eggnog):
        try:
            eggnog_df = pd.read_csv(args.eggnog, sep='\t', comment='#', header=None, usecols=[0, 6, 7, 9, 11, 12], names=['qseqid_raw', 'COG_category', 'EggNOG_Description', 'GOs', 'KEGG_ko', 'KEGG_Pathway'])
            eggnog_df = add_transcript_id(eggnog_df).set_index('transcript_id')
            final_df = final_df.join(eggnog_df[['COG_category', 'EggNOG_Description', 'GOs', 'KEGG_ko', 'KEGG_Pathway']])
        except Exception as e:
            logging.warning(f"Failed to process eggNOG file {args.eggnog}: {e}")

    if file_is_valid(args.interpro):
        try:
            interpro_df_raw = pd.read_csv(args.interpro, sep='\t', header=None, usecols=[0, 5, 11, 12, 13], names=['qseqid_raw', 'signature_desc', 'interpro_id', 'interpro_desc', 'go_terms'])
            interpro_df_raw = add_transcript_id(interpro_df_raw)
            interpro_df = interpro_df_raw.groupby('transcript_id').agg({'interpro_id': lambda x: ';'.join(x.dropna().astype(str).unique()), 'interpro_desc': lambda x: ';'.join(x.dropna().astype(str).unique()), 'signature_desc': 'first', 'go_terms': lambda x: '|'.join(x.dropna().astype(str).unique())})
            interpro_df['InterPro_Hits'] = interpro_df['interpro_id'] + " (" + interpro_df['interpro_desc'] + ")"
            final_df = final_df.join(interpro_df[['InterPro_Hits', 'signature_desc', 'go_terms']])
        except Exception as e:
            logging.warning(f"Failed to process InterPro file {args.interpro}: {e}")

    if file_is_valid(args.rfam):
        try:
            rfam_df = parse_rfam_tblout(args.rfam).set_index('transcript_id')
            final_df = final_df.join(rfam_df)
        except Exception as e:
            logging.warning(f"Failed to process Rfam file {args.rfam}: {e}")
    
    # --- Apply the Expert Prioritization Logic row-by-row ---
    final_df.reset_index(inplace=True) 
    logging.info("Applying expert prioritization logic to each transcript.")
    final_annotations = final_df.apply(apply_prioritization, axis=1, result_type='expand')
    final_df[['Best_Annotation_Source', 'Functional_Description']] = final_annotations

    # --- Final Cleanup and GO Term Merging ---
    logging.info("Performing final cleanup and formatting.")
    final_df['eggNOG_GOs'] = final_df.get('GOs', pd.Series(dtype=str)).fillna('').str.replace(',', '|')
    final_df['InterPro_GOs'] = final_df.get('go_terms', pd.Series(dtype=str)).fillna('')
    go_series = final_df['eggNOG_GOs'] + '|' + final_df['InterPro_GOs']
    final_df['Final_GO_Terms'] = go_series.apply(lambda x: ';'.join(sorted(list(set(term for term in x.split('|') if term and term != '-')))))
    
    final_columns = [
        'Transcript_ID', 'Best_Annotation_Source', 'Functional_Description', 'InterPro_Hits', 
        'Final_GO_Terms', 'COG_category', 'KEGG_ko', 'KEGG_Pathway'
    ]
    for col in final_columns:
        if col not in final_df.columns:
            final_df[col] = '-'
    
    final_df = final_df[final_columns].fillna('-')
    
    logging.info(f"Writing final annotations to {args.output}")
    final_df.to_csv(args.output, sep='\t', index=False, na_rep='-')

    # --- Generate Summary ---
    logging.info("Generating annotation summary log.")
    annotated_count = len(final_df[final_df['Best_Annotation_Source'] != 'Unknown'])
    percentage = (annotated_count / total_queries) if total_queries > 0 else 0
    
    summary = f"""
    Annotation Summary:
    ---------------------------------
    Total Input Sequences: {total_queries}
    Total Annotated Sequences: {annotated_count} ({percentage:.2%})
    
    Breakdown by Best Annotation Source:
    ------------------------------------
    """
    summary += final_df['Best_Annotation_Source'].value_counts().to_string()
    
    with open(args.summary, 'w') as f:
        f.write(summary)
    logging.info(f"Summary written to {args.summary}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Integrate annotation results from a Nextflow pipeline.")
    parser.add_argument('--query_ids', required=True, help='Path to file with all original transcript IDs.')
    parser.add_argument('--sprot', help='Path to Swiss-Prot DIAMOND results')
    parser.add_argument('--trembl', help='Path to TrEMBL DIAMOND results')
    parser.add_argument('--uniref90', help='Path to UniRef90 DIAMOND results')
    parser.add_argument('--interpro', help='Path to InterProScan results')
    parser.add_argument('--eggnog', help='Path to eggNOG-mapper annotations')
    parser.add_argument('--rfam', help='Path to Rfam scan results')
    parser.add_argument('--output', required=True, help='Path to the final output TSV file')
    parser.add_argument('--summary', required=True, help='Path to the summary log file')
    
    args = parser.parse_args()
    main(args)
