#!/usr/bin/env python
import duckdb
import pandas as pd
import argparse
import logging
import os

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] - %(message)s")

def file_is_valid(path):
    """Check if a file path is not 'null_file', exists, and is not empty."""
    return path is not None and path != 'null_file' and os.path.exists(path) and os.path.getsize(path) > 0

def clean_protein_id(df, column='qseqid_raw'):
    """Remove .pX suffix from TransDecoder protein IDs."""
    df['qseqid'] = df[column].str.replace(r'\.p\d+$', '', regex=True)
    return df

def main(args):
    logging.info("Initializing in-memory DuckDB for smart annotation integration.")
    con = duckdb.connect(database=':memory:', read_only=False)

    # --- Load all data into staging tables ---
    logging.info(f"Loading all query IDs from {args.query_ids}")
    query_ids_df = pd.read_csv(args.query_ids, header=None, names=['qseqid'])
    con.register('base_query_ids', query_ids_df)
    total_queries = len(query_ids_df)
    logging.info(f"Total sequences to annotate: {total_queries}")

    # --- Protein Homology Hits (Hierarchical) ---
    if file_is_valid(args.sprot):
        logging.info(f"Loading Swiss-Prot hits from {args.sprot}")
        sprot_df = pd.read_csv(args.sprot, sep='\t', header=None, names=['qseqid_raw', 'sseqid', 'pident', 'evalue', 'bitscore', 'stitle'])
        sprot_df = clean_protein_id(sprot_df)
        con.register('stg_sprot', sprot_df)

    if file_is_valid(args.trembl):
        logging.info(f"Loading TrEMBL hits from {args.trembl}")
        trembl_df = pd.read_csv(args.trembl, sep='\t', header=None, names=['qseqid_raw', 'sseqid', 'pident', 'evalue', 'bitscore', 'stitle'])
        trembl_df = clean_protein_id(trembl_df)
        con.register('stg_trembl', trembl_df)

    if file_is_valid(args.uniref90):
        logging.info(f"Loading UniRef90 hits from {args.uniref90}")
        uniref90_df = pd.read_csv(args.uniref90, sep='\t', header=None, names=['qseqid_raw', 'sseqid', 'pident', 'evalue', 'bitscore', 'stitle'])
        uniref90_df = clean_protein_id(uniref90_df)
        con.register('stg_uniref90', uniref90_df)

    # --- Other Annotation Sources ---
    if file_is_valid(args.eggnog):
        logging.info(f"Loading eggNOG annotations from {args.eggnog}")
        eggnog_df = pd.read_csv(args.eggnog, sep='\t', comment='#', header=None,
                                names=['qseqid_raw', 'seed_ortholog', 'evalue', 'score', 'eggNOG_OGs', 'max_annot_lvl', 
                                       'COG_category', 'Description', 'Preferred_name', 'GOs', 'EC', 'KEGG_ko', 
                                       'KEGG_Pathway', 'KEGG_Module', 'KEGG_Reaction', 'KEGG_rclass', 'BRITE', 
                                       'KEGG_TC', 'CAZy', 'BiGG_Reaction', 'PFAMs'])
        eggnog_df = clean_protein_id(eggnog_df)
        con.register('stg_eggnog', eggnog_df)

    if file_is_valid(args.interpro):
        logging.info(f"Loading InterProScan results from {args.interpro}")
        interpro_df = pd.read_csv(args.interpro, sep='\t', header=None,
                                  names=['qseqid_raw', 'md5', 'len', 'analysis', 'signature_id', 'description',
                                         'start', 'stop', 'score', 'status', 'date', 'interpro_id',
                                         'interpro_desc', 'go_terms', 'pathway_info'],
                                  usecols=['qseqid_raw', 'interpro_id', 'interpro_desc', 'go_terms', 'pathway_info'])
        interpro_df = clean_protein_id(interpro_df)
        con.register('stg_interpro', interpro_df)
        
    if file_is_valid(args.rfam):
        logging.info(f"Loading Rfam hits from {args.rfam}")
        rfam_df = pd.read_csv(args.rfam, delim_whitespace=True, comment='#', header=None,
                              usecols=[0, 2], names=['qseqid', 'rfam_hit'])
        con.register('stg_rfam', rfam_df)

    # --- The "Golden" Integration Query using Common Table Expressions (CTEs) ---
    logging.info("Executing final integration query. This may take a moment.")
    final_query = """
    WITH
    -- 1. Combine all diamond hits with a source priority
    diamond_hits AS (
        SELECT qseqid, stitle, 'Swiss-Prot' as source, bitscore, evalue FROM stg_sprot
        UNION ALL
        SELECT qseqid, stitle, 'TrEMBL' as source, bitscore, evalue FROM stg_trembl
        UNION ALL
        SELECT qseqid, stitle, 'UniRef90' as source, bitscore, evalue FROM stg_uniref90
    ),
    -- 2. Select the single best diamond hit based on our priority order and bitscore
    diamond_best AS (
        SELECT DISTINCT ON (qseqid) *
        FROM diamond_hits
        ORDER BY qseqid,
                 CASE source
                     WHEN 'Swiss-Prot' THEN 1
                     WHEN 'TrEMBL' THEN 2
                     WHEN 'UniRef90' THEN 3
                     ELSE 4
                 END,
                 bitscore DESC
    ),
    -- 3. Aggregate InterProScan results
    interpro_agg AS (
        SELECT
            qseqid,
            STRING_AGG(DISTINCT interpro_id || ' (' || interpro_desc || ')', ';') AS InterPro_Hits,
            STRING_AGG(DISTINCT go_terms, '|') AS InterPro_GOs
        FROM stg_interpro
        WHERE interpro_id IS NOT NULL AND interpro_id != '-'
        GROUP BY qseqid
    ),
    -- 4. Select key eggNOG fields
    eggnog_agg AS (
        SELECT
            qseqid,
            Description AS EggNOG_Description,
            COG_category,
            KEGG_ko,
            KEGG_Pathway,
            REPLACE(GOs, ',', '|') AS EggNOG_GOs
        FROM stg_eggnog
    ),
    -- 5. Combine GO terms from both InterPro and eggNOG and deduplicate
    combined_go AS (
        SELECT
            qseqid,
            STRING_AGG(DISTINCT go_term, ';') AS Final_GO_Terms
        FROM (
            SELECT qseqid, UNNEST(STRING_SPLIT(InterPro_GOs, '|')) AS go_term FROM interpro_agg WHERE InterPro_GOs IS NOT NULL
            UNION ALL
            SELECT qseqid, UNNEST(STRING_SPLIT(EggNOG_GOs, '|')) AS go_term FROM eggnog_agg WHERE EggNOG_GOs IS NOT NULL
        )
        WHERE go_term != '' AND go_term IS NOT NULL
        GROUP BY qseqid
    ),
    -- 6. Get best Rfam hit for non-coding RNAs
    rfam_best AS (
        SELECT DISTINCT ON (qseqid) qseqid, rfam_hit FROM stg_rfam
    )
    -- Final SELECT to build the report from the base list of all query IDs
    SELECT
        b.qseqid AS Transcript_ID,
        -- Annotation Source Logic (more detailed)
        CASE
            WHEN rf.rfam_hit IS NOT NULL THEN 'Non-coding RNA (Rfam)'
            WHEN db.source IS NOT NULL THEN 'Protein (' || db.source || ' Hit)'
            WHEN ip.InterPro_Hits IS NOT NULL THEN 'Protein (InterProScan Hit)'
            WHEN eg.EggNOG_Description IS NOT NULL THEN 'Protein (eggNOG Hit)'
            ELSE 'Unknown'
        END AS Best_Annotation_Source,
        -- Prioritized Functional Description
        COALESCE(db.stitle, eg.EggNOG_Description, rf.rfam_hit, 'No description found') AS Functional_Description,
        ip.InterPro_Hits,
        go.Final_GO_Terms,
        eg.COG_category,
        eg.KEGG_ko,
        eg.KEGG_Pathway,
        db.source AS Diamond_Hit_Source,
        db.bitscore AS Diamond_Bitscore,
        db.evalue AS Diamond_Evalue
    FROM base_query_ids b
    LEFT JOIN diamond_best db ON b.qseqid = db.qseqid
    LEFT JOIN interpro_agg ip ON b.qseqid = ip.qseqid
    LEFT JOIN eggnog_agg eg ON b.qseqid = eg.qseqid
    LEFT JOIN combined_go go ON b.qseqid = go.qseqid
    LEFT JOIN rfam_best rf ON b.qseqid = rf.qseqid
    ORDER BY b.qseqid;
    """
    
    # Execute query and save results
    result_df = con.execute(final_query).fetchdf()
    logging.info(f"Integration complete. Writing {len(result_df)} annotations to {args.output}")
    result_df.to_csv(args.output, sep='\t', index=False, na_rep='-')

    # --- Generate Summary ---
    logging.info("Generating annotation summary log.")
    annotated_count = len(result_df[result_df['Best_Annotation_Source'] != 'Unknown'])
    summary = f"""
    Annotation Summary:
    ---------------------------------
    Total Input Sequences: {total_queries}
    Total Annotated Sequences: {annotated_count} ({annotated_count/total_queries:.2%})
    
    Breakdown by Best Annotation Source:
    ------------------------------------
    """
    summary += result_df['Best_Annotation_Source'].value_counts().to_string()
    
    with open(args.summary, 'w') as f:
        f.write(summary)
    logging.info(f"Summary written to {args.summary}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Integrate annotation results from a Nextflow pipeline.")
    parser.add_argument('--query_ids', required=True, help='Path to file with all original query IDs.')
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
