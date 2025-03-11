import os
import subprocess
import argparse
import yaml
import pandas as pd
import gzip
import shutil
from Bio import SeqIO
import logging
import duckdb
from pathlib import Path
import hashlib
import time

# logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# default config
DEFAULT_CONFIG = {
    "fasta_file": "input.fasta",
    "database": {
        "source_dir": "./databases",
        "uniprot_sprot": "uniprot_sprot.fasta.gz",
        "uniprot_trembl": "uniprot_trembl.fasta.gz",
        "uniref90": "uniref90.fasta.gz",
        "ncRNA": "rnacentral_active.fasta.gz",
        "idmapping": "idmapping_selected.tab.gz",
        "cdd": "Cdd_LE.tar.gz",
        "uniprot_sprot_dat": "uniprot_sprot.dat.gz",
        "uniprot_trembl_dat": "uniprot_trembl.dat.gz",
        "pfam": "Pfam-A.hmm.gz",
        "pfam2go": "pfam2go.txt"
    },
    "duckdb": {
        "path": "annotations.duckdb"
    },
    "diamond": {
        "evalue": 1e-5,
        "outfmt": 6
    },
    "blastn": {
        "evalue": 1e-5,
        "outfmt": 6
    },
    "transdecoder": {
        "min_length": 50
    },
    "output": {
        "dir": "./output",
        "sprot_out": "sprot.out",
        "trembl_out": "trembl.out",
        "uniref90_out": "uniref90.out",
        "ncRNA_out": "ncRNA.out",
        "annotations": "annotations.tsv"
    }
}

def generate_config_template(user_dir):
    config_path = os.path.join(user_dir, "config.yaml")
    if not os.path.exists(config_path):
        with open(config_path, "w") as f:
            yaml.dump(DEFAULT_CONFIG, f, default_flow_style=False)
        logging.info(f"Generated config template at {config_path}. Please edit it and rerun.")
        exit(0)

def load_config(user_dir):
    config_path = os.path.join(user_dir, "config.yaml")
    generate_config_template(user_dir)
    with open(config_path) as f:
        config = yaml.safe_load(f)
    required_keys = set(DEFAULT_CONFIG.keys())
    missing_keys = required_keys - set(config.keys())
    if missing_keys:
        logging.error(f"Config file missing required keys: {missing_keys}")
        raise KeyError(f"Invalid config file: missing {missing_keys}")
    logging.info(f"Loaded configuration from {config_path}")
    return config

def parse_arguments():
    parser = argparse.ArgumentParser(description="Plant Transcriptome Annotation Tool")
    parser.add_argument("--fasta", help="Input FASTA file", default=DEFAULT_CONFIG["fasta_file"])
    parser.add_argument("--db_dir", help="Database directory", default=DEFAULT_CONFIG["database"]["source_dir"])
    parser.add_argument("--threads", type=int, help="Number of CPU threads", default=4)
    args = parser.parse_args()
    logging.info(f"Arguments parsed: fasta={args.fasta}, db_dir={args.db_dir}, threads={args.threads}")
    return args

def check_files(config, args):
    db_dir = os.path.abspath(args.db_dir)
    if not os.path.exists(args.fasta):
        logging.error(f"Input FASTA file not found: {args.fasta}")
        raise FileNotFoundError(f"Input FASTA file not found: {args.fasta}")
    required_files = ["uniprot_sprot", "uniprot_trembl", "uniref90", "ncRNA", "idmapping", "cdd", "uniprot_sprot_dat", "uniprot_trembl_dat", "pfam", "pfam2go"]
    missing_files = [f for f in required_files if not os.path.exists(os.path.join(db_dir, config["database"][f]))]
    if missing_files:
        logging.error(f"Missing files in {db_dir}: {', '.join(missing_files)}")
        raise FileNotFoundError(f"Missing files in {db_dir}: {', '.join(missing_files)}")
    logging.info("All required database files found.")

def get_file_hash(file_path):
    logging.info(f"Computing MD5 hash for file: {file_path}")
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    hash_value = hash_md5.hexdigest()
    logging.info(f"MD5 hash computed: {hash_value}")
    return hash_value

def build_diamond_index(db_file, output_db, threads):
    output_db = os.path.abspath(output_db)
    dmnd_file = f"{output_db}.dmnd"
    if os.path.exists(dmnd_file):
        logging.info(f"DIAMOND index already exists and will be reused: {dmnd_file}")
        return output_db
    if db_file.endswith(".gz"):
        uncompressed = output_db + ".fasta"
        if not os.path.exists(uncompressed):
            with gzip.open(db_file, "rb") as f_in, open(uncompressed, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            logging.info(f"Decompressed {db_file} to {uncompressed}")
        db_file = uncompressed
    cmd = ["diamond", "makedb", "--in", db_file, "--db", output_db, "--threads", str(threads)]
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"DIAMOND index built: {dmnd_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to build DIAMOND index: {e}")
        raise
    return output_db

def run_diamond(query, db_index, output_file, config, threads, checkpoint_file):
    output_hash_file = checkpoint_file + ".hash"
    if os.path.exists(checkpoint_file) and os.path.exists(output_file):
        with open(output_hash_file, "r") as f:
            stored_hash = f.read().strip()
        current_hash = get_file_hash(output_file)
        if stored_hash == current_hash:
            logging.info(f"DIAMOND alignment already completed and verified: {output_file}")
            return output_file
    cmd = [
        "diamond", "blastx", "--query", query, "--db", db_index, "--out", output_file,
        "--evalue", str(config["diamond"]["evalue"]), "--threads", str(threads),
        "--outfmt", str(config["diamond"]["outfmt"])
    ]
    try:
        subprocess.run(cmd, check=True)
        current_hash = get_file_hash(output_file)
        with open(checkpoint_file, "w") as f:
            f.write("done")
        with open(output_hash_file, "w") as f:
            f.write(current_hash)
        logging.info(f"DIAMOND alignment completed: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"DIAMOND alignment failed: {e}")
        raise
    return output_file

def filter_unmatched(fasta_file, hits, output_fasta):
    matched_ids = set(hits["qseqid"]) if not hits.empty else set()
    total_sequences = 0
    unmatched_sequences = 0
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            total_sequences += 1
            if record.id not in matched_ids:
                SeqIO.write(record, out, "fasta")
                unmatched_sequences += 1
    if total_sequences == 0:
        logging.warning(f"Input FASTA file {fasta_file} is empty")
    logging.info(f"Processed {total_sequences} sequences, {unmatched_sequences} unmatched written to {output_fasta}")
    return output_fasta

def build_blast_index(db_file, output_db):
    output_db = os.path.abspath(output_db)
    nhr_file = f"{output_db}.00.nhr"
    if os.path.exists(nhr_file):
        logging.info(f"BLAST index already exists and will be reused: {nhr_file}")
        return output_db
    if db_file.endswith(".gz"):
        uncompressed = output_db + ".fasta"
        if not os.path.exists(uncompressed):
            with gzip.open(db_file, "rb") as f_in, open(uncompressed, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            logging.info(f"Decompressed {db_file} to {uncompressed}")
        db_file = uncompressed
    cmd = ["makeblastdb", "-in", db_file, "-dbtype", "nucl", "-out", output_db]
    try:
        subprocess.run(cmd, check=True)
        logging.info(f"BLAST index built: {output_db}")
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to build BLAST index: {e}")
        raise
    return output_db

def run_blastn(query, db_index, output_file, config, threads, checkpoint_file):
    output_hash_file = checkpoint_file + ".hash"
    if os.path.exists(checkpoint_file) and os.path.exists(output_file):
        with open(output_hash_file, "r") as f:
            stored_hash = f.read().strip()
        current_hash = get_file_hash(output_file)
        if stored_hash == current_hash:
            logging.info(f"BLASTN alignment already completed and verified: {output_file}")
            return output_file
    cmd = [
        "blastn", "-query", query, "-db", db_index, "-out", output_file,
        "-evalue", str(config["blastn"]["evalue"]), "-num_threads", str(threads),
        "-outfmt", str(config["blastn"]["outfmt"])
    ]
    try:
        subprocess.run(cmd, check=True)
        current_hash = get_file_hash(output_file)
        with open(checkpoint_file, "w") as f:
            f.write("done")
        with open(output_hash_file, "w") as f:
            f.write(current_hash)
        logging.info(f"BLASTN alignment completed: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"BLASTN alignment failed: {e}")
        raise
    return output_file

def run_transdecoder(fasta_file, config, output_dir):
    checkpoint_file = os.path.join(output_dir, "transdecoder.done")
    pep_file = os.path.join(output_dir, f"{os.path.basename(fasta_file)}.transdecoder.pep")
    output_hash_file = checkpoint_file + ".hash"
    if os.path.exists(checkpoint_file) and os.path.exists(pep_file):
        with open(output_hash_file, "r") as f:
            stored_hash = f.read().strip()
        current_hash = get_file_hash(pep_file)
        if stored_hash == current_hash:
            logging.info(f"TransDecoder already completed and verified: {pep_file}")
            return pep_file
    os.makedirs(output_dir, exist_ok=True)
    try:
        subprocess.run([
            "TransDecoder.LongOrfs", "-t", fasta_file, "-m", str(config["transdecoder"]["min_length"]),
            "--output_dir", output_dir
        ], check=True)
        subprocess.run([
            "TransDecoder.Predict", "-t", fasta_file, "--output_dir", output_dir, "--single_best_only"
        ], check=True)
        current_hash = get_file_hash(pep_file)
        with open(checkpoint_file, "w") as f:
            f.write("done")
        with open(output_hash_file, "w") as f:
            f.write(current_hash)
        logging.info(f"TransDecoder completed: input={fasta_file}, output={pep_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"TransDecoder failed: {e}")
        raise
    return pep_file

def run_hmmscan(query, hmm_file, output_file, threads, checkpoint_file):
    output_hash_file = checkpoint_file + ".hash"
    if os.path.exists(checkpoint_file) and os.path.exists(output_file):
        with open(output_hash_file, "r") as f:
            stored_hash = f.read().strip()
        current_hash = get_file_hash(output_file)
        if stored_hash == current_hash:
            logging.info(f"HMMSCAN already completed and verified: {output_file}")
            return output_file

    uncompressed_hmm = hmm_file.replace(".gz", "")
    if not os.path.exists(uncompressed_hmm):
        with gzip.open(hmm_file, "rb") as f_in, open(uncompressed_hmm, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        logging.info(f"Decompressed {hmm_file} to {uncompressed_hmm}")

    hmm_extensions = [".h3m", ".h3f", ".h3i", ".h3p"]
    index_complete = all(os.path.exists(uncompressed_hmm + ext) for ext in hmm_extensions)
    if not index_complete:
        logging.info(f"Running hmmpress on {uncompressed_hmm}")
        try:
            subprocess.run(["hmmpress", "-f", uncompressed_hmm], check=True)
            logging.info(f"Completed hmmpress for {uncompressed_hmm}")
        except subprocess.CalledProcessError as e:
            logging.error(f"hmmpress failed: {e}")
            raise

    cmd = ["hmmscan", "--cpu", str(threads), "-o", output_file, uncompressed_hmm, query]
    try:
        subprocess.run(cmd, check=True)
        current_hash = get_file_hash(output_file)
        with open(checkpoint_file, "w") as f:
            f.write("done")
        with open(output_hash_file, "w") as f:
            f.write(current_hash)
        logging.info(f"HMMSCAN completed: {output_file}")
    except subprocess.CalledProcessError as e:
        logging.error(f"HMMSCAN failed: {e}")
        raise
    return output_file

def parse_hmmscan_output(hmmscan_file):
    logging.info(f"Parsing HMMSCAN output from {hmmscan_file}")
    domains = {}
    unparsed_lines = 0
    with open(hmmscan_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip() or "E-value" in line or "best" in line:
                continue
            parts = line.split()
            if len(parts) >= 9:
                try:
                    qseqid, domain, score = parts[2], parts[0], float(parts[5])
                    if qseqid not in domains or domains[qseqid][1] < score:
                        domains[qseqid] = (domain, score)
                except (ValueError, IndexError):
                    unparsed_lines += 1
                    continue
            else:
                unparsed_lines += 1
    logging.info(f"Parsed {len(domains)} domains, skipped {unparsed_lines} unparsable lines")
    return domains

def parse_uniprot_dat(dat_file):
    logging.info(f"Parsing UniProt DAT file: {dat_file}")
    descriptions = {}
    with gzip.open(dat_file, "rt") as f:
        uniprot_id = None
        desc_lines = []
        for line in f:
            if line.startswith("AC"):
                ac_line = line[5:].strip().rstrip(";")
                uniprot_id = ac_line.split(";")[0].strip()
            elif line.startswith("DE"):
                desc_lines.append(line[5:].strip(";").strip())
            elif line.startswith("//"):
                if uniprot_id and desc_lines:
                    rec_name = next((d.split("Full=")[1] for d in desc_lines if "RecName: Full=" in d), None)
                    descriptions[uniprot_id] = rec_name if rec_name else " ".join(desc_lines)
                uniprot_id, desc_lines = None, []
    logging.info(f"Parsed {len(descriptions)} descriptions from {dat_file}")
    return descriptions

def parse_pfam2go(pfam2go_file):
    logging.info(f"Parsing Pfam2GO file: {pfam2go_file}")
    pfam2go = {}
    with open(pfam2go_file, "r") as f:
        for line in f:
            if line.startswith("!"):
                continue
            parts = line.strip().split(" > ")
            if len(parts) != 2:
                continue
            pfam_part = parts[0].split()[0]
            pfam_id = pfam_part.split(":")[1]
            go_part = parts[1]
            go_ids = [go.split("; ")[1].strip() for go in go_part.split(" ; ") if "; " in go]
            pfam2go[pfam_id] = ";".join(go_ids)
    logging.info(f"Parsed {len(pfam2go)} Pfam to GO mappings from {pfam2go_file}")
    return pfam2go

def initialize_database(config):
    conn = duckdb.connect(config["duckdb"]["path"])
    try:
        conn.execute("CREATE TABLE IF NOT EXISTS hits (qseqid VARCHAR, sseqid VARCHAR, evalue DOUBLE, bitscore DOUBLE, UniProtKB_AC VARCHAR, GeneID VARCHAR, GO VARCHAR, Description VARCHAR)")
        logging.info("Created 'hits' table in DuckDB")
        conn.execute("CREATE TABLE IF NOT EXISTS domains (qseqid VARCHAR, domain VARCHAR, GO VARCHAR)")
        logging.info("Created 'domains' table in DuckDB")
        conn.execute("CREATE TABLE IF NOT EXISTS ncRNA (qseqid VARCHAR, ncRNA VARCHAR)")
        logging.info("Created 'ncRNA' table in DuckDB")
        conn.execute("CREATE TABLE IF NOT EXISTS unknown (qseqid VARCHAR, function VARCHAR)")
        logging.info("Created 'unknown' table in DuckDB")
        conn.execute("CREATE TABLE IF NOT EXISTS descriptions (UniProtKB_AC VARCHAR, Description VARCHAR)")
        logging.info("Created 'descriptions' table in DuckDB")
    except duckdb.Error as e:
        logging.error(f"Failed to initialize DuckDB tables: {e}")
        raise
    return conn

def main():
    start_time = time.time()
    logging.info("Starting Plant Transcriptome Annotation Tool")
    
    args = parse_arguments()
    config = load_config(os.getcwd())
    check_files(config, args)
    output_dir = config["output"]["dir"]
    os.makedirs(output_dir, exist_ok=True)
    conn = initialize_database(config)

    pfam2go = parse_pfam2go(os.path.join(args.db_dir, config["database"]["pfam2go"]))

    # Step 1: Diamond BLASTX against uniprot_sprot
    logging.info("Starting Step 1: Diamond BLASTX against UniProt Swiss-Prot")
    step_start = time.time()
    sprot_db = build_diamond_index(os.path.join(args.db_dir, config["database"]["uniprot_sprot"]),
                                   os.path.join(args.db_dir, "sprot"), args.threads)
    sprot_out = run_diamond(args.fasta, sprot_db, os.path.join(output_dir, config["output"]["sprot_out"]),
                            config, args.threads, os.path.join(output_dir, "sprot.done"))
    sprot_hits = pd.read_csv(sprot_out, sep="\t", header=None,
                             names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    if "bitscore" not in sprot_hits.columns:
        logging.error("Bitscore column missing in Diamond output")
        raise ValueError("Invalid Diamond output format")
    sprot_best_chunk = sprot_hits.loc[sprot_hits.groupby("qseqid")["bitscore"].idxmax()]
    conn.register("sprot_best_chunk_temp", sprot_best_chunk)
    conn.execute("INSERT INTO hits (qseqid, sseqid, evalue, bitscore) SELECT qseqid, sseqid, evalue, bitscore FROM sprot_best_chunk_temp")
    conn.unregister("sprot_best_chunk_temp")
    logging.info(f"Step 1: Inserted {len(sprot_best_chunk)} records into 'hits' table")
    logging.info(f"Completed Step 1 in {time.time() - step_start:.2f} seconds")

    # Step 2: Filter unmatched and BLASTX against uniprot_trembl
    logging.info("Starting Step 2: Filter unmatched sequences and BLASTX against UniProt TrEMBL")
    step_start = time.time()
    unmatched_fasta = os.path.join(output_dir, "unmatched_sprot.fasta")
    sprot_best = conn.execute("SELECT DISTINCT qseqid FROM hits").fetchdf()
    filter_unmatched(args.fasta, sprot_best, unmatched_fasta)
    trembl_db = build_diamond_index(os.path.join(args.db_dir, config["database"]["uniprot_trembl"]),
                                    os.path.join(args.db_dir, "trembl"), args.threads)
    trembl_out = run_diamond(unmatched_fasta, trembl_db, os.path.join(output_dir, config["output"]["trembl_out"]),
                            config, args.threads, os.path.join(output_dir, "trembl.done"))
    trembl_hits = pd.read_csv(trembl_out, sep="\t", header=None,
                              names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    trembl_best_chunk = trembl_hits.loc[trembl_hits.groupby("qseqid")["bitscore"].idxmax()]
    conn.register("trembl_best_chunk_temp", trembl_best_chunk)
    conn.execute("INSERT INTO hits (qseqid, sseqid, evalue, bitscore) SELECT qseqid, sseqid, evalue, bitscore FROM trembl_best_chunk_temp")
    conn.unregister("trembl_best_chunk_temp")
    logging.info(f"Step 2: Inserted {len(trembl_best_chunk)} records into 'hits' table")
    logging.info(f"Completed Step 2 in {time.time() - step_start:.2f} seconds")

    # Step 3: Filter unmatched and BLASTX against uniref90
    logging.info("Starting Step 3: Filter unmatched sequences and BLASTX against UniRef90")
    step_start = time.time()
    unmatched_fasta2 = os.path.join(output_dir, "unmatched_trembl.fasta")
    trembl_best = conn.execute("SELECT DISTINCT qseqid FROM hits WHERE qseqid NOT IN (SELECT qseqid FROM hits WHERE evalue IS NOT NULL AND sseqid LIKE '%|sp|%')").fetchdf()
    filter_unmatched(unmatched_fasta, trembl_best, unmatched_fasta2)
    uniref90_db = build_diamond_index(os.path.join(args.db_dir, config["database"]["uniref90"]),
                                      os.path.join(args.db_dir, "uniref90"), args.threads)
    uniref90_out = run_diamond(unmatched_fasta2, uniref90_db, os.path.join(output_dir, config["output"]["uniref90_out"]),
                               config, args.threads, os.path.join(output_dir, "uniref90.done"))
    uniref90_hits = pd.read_csv(uniref90_out, sep="\t", header=None,
                                names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    uniref90_best_chunk = uniref90_hits.loc[uniref90_hits.groupby("qseqid")["bitscore"].idxmax()]
    conn.register("uniref90_best_chunk_temp", uniref90_best_chunk)
    conn.execute("INSERT INTO hits (qseqid, sseqid, evalue, bitscore) SELECT qseqid, sseqid, evalue, bitscore FROM uniref90_best_chunk_temp")
    conn.unregister("uniref90_best_chunk_temp")
    logging.info(f"Step 3: Inserted {len(uniref90_best_chunk)} records into 'hits' table")
    logging.info(f"Completed Step 3 in {time.time() - step_start:.2f} seconds")

    # Step 4: Annotate with idmapping and extract Accession Number from sseqid
    logging.info("Starting Step 4: Annotate hits with UniProt ID mapping")
    step_start = time.time()
    idmapping_file = os.path.join(args.db_dir, config["database"]["idmapping"])
    conn.execute("CREATE TABLE IF NOT EXISTS idmapping (UniProtKB_AC VARCHAR, GeneID VARCHAR, GO VARCHAR)")
    logging.info("Created 'idmapping' table in DuckDB")
    for chunk in pd.read_csv(idmapping_file, sep="\t", header=None, usecols=[0, 2, 6], names=["UniProtKB_AC", "GeneID", "GO"], compression="gzip", chunksize=10000):
        conn.register("idmapping_chunk_temp", chunk)
        conn.execute("INSERT INTO idmapping SELECT * FROM idmapping_chunk_temp")
        conn.unregister("idmapping_chunk_temp")
    conn.execute("""
        UPDATE hits 
        SET UniProtKB_AC = REGEXP_EXTRACT(sseqid, '\|([A-Za-z0-9]+)\|', 1)
        WHERE UniProtKB_AC IS NULL
    """)
    conn.execute("""
        UPDATE hits 
        SET GeneID = (SELECT GeneID FROM idmapping WHERE idmapping.UniProtKB_AC = hits.UniProtKB_AC LIMIT 1),
            GO = (SELECT GO FROM idmapping WHERE idmapping.UniProtKB_AC = hits.UniProtKB_AC LIMIT 1)
    """)
    logging.info(f"Step 4: Updated UniProtKB_AC, GeneID, and GO for hits")
    logging.info(f"Completed Step 4 in {time.time() - step_start:.2f} seconds")

    # Step 5: Add description from uniprot dat files
    logging.info("Starting Step 5: Add descriptions from UniProt DAT files")
    step_start = time.time()
    sprot_desc = parse_uniprot_dat(os.path.join(args.db_dir, config["database"]["uniprot_sprot_dat"]))
    trembl_desc = parse_uniprot_dat(os.path.join(args.db_dir, config["database"]["uniprot_trembl_dat"]))
    conn.execute("DELETE FROM descriptions")
    desc_df = pd.DataFrame([(k, v) for d in [sprot_desc, trembl_desc] for k, v in d.items()], columns=["UniProtKB_AC", "Description"])
    conn.register("desc_df_temp", desc_df)
    conn.execute("INSERT INTO descriptions SELECT * FROM desc_df_temp")
    conn.unregister("desc_df_temp")
    conn.execute("UPDATE hits SET Description = (SELECT Description FROM descriptions WHERE descriptions.UniProtKB_AC = hits.UniProtKB_AC LIMIT 1)")
    conn.execute("UPDATE hits SET Description = 'N/A' WHERE Description IS NULL")
    desc_count = conn.execute("SELECT COUNT(*) FROM hits WHERE Description != 'N/A'").fetchone()[0]
    logging.info(f"Step 5: {desc_count} hits updated with descriptions")
    logging.info(f"Completed Step 5 in {time.time() - step_start:.2f} seconds")

    # Step 6: Process ALL sequences with TransDecoder and Pfam domain search, add GO from pfam2go
    logging.info("Starting Step 6: Process all sequences with TransDecoder and Pfam domain search")
    step_start = time.time()
    pep_file = run_transdecoder(args.fasta, config, os.path.join(output_dir, "transdecoder_all"))
    hmmscan_out = run_hmmscan(pep_file, os.path.join(args.db_dir, config["database"]["pfam"]),
                              os.path.join(output_dir, "pfam_all.out"), args.threads, os.path.join(output_dir, "pfam_all.done"))
    domains = parse_hmmscan_output(hmmscan_out)
    domain_df = pd.DataFrame([(k, v[0], pfam2go.get(v[0], "")) for k, v in domains.items()], 
                             columns=["qseqid", "domain", "GO"])
    conn.execute("DELETE FROM domains")
    conn.register("domain_df_temp", domain_df)
    conn.execute("INSERT INTO domains SELECT * FROM domain_df_temp")
    conn.unregister("domain_df_temp")
    logging.info(f"Step 6: Inserted {len(domain_df)} domains with GO annotations into 'domains' table")
    logging.info(f"Completed Step 6 in {time.time() - step_start:.2f} seconds")

    # Step 7: BLASTN against ncRNA for remaining unmatched
    logging.info("Starting Step 7: BLASTN against ncRNA for remaining unmatched sequences")
    step_start = time.time()
    unmatched_fasta4 = os.path.join(output_dir, "unmatched_domains.fasta")
    all_annotated = conn.execute("SELECT DISTINCT qseqid FROM hits UNION SELECT qseqid FROM domains").fetchdf()
    filter_unmatched(args.fasta, all_annotated, unmatched_fasta4)
    ncRNA_db = build_blast_index(os.path.join(args.db_dir, config["database"]["ncRNA"]),
                                 os.path.join(args.db_dir, "ncRNA"))
    ncRNA_out = run_blastn(unmatched_fasta4, ncRNA_db, os.path.join(output_dir, config["output"]["ncRNA_out"]),
                           config, args.threads, os.path.join(output_dir, "ncRNA.done"))
    ncRNA_hits = pd.read_csv(ncRNA_out, sep="\t", header=None,
                             names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    ncRNA_best_chunk = ncRNA_hits.loc[ncRNA_hits.groupby("qseqid")["bitscore"].idxmax()]
    ncRNA_best_chunk["ncRNA"] = "Yes"
    conn.register("ncRNA_best_chunk_temp", ncRNA_best_chunk)
    conn.execute("INSERT INTO ncRNA (qseqid, ncRNA) SELECT qseqid, 'Yes' FROM ncRNA_best_chunk_temp")
    conn.unregister("ncRNA_best_chunk_temp")
    logging.info(f"Step 7: Inserted {len(ncRNA_best_chunk)} records into 'ncRNA' table")
    logging.info(f"Completed Step 7 in {time.time() - step_start:.2f} seconds")

    # Step 8: Mark remaining as unknown
    logging.info("Starting Step 8: Mark remaining sequences as unknown")
    step_start = time.time()
    all_transcripts = set(record.id for record in SeqIO.parse(args.fasta, "fasta"))
    annotated_ids = set(conn.execute("SELECT qseqid FROM hits UNION SELECT qseqid FROM domains UNION SELECT qseqid FROM ncRNA").fetchdf()["qseqid"])
    unknown_ids = all_transcripts - annotated_ids
    unknown_df = pd.DataFrame({"qseqid": list(unknown_ids), "function": "unknown"})
    conn.register("unknown_df_temp", unknown_df)
    conn.execute("INSERT INTO unknown SELECT * FROM unknown_df_temp")
    conn.unregister("unknown_df_temp")
    logging.info(f"Step 8: Marked {len(unknown_ids)} sequences as unknown")
    logging.info(f"Completed Step 8 in {time.time() - step_start:.2f} seconds")

    # Step 9: Export final results with merged and deduplicated GO information
    logging.info("Starting Step 9: Export final annotation results")
    step_start = time.time()
    output_file = os.path.join(output_dir, config["output"]["annotations"])
    conn.execute(f"""
        COPY (
            SELECT 
                h.qseqid, 
                h.sseqid, 
                h.evalue, 
                h.bitscore, 
                h.UniProtKB_AC, 
                h.GeneID, 
                STRING_AGG(DISTINCT go_item, ';') AS GO,
                h.Description, 
                d.domain AS Pfam_Domain
            FROM hits h
            LEFT JOIN domains d ON h.qseqid = d.qseqid
            LEFT JOIN UNNEST(SPLIT(COALESCE(h.GO, '') || ';' || COALESCE(d.GO, ''), ';')) AS t(go_item) ON go_item != ''
            GROUP BY h.qseqid, h.sseqid, h.evalue, h.bitscore, h.UniProtKB_AC, h.GeneID, h.Description, d.domain
            UNION ALL
            SELECT 
                d.qseqid, 
                NULL, NULL, NULL, NULL, NULL, 
                d.GO, 
                NULL, 
                d.domain
            FROM domains d
            WHERE d.qseqid NOT IN (SELECT qseqid FROM hits)
            UNION ALL
            SELECT 
                n.qseqid, 
                NULL, NULL, NULL, NULL, NULL, 
                NULL, 
                NULL, 
                NULL
            FROM ncRNA n
            WHERE n.qseqid NOT IN (SELECT qseqid FROM hits UNION SELECT qseqid FROM domains)
            UNION ALL
            SELECT 
                u.qseqid, 
                NULL, NULL, NULL, NULL, NULL, 
                NULL, 
                NULL, 
                NULL
            FROM unknown u
            WHERE u.qseqid NOT IN (SELECT qseqid FROM hits UNION SELECT qseqid FROM domains UNION SELECT qseqid FROM ncRNA)
        ) TO '{output_file}' (DELIMITER '\t', HEADER)
    """)
    logging.info(f"Step 9: Exported annotations to {output_file} with merged GO terms")
    logging.info(f"Completed Step 9 in {time.time() - step_start:.2f} seconds")

    conn.close()
    logging.info(f"Finished Plant Transcriptome Annotation Tool in {time.time() - start_time:.2f} seconds")

if __name__ == "__main__":
    main()
