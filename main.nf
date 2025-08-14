#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Log pipeline parameters at startup
log.info """
    ================================================================
     Gene Annotator Pro - Offline High-Precision Annotation
    ================================================================
    Input FASTA      : ${params.input}
    Output Directory : ${params.outdir}
    Prepared DB Dir  : ${params.prepared_db_dir}
    ================================================================
    """

// Helper function to check for input file
if (!params.input) {
    exit 1, "[ERROR] Input FASTA file not specified! Please run with --input <your_fasta_file>"
}

// ====================================================================================
// ==                                 Main Workflow                                  ==
// ====================================================================================
workflow {

    ch_input = Channel.fromPath(params.input)
    
    // ------------------------------------------------------------------
    // STAGE 1: Sequence Processing & Annotation
    // ------------------------------------------------------------------
    
    // Step 1.1: Predict proteins
    TRANSDECODER(ch_input)
    
    // Step 1.2: Hierarchical DIAMOND Search
    DIAMOND_SPROT(TRANSDECODER.out.pep)
    FILTER_FASTA_SPROT(TRANSDECODER.out.pep, DIAMOND_SPROT.out.tsv.collect().ifEmpty([]))
    
    DIAMOND_TREMBL(FILTER_FASTA_SPROT.out.unmatched_fasta)
    FILTER_FASTA_TREMBL(FILTER_FASTA_SPROT.out.unmatched_fasta, DIAMOND_TREMBL.out.tsv.collect().ifEmpty([]))

    DIAMOND_UNIREF90(FILTER_FASTA_TREMBL.out.unmatched_fasta)

    // Step 1.3: Run other protein annotations in parallel on ALL predicted proteins
    INTERPROSCAN(TRANSDECODER.out.pep)
    EGGNOG_MAPPER(TRANSDECODER.out.pep)
    
    // Step 1.4: Annotate non-coding candidates
    RFAM_SCAN(TRANSDECODER.out.ncrna_candidates)
    
    // ------------------------------------------------------------------
    // STAGE 2: FINAL INTEGRATION
    // ------------------------------------------------------------------
    ch_for_integration = Channel.empty()
        .mix(DIAMOND_SPROT.out.tsv.ifEmpty(null))
        .mix(DIAMOND_TREMBL.out.tsv.ifEmpty(null))
        .mix(DIAMOND_UNIREF90.out.tsv.ifEmpty(null))
        .mix(INTERPROSCAN.out.tsv.ifEmpty(null))
        .mix(EGGNOG_MAPPER.out.annotations.ifEmpty(null))
        .mix(RFAM_SCAN.out.tblout.ifEmpty(null))
        .toList()

    INTEGRATE_RESULTS(ch_for_integration, TRANSDECODER.out.all_ids)
}


// ====================================================================================
// ==                                Process Definitions                             ==
// ====================================================================================

// (All PREPARE_* and CLEAN_DBS processes have been removed)

process TRANSDECODER {
    label 'medium_cpu'
    publishDir "${params.outdir}/1_transdecoder", mode: 'copy', pattern: '*.{pep,fasta,txt}'
//    conda 'bioconda::transdecoder=5.7.0 bioconda::seqtk'
    input: path fasta
    output:
        path("*.transdecoder.pep"), emit: pep
        path("ncrna_candidates.fasta"), emit: ncrna_candidates
        path("all_query_ids.txt"), emit: all_ids
    script:
        def base = fasta.baseName
        """
        TransDecoder.LongOrfs -t ${fasta} -m ${params.transdecoder_min_len}
        TransDecoder.Predict -t ${fasta} --single_best_only
        grep '^>' ${base}.transdecoder.pep | sed 's/>//' | cut -d ' ' -f 1 | sed 's/\\.p[0-9]*\$//' | sort -u > coding_ids.txt
        grep '^>' ${fasta} | sed 's/>//' | sort -u > all_query_ids.txt
        grep -v -w -F -f coding_ids.txt all_query_ids.txt > ncrna_ids.txt
        seqtk subseq ${fasta} ncrna_ids.txt > ncrna_candidates.fasta
        """
}

process DIAMOND_SPROT {
    label 'high_cpu'
    publishDir "${params.outdir}/2_diamond/sprot", mode: 'copy'
//    conda 'bioconda::diamond=2.1.8'
    input: path pep_fasta
    output: path("*.sprot.tsv"), emit: tsv
    script:
        def base = pep_fasta.baseName.replaceAll("\\.transdecoder\\.pep", "")
        "diamond blastp --query ${pep_fasta} --db ${params.db_sprot_dmnd} --out ${base}.sprot.tsv --evalue ${params.diamond_evalue} --threads ${task.cpus} --max-target-seqs 1 --outfmt 6 qseqid sseqid pident evalue bitscore stitle"
}

process FILTER_FASTA_SPROT {
    label 'low_cpu'
    publishDir "${params.outdir}/1_transdecoder", mode: 'copy', pattern: "*.unmatched_sprot.fasta"
//    conda 'bioconda::seqtk'
    input: path original_fasta; path tsv_results
    output: path("*.unmatched_sprot.fasta"), emit: unmatched_fasta
    script:
        def base = original_fasta.baseName.replaceAll("\\.transdecoder\\.pep", "")
        """
        cut -f 1 ${tsv_results} | sed 's/\\.p[0-9]*\$//' | sort -u > matched_ids.txt
        grep '^>' ${original_fasta} | sed 's/>//' | sed 's/\\.p[0-9]*\$//' | sort -u > all_ids.txt
        grep -v -w -F -f matched_ids.txt all_ids.txt > unmatched_ids.txt
        seqtk subseq ${original_fasta} unmatched_ids.txt > ${base}.unmatched_sprot.fasta
        """
}

process DIAMOND_TREMBL {
    label 'high_cpu'
    publishDir "${params.outdir}/2_diamond/trembl", mode: 'copy'
//    conda 'bioconda::diamond=2.1.8'
    input: path pep_fasta
    output: path("*.trembl.tsv"), emit: tsv
    script:
        def base = pep_fasta.baseName.replaceAll("\\.unmatched_sprot", "")
        """
        if [ -s ${pep_fasta} ]; then
            diamond blastp --query ${pep_fasta} --db ${params.db_trembl_dmnd} --out ${base}.trembl.tsv --evalue ${params.diamond_evalue} --threads ${task.cpus} --max-target-seqs 1 --outfmt 6 qseqid sseqid pident evalue bitscore stitle
        else
            touch ${base}.trembl.tsv
        fi
        """
}

process FILTER_FASTA_TREMBL {
    label 'low_cpu'
    publishDir "${params.outdir}/1_transdecoder", mode: 'copy', pattern: "*.unmatched_trembl.fasta"
//    conda 'bioconda::seqtk'
    input: path original_fasta; path tsv_results
    output: path("*.unmatched_trembl.fasta"), emit: unmatched_fasta
    script:
        def base = original_fasta.baseName.replaceAll("\\.unmatched_sprot", "")
        """
        cut -f 1 ${tsv_results} | sed 's/\\.p[0-9]*\$//' | sort -u > matched_ids.txt
        grep '^>' ${original_fasta} | sed 's/>//' | sed 's/\\.p[0-9]*\$//' | sort -u > all_ids.txt
        grep -v -w -F -f matched_ids.txt all_ids.txt > unmatched_ids.txt
        seqtk subseq ${original_fasta} unmatched_ids.txt > ${base}.unmatched_trembl.fasta
        """
}

process DIAMOND_UNIREF90 {
    label 'high_cpu'
    publishDir "${params.outdir}/2_diamond/uniref90", mode: 'copy'
//    conda 'bioconda::diamond=2.1.8'
    input: path pep_fasta
    output: path("*.uniref90.tsv"), emit: tsv
    script:
        def base = pep_fasta.baseName.replaceAll("\\.unmatched_trembl", "")
        """
        if [ -s ${pep_fasta} ]; then
            diamond blastp --query ${pep_fasta} --db ${params.db_uniref90_dmnd} --out ${base}.uniref90.tsv --evalue ${params.diamond_evalue} --threads ${task.cpus} --max-target-seqs 1 --outfmt 6 qseqid sseqid pident evalue bitscore stitle
        else
            touch ${base}.uniref90.tsv
        fi
        """
}

process INTERPROSCAN {
    label 'super_high_cpu'
    publishDir "${params.outdir}/3_interproscan", mode: 'copy'
//    conda 'bioconda::interproscan' // Let conda find the best version

    input:
    path pep_fasta

    output:
    path("*.interpro.tsv"), emit: tsv

    script:
    def base = pep_fasta.baseName.replaceAll("\\.transdecoder\\.pep", "")
    def cleaned_pep = "cleaned.pep"
    """
    # Step 1: Clean the protein FASTA file by removing all asterisk (*) characters.
    # This is required because InterProScan does not accept stop codons in sequences.
    sed 's/\\*//g' ${pep_fasta} > ${cleaned_pep}

    # Step 2: Run interproscan.sh on the CLEANED file.
    interproscan.sh \\
        -i ${cleaned_pep} \\
        -d . \\
        -T . \\
        -appl Pfam,SMART,CDD,SUPERFAMILY,PRINTS \\
        --cpu ${task.cpus} \\
        --goterms \\
        --pathways \\
        --formats tsv \\
        --disable-precalc
    
    # Step 3: Rename the final output file.
    mv *.tsv ${base}.interpro.tsv
    """
}


process EGGNOG_MAPPER {
    label 'high_cpu'
    publishDir "${params.outdir}/4_eggnog", mode: 'copy'
//    conda 'bioconda::eggnog-mapper=2.1.9'
    input: path pep_fasta
    output: path("emapper.annotations"), emit: annotations
    script:
        def base = pep_fasta.baseName.replaceAll("\\.transdecoder\\.pep", "")
        """
        emapper.py \\
            -i ${pep_fasta} \\
            -o ${base} \\
            --output_dir . \\
            --data_dir ${params.db_eggnog_dir} \\
            --cpu ${task.cpus}
        
        mv ${base}.emapper.annotations emapper.annotations
        """
}

process RFAM_SCAN {
    label 'high_cpu'
    publishDir "${params.outdir}/5_rfam", mode: 'copy'
//    conda 'bioconda::infernal=1.1.4'
    input: path ncrna_fasta
    output: path("rfam.tblout"), emit: tblout
    script: "cmscan --cpu ${task.cpus} --tblout rfam.tblout ${params.db_rfam_cm} ${ncrna_fasta} > rfam.log"
}

process INTEGRATE_RESULTS {
    label 'low_cpu'
    publishDir "${params.outdir}/final_report", mode: 'copy'
    conda 'bioconda::pandas=2.1.0 conda-forge::duckdb=0.9.2'
    input: path results_list; path all_ids_file
    output: path "final_annotations.tsv"; path "annotation_summary.log"
    script:
        def sprot_file = results_list.find { it?.name.endsWith('.sprot.tsv') } ?: 'null_file'
        def trembl_file = results_list.find { it?.name.endsWith('.trembl.tsv') } ?: 'null_file'
        def uniref90_file = results_list.find { it?.name.endsWith('.uniref90.tsv') } ?: 'null_file'
        def interpro_file = results_list.find { it?.name.endsWith('.interpro.tsv') } ?: 'null_file'
        def eggnog_file = results_list.find { it?.name.endsWith('emapper.annotations') } ?: 'null_file'
        def rfam_file = results_list.find { it?.name.endsWith('.tblout') } ?: 'null_file'
        """
        python ${baseDir}/bin/integrate_results.py \\
            --query_ids ${all_ids_file} \\
            --sprot ${sprot_file} \\
            --trembl ${trembl_file} \\
            --uniref90 ${uniref90_file} \\
            --interpro ${interpro_file} \\
            --eggnog ${eggnog_file} \\
            --rfam ${rfam_file} \\
            --output final_annotations.tsv \\
            --summary annotation_summary.log
        """
}
