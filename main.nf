#!/usr/bin/env nextflow

/*
  Geneshot: A pipeline to robustly identify which alleles (n.e.e peptide coding sequences)
  are present in a microbial community.

  There is an optional preprocessing submodule, 
  followed by assembly (spades) + extraction of peptide coding sequences from the contigs
  The short reads are then aligned against the assembled peptides plus uniref100.
  We use the FAMLI algorithm to adjuticate these alignments.
  Annotations can follow. 
*/

// Using DSL-2
nextflow.enable.dsl=2

// Default values for boolean flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below

params.help = false
params.output = './results/'
params.manifest = null

// Flow control
params.nopreprocess = false
params.nocomposition = false

// Preprocessing options
params.host_index_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz'
params.host_index = false
params.min_host_align_score = 30
params.savereads = false

// Assembly-based Allele catalog options
params.phred_offset = 33 // spades

// Quantification via Alignment options
params.dmnd_min_identity = 90 // DIAMOND
params.dmnd_min_coverage = 80 // DIAMOND
params.dmnd_top_pct = 1 // DIAMOND
params.dmnd_min_score = 20 // DIAMOND
params.gencode = 11 //DIAMOND
params.sd_mean_cutoff = 3.0 // FAMLI

// Annotation options
params.noannot = false
params.taxonomic_dmnd = false
params.ncbi_taxdump = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
params.eggnog_db = false
params.eggnog_dmnd = false

// CAG options
params.distance_threshold = 0.25
params.distance_metric = "cosine"
params.linkage_type = "average"
params.cag_batchsize = 10000

// Statistical analysis options
params.formula = false
params.fdr_method = "fdr_bh"
params.corncob_batches = 10


// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run Golob-Minot/geneshot <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)

    
    Flow control options:
      --nopreprocess        If specified, omit the preprocessing steps (removing adapters and human sequences). Assume manifest is QCed
      --nocompostion        If specified, will skip the metaphlan2 compositional analysis steps.

    Options:
      --output              Folder to place analysis outputs (default ./results)
      -w                    Working directory. Defaults to `./work`

    For preprocessing:
      --host_index_url        URL for host genome index, defaults to current Human Genome
      --host_index            or Cached copy of the bwa indexed human genome, TGZ format
      --min_host_align_score  Minimum alignment score for human genome (default 30)
      --savereads             If provided, save the preprocessed reads to the qc/ subdirectory.

    For Allele Catalog via Assembly:
      --phred_offset        for spades. Default 33.

    For Quantification via Alignment:
      --dmnd_min_identity   Amino acid identity cutoff used to align short reads (default: 90) (DIAMOND)
      --dmnd_min_coverage   Query coverage cutoff used to align short reads (default: 80) (DIAMOND)
      --dmnd_top_pct        Keep top X% of alignments for each short read (default: 1) (DIAMOND)
      --dmnd_min_score      Minimum score for short read alignment (default: 20) (DIAMOND)
      --gencode             Genetic code used for conceptual translation (default: 11) (DIAMOND)
      --sd_mean_cutoff      Ratio of standard deviation / mean depth of sequencing used to filter genes (default: 3.0) (FAMLI)

    Batchfile:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This can be repeated. 
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
      If index reads are provided, the column titles should be 'I1' and 'I2'

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.manifest == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Make sure that --output ends with trailing "/" characters
if (!params.output.endsWith("/")){
    params.output_folder = params.output.concat("/")
} else {
    params.output_folder = params.output
}

// Import the preprocess_wf module
include { Preprocess_wf } from './modules/preprocess' params(
    manifest: params.manifest,
    host_index: params.host_index,
    host_index_url: params.host_index_url,
    min_host_align_score: params.min_host_align_score,
    savereads: params.savereads,
    output: params.output_folder,

)
include { Read_manifest } from './modules/general'

include {  Metaphlan2_wf } from './modules/composition' params(
    manifest: params.manifest,
    output: params.output_folder,
)

// Import the workflows used for assembly-based allele-catalog
include { Allele_catalog } from './modules/allele_catalog' params(
    output: params.output_folder,
    phred_offset: params.phred_offset,
)

// Import the workflow responsible for clustering alleles into 'genes'
include { Allele_clustering } from './modules/allele_clustering' params(
    output: params.output_folder,
)

include { Alignment_wf } from './modules/quantify' params (
    output: params.output_folder,
    dmnd_min_identity: params.dmnd_min_identity,
    dmnd_min_coverage: params.dmnd_min_coverage,
    dmnd_top_pct: params.dmnd_top_pct,
    dmnd_min_score: params.dmnd_min_score,
    gencode: params.gencode,
    sd_mean_cutoff: params.sd_mean_cutoff,

)
/*



// ---


/*
// Import the workflows used for annotation
include { Annotation_wf } from './modules/annotation' params(
    output_folder: output_folder,
    phred_offset: params.phred_offset,
    min_identity: params.min_identity,
    min_coverage: params.min_coverage,
    noannot: params.noannot,
    eggnog_db: params.eggnog_db,
    eggnog_dmnd: params.eggnog_dmnd,
    taxonomic_dmnd: params.taxonomic_dmnd,
    gencode: params.gencode,
)


*/



workflow {
    main:

    // ##########################
    // # PREPROCESSING          #
    // ##########################
    manifest_file = Channel.from(file(params.manifest))
    manifest_qced = Read_manifest(manifest_file)

    // Phase I: Preprocessing
    if (!params.nopreprocess) {

        // Run the entire preprocessing workflow
        Preprocess_wf(
            manifest_qced.valid_paired_indexed,
             manifest_qced.valid_paired
        )

        combined_reads_pe = Preprocess_wf.out

    } else {
        // If the user specified --nopreprocess, then just read in the manifest assuming these are already QCed and normalized.
        combined_reads_pe = manifest_qced.valid_paired.mix(manifest_qced.valid_paired_indexed)
            .map { 
                r -> [r.specimen, file(r.R1), file(r.R2)]
            }
    }

    // ##########################
    // # COMPOSITIONAL ANALYSIS #
    // ##########################
    if (!params.nocomposition) {
        Metaphlan2_wf(
            combined_reads_pe,
            Channel.from([])
        )
    }


    // #########################################
    // # ALLELE CATALOG FROM  DE NOVO ASSEMBLY #
    // #########################################

    Allele_catalog(
        combined_reads_pe,
    )

    // #########################################
    // # ALLELE CLUSTERING INTO "GENES"        #
    // #########################################

    Allele_clustering(
        Allele_catalog.out.alleles,
        Allele_catalog.out.allele_info
    )

    // ##################################
    // # ALIGNMENT-BASED QUANTIFICATION #
    // ##################################

    Alignment_wf(
        Allele_catalog.out.alleles,
        Allele_catalog.out.alleles_dmdb,
        combined_reads_pe,
    )
    Alignment_wf.out.specimen_allele_quant
    
    // */

}
