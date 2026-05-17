import rich_click as click
from dataclasses import fields
from typing import Any, Dict, Optional
import psutil



# Utility to create config dataclass from kwargs with optional namespacing
def create_config(config_class, namespace: Optional[str] = None, **kwargs):
    """
    Create a configuration object from a dataclass, using kwargs with optional namespacing.
    If namespace is provided, parameters with the namespace suffix are used preferentially.
    For example, for namespace="tc", "e_value_threshold_tc" in kwargs will override "e_value_threshold".
    1. namespaced parameter (e.g., e_value_threshold_tc) if exists
    2. general parameter (e.g., e_value_threshold)
    3. default value in dataclass
    """
    field_names = {f.name for f in fields(config_class)}
    resolved: Dict[str, Any] = {}
    for fname in field_names:
        ns_key = f"{fname}_{namespace}" if namespace else None
        if ns_key and ns_key in kwargs:
            resolved[fname] = kwargs[ns_key]
        elif fname in kwargs:
            resolved[fname] = kwargs[fname]
    return config_class.from_dict(config_class, resolved)



output_dir_option = click.option(
    '--output_dir',
    required=True,
    type=click.Path(file_okay=False, dir_okay=True),
    help='Output directory (created if needed). All run_dbcan result files are written here.',
)
threads_option = click.option(
    '--threads',
    type=int,
    default=psutil.cpu_count() or 1,
    show_default=True,
    help='Parallel worker count for supported tools (DIAMOND, pyhmmer CPU threads, etc.).',
)


class MethodsType(click.ParamType):
    """Click parameter type for comma-separated methods string."""
    name = "methods"
    
    def convert(self, value, param, ctx):
        if value is None or not value.strip():
            return ('diamond', 'hmm', 'dbCANsub')
        
        valid_methods = {'diamond', 'hmm', 'dbCANsub'}
        methods = [m.strip().lower() for m in value.split(',') if m.strip()]
        
        if not methods:
            return ('diamond', 'hmm', 'dbCANsub')
        
        normalized = []
        for m in methods:
            if m == 'dbcansub':
                normalized.append('dbCANsub')
            elif m in valid_methods:
                normalized.append(m)
            else:
                self.fail(
                    f"Invalid method '{m}'. Valid options are: diamond, hmm, dbCANsub",
                    param=param,
                    ctx=ctx
                )
        
        return tuple(normalized)


methods_option = click.option(
    '--methods',
    default='diamond,hmm,dbCANsub',
    type=MethodsType(),
    help=(
        'CAZyme annotation modules for CAZyme_annotation / easy_* step 1, comma-separated. '
        'Choices: diamond (CAZy DIAMOND), hmm (dbCAN HMM), dbCANsub (dbCAN-sub HMM). '
        'Example: --methods diamond,hmm'
    ),
    show_default=True,
)



# Define group options
def general_options(func):
    func = click.option(
        '--input_raw_data',
        required=True,
        help='Path to input sequences (FASTA/FASTA.gz): nucleotide for prok/meta modes, proteins when --mode=protein.',
    )(func)
    func = output_dir_option(func)
    func = click.option(
        '--mode',
        default='prok',
        required=True,
        show_default=True,
        help='Input type: prok (prokaryote DNA), meta (metagenome DNA), protein (amino-acid FASTA).',
    )(func)

    return func

def database_options(func):
    """Common database options shared by most commands (no CGC switch)."""
    func = click.option(
        '--db_dir',
        required=True,
        type=click.Path(file_okay=False, dir_okay=True),
        help='Directory containing dbCAN database files (HMM, DIAMOND, etc.).',
    )(func)
    return func

def database_download_options(func):
    """Database downloader-only options (adds CGC switch)."""
    func = database_options(func)
    func = click.option(
        '--cgc/--no-cgc',
        is_flag=True,
        default=True,
        show_default=True,
        help='With --cgc (default): download CGC-related DB assets. With --no-cgc: skip them (database subcommand only).',
    )(func)
    func = click.option(
        '--aws_s3',
        is_flag=True,
        default=False,
        show_default=True,
        help='Download from the pinned AWS S3 release; omit for HTTP db_current (moving snapshot).',
    )(func)
    func = click.option(
        '--timeout',
        type=int,
        default=30,
        show_default=True,
        help='HTTP(S) request timeout in seconds per download attempt.',
    )(func)
    func = click.option(
        '--retries',
        type=int,
        default=3,
        show_default=True,
        help='Retries for transient HTTP/S3 download failures.',
    )(func)
    func = click.option(
        '--resume/--no-resume',
        is_flag=True,
        default=True,
        show_default=True,
        help='Resume partial downloads when supported (--no-resume: always fetch from scratch).',
    )(func)
    func = click.option(
        '--no-overwrite',
        is_flag=True,
        default=False,
        show_default=True,
        help='If set, skip downloading files that already exist in db_dir.',
    )(func)
    func = click.option(
        '--verify-ssl/--no-verify-ssl',
        is_flag=True,
        default=True,
        show_default=True,
        help='Verify TLS certificates for HTTPS downloads (--no-verify-ssl: insecure, not recommended).',
    )(func)
    return func

def diamond_options(func):
    func = click.option(
        '--e_value_threshold',
        type=float,
        help='Maximum E-value for DIAMOND CAZy hits (float; stricter = smaller).',
        default=1e-102,
        show_default=True,
    )(func)
    func = click.option(
        '--verbose_option',
        is_flag=True,
        default=False,
        show_default=True,
        help='Pass DIAMOND --verbose (more DIAMOND stderr output).',
    )(func)
    return func

def diamond_tc_options(func):
    func = click.option(
        '--e_value_threshold_tc',
        type=float,
        help='Maximum E-value for TC (transporter) DIAMOND vs TCDB.',
        default=1e-4,
        show_default=True,
    )(func)
    func = click.option(
        '--coverage_threshold_tc',
        type=float,
        help='DIAMOND --query-cover for TCDB search, as percent of query length (0–100; default 35 = 35%%).',
        default=35,
        show_default=True,
    )(func)
    return func


def diamond_tf_options(func):
    func = click.option(
        '--e_value_threshold_tf_diamond',
        type=float,
        help='Maximum E-value for TF DIAMOND (prokaryotic TFDB path).',
        default=1e-4,
        show_default=True,
    )(func)
    func = click.option(
        '--coverage_threshold_tf_diamond',
        type=float,
        help='DIAMOND --query-cover for TF DIAMOND: minimum query coverage in percent (0–100; same semantics as --coverage_threshold_tc, default 35%%).',
        default=35,
        show_default=True,
    )(func)
    func = click.option(
        '--prokaryotic/--no-prokaryotic',
        is_flag=True,
        default=True,
        show_default=True,
        help='Run prokaryotic TF DIAMOND step (--no-prokaryotic: skip it; use fungi TF HMM with --fungi).',
    )(func)
    return func

def pyhmmer_dbcan_options(func):
    _mem_note = (
        " Applies to pyhmmer TF/STP in gff_process (easy_* step 2) as well—no separate TF/STP memory options."
    )
    func = click.option(
        '--e_value_threshold_dbcan',
        type=float,
        help='Maximum domain-independent E-value to keep a dbCAN HMM hit (pyhmmer; CAZyme_annotation / easy_* step 1 only).',
        default=1e-15,
        show_default=True,
    )(func)
    func = click.option(
        '--coverage_threshold_dbcan',
        type=float,
        help='Minimum hit coverage (0.0–1.0 fraction of query length) for dbCAN HMM (step 1 only).',
        default=0.35,
        show_default=True,
    )(func)
    func = click.option(
        '--csv_buffer_size',
        type=int,
        default=5000,
        show_default=True,
        help=(
            'Buffer this many HMM hit rows before flushing to TSV (larger = fewer writes, slightly more RAM).'
            + _mem_note
        ),
    )(func)
    func = click.option(
        '--batch_size',
        type=int,
        default=None,
        show_default=True,
        help=(
            'Protein sequences per pyhmmer batch for dbCAN HMM; omit for automatic sizing from free RAM.'
            + _mem_note
        ),
    )(func)
    func = click.option(
        '--max_memory_usage',
        type=float,
        default=0.8,
        show_default=True,
        help=(
            'If system memory usage ratio exceeds this (0.0–1.0), emit warnings and tighten batching.'
            + _mem_note
        ),
    )(func)
    func = click.option(
        '--memory_safety_factor',
        type=float,
        default=0.5,
        show_default=True,
        help=(
            'Fraction of available RAM used when estimating automatic batch_size (0.0–1.0; lower = smaller batches).'
            + _mem_note
        ),
    )(func)
    func = click.option(
        '--max_retries',
        type=int,
        default=3,
        show_default=True,
        help='Retries after MemoryError during pyhmmer HMM search, halving batch each retry.' + _mem_note,
    )(func)
    func = click.option(
        '--enable_memory_monitoring/--no-enable_memory_monitoring',
        is_flag=True,
        default=True,
        show_default=True,
        help='Track RAM and adapt pyhmmer batching; use --no-enable_memory_monitoring to disable.' + _mem_note,
    )(func)
    func = click.option(
        '--large/--no-large',
        'large_mode',
        is_flag=True,
        default=False,
        show_default=True,
        help=(
            'Streaming-safe pyhmmer mode for huge inputs (less preload, lower OOM risk).'
            + _mem_note
        ),
    )(func)
    func = click.option(
        '--large_input_threshold_mb',
        type=int,
        default=5000,
        show_default=True,
        help=(
            'If the input FASTA for the HMM step exceeds this size (MB), enable large mode automatically.'
            + _mem_note
        ),
    )(func)
    return func

def dbcansub_options(func):
    func = click.option(
        '--e_value_threshold_dbsub',
        type=float,
        help='Maximum E-value to keep a dbCAN-sub HMM hit (pyhmmer; step 1 only).',
        default=1e-15,
        show_default=True,
    )(func)
    func = click.option(
        '--coverage_threshold_dbsub',
        type=float,
        help='Minimum hit coverage (0.0–1.0) for dbCAN-sub HMM (step 1 only).',
        default=0.35,
        show_default=True,
    )(func)
    func = click.option(
        '--csv_buffer_size_dbsub',
        type=int,
        default=5000,
        show_default=True,
        help='dbCAN-sub: buffer this many HMM hit rows before flushing to disk.',
    )(func)
    func = click.option(
        '--batch_size_dbsub',
        type=int,
        default=None,
        show_default=True,
        help='dbCAN-sub: sequences per pyhmmer batch; omit for automatic sizing (independent from --batch_size).',
    )(func)
    func = click.option(
        '--max_memory_usage_dbsub',
        type=float,
        default=0.8,
        show_default=True,
        help='dbCAN-sub: system memory usage ratio threshold (0.0–1.0) for warnings / throttling.',
    )(func)
    func = click.option(
        '--memory_safety_factor_dbsub',
        type=float,
        default=0.5,
        show_default=True,
        help='dbCAN-sub: RAM fraction used in automatic batch_size estimate (0.0–1.0).',
    )(func)
    func = click.option(
        '--max_retries_dbsub',
        type=int,
        default=3,
        show_default=True,
        help='dbCAN-sub: max pyhmmer retries after MemoryError.',
    )(func)
    func = click.option(
        '--enable_memory_monitoring_dbsub/--no-enable_memory_monitoring_dbsub',
        is_flag=True,
        default=True,
        show_default=True,
        help='dbCAN-sub: enable RAM monitoring and adaptive batching.',
    )(func)
    func = click.option(
        '--large_dbsub/--no-large_dbsub',
        'large_mode_dbsub',
        is_flag=True,
        default=False,
        show_default=True,
        help='dbCAN-sub: force streaming-safe pyhmmer mode.',
    )(func)
    func = click.option(
        '--large_input_threshold_mb_dbsub',
        type=int,
        default=5000,
        show_default=True,
        help='dbCAN-sub: auto large mode when input FASTA exceeds this size (MB).',
    )(func)
    return func

def pyhmmer_tf(func):
    func = click.option(
        '--e_value_threshold_tf',
        type=float,
        help='Maximum E-value for TF pyhmmer hits (gff_process / easy_* step 2; fungi mode only).',
        default=1e-4,
        show_default=True,
    )(func)
    func = click.option(
        '--coverage_threshold_tf',
        type=float,
        help='Minimum hit coverage (0.0–1.0) for TF pyhmmer (step 2; requires --fungi).',
        default=0.35,
        show_default=True,
    )(func)
    func = click.option(
        '--fungi/--no-fungi',
        is_flag=True,
        default=False,
        show_default=True,
        help='Run TF HMM search (--fungi); default skips TF pyhmmer (prokaryotes use TF DIAMOND instead).',
    )(func)
    return func

def pyhmmer_stp(func):
    func = click.option(
        '--e_value_threshold_stp',
        type=float,
        help='Maximum E-value for STP pyhmmer hits (gff_process / easy_* step 2).',
        default=1e-4,
        show_default=True,
    )(func)
    func = click.option(
        '--coverage_threshold_stp',
        type=float,
        help='Minimum hit coverage (0.0–1.0) for STP pyhmmer (step 2).',
        default=0.35,
        show_default=True,
    )(func)
    return func

def diamond_sulfatase_options(func):
    func = click.option(
        '--e_value_threshold_sulfatase',
        type=float,
        help='Maximum E-value for Sulfatase DIAMOND.',
        default=1e-4,
        show_default=True,
    )(func)
    func = click.option(
        '--coverage_threshold_sulfatase',
        type=float,
        help='DIAMOND --query-cover for Sulfatase: minimum query coverage in percent (0–100).',
        default=35,
        show_default=True,
    )(func)
    return func

def diamond_peptidase_options(func):
    func = click.option(
        '--e_value_threshold_peptidase',
        type=float,
        help='Maximum E-value for Peptidase DIAMOND.',
        default=1e-4,
        show_default=True,
    )(func)
    func = click.option(
        '--coverage_threshold_peptidase',
        type=float,
        help='DIAMOND --query-cover for Peptidase: minimum query coverage in percent (0–100).',
        default=35,
        show_default=True,
    )(func)
    return func

def pyhmmer_pfam(func):
    func = click.option(
        '--run_pfam',
        is_flag=True,
        default=False,
        show_default=True,
        help='Run Pfam pyhmmer on null genes (Pfam_null_cgc / CGC null annotation).',
    )(func)
    func = click.option(
        '--e_value_threshold_pfam',
        type=float,
        help='Maximum domain-independent E-value for Pfam hits.',
        default=1e-4,
        show_default=True,
    )(func)
    func = click.option(
        '--coverage_threshold_pfam',
        type=float,
        help='Minimum HMM alignment coverage (0.0–1.0 fraction of HMM length) for Pfam hits.',
        default=0.35,
        show_default=True,
    )(func)
    func = click.option(
        '--null_from_gff',
        is_flag=True,
        default=False,
        show_default=True,
        help='Extract null-gene proteins from cgc.gff instead of cgc_standard_out.tsv.',
    )(func)
    return func

def cgc_gff_option(func):
    func = click.option(
        '--input_gff',
        required=False,
        default=None,
        type=click.Path(),
        help='Annotation GFF. Required when --mode=protein; otherwise defaults to <output_dir>/uniInput.gff.',
    )(func)
    func = click.option(
        '--gff_type',
        required=False,
        default=None,
        type=str,
        help="GFF dialect label (e.g. prodigal, NCBI_prok). When --mode!=protein defaults to 'prodigal' if unset.",
    )(func)
    func = click.option(
        '--prodigal-gff-streaming',
        'prodigal_gff_streaming',
        type=click.Choice(['auto', 'on', 'off'], case_sensitive=False),
        default='auto',
        show_default=True,
        help=(
            "Prodigal GFF annotation: 'auto' stream-annotate when input exceeds threshold (faster for huge GFF); "
            "'on' always use streaming; 'off' always use BCBio GFF.parse."
        ),
    )(func)
    func = click.option(
        '--prodigal-streaming-threshold-mb',
        'prodigal_streaming_threshold_mb',
        type=int,
        default=50,
        show_default=True,
        help='When --prodigal-gff-streaming=auto, stream if GFF size exceeds this (decimal MB). Use 0 to stream any non-empty file.',
    )(func)
    return func

def cgc_options(func):
    func = click.option(
        '--additional_genes',
        multiple=True,
        default=["TC"],
        show_default=True,
        help='Gene class tags required alongside CAZyme for CGC signatures. Repeat: --additional_genes TC --additional_genes TF. Choices include TC, TF, STP.',
    )(func)
    func = click.option(
        '--additional_logic',
        type=click.Choice(['all', 'any']),
        default='all',
        show_default=True,
        help="How to combine --additional_genes: 'all' = every listed class must be present; 'any' = at least --additional_min_categories classes.",
    )(func)
    func = click.option(
        '--additional_min_categories',
        type=int,
        default=1,
        show_default=True,
        help='When --additional_logic=any, minimum number of distinct additional gene classes that must match.',
    )(func)
    func = click.option(
        '--num_null_gene',
        type=int,
        default=2,
        show_default=True,
        help='Max number of intervening non-signature (null) genes allowed between core CGC genes.',
    )(func)
    func = click.option(
        '--base_pair_distance',
        type=int,
        default=15000,
        show_default=True,
        help='Max distance (bp) between CGC signature genes when --use_distance is enabled.',
    )(func)
    func = click.option(
        '--use_null_genes/--no-use_null_genes',
        is_flag=True,
        default=True,
        show_default=True,
        help='Allow null genes between CGC signatures (--no-use_null_genes: disable).',
    )(func)
    func = click.option(
        '--use_distance',
        is_flag=True,
        default=False,
        show_default=True,
        help='Also require signature genes to fall within --base_pair_distance bp.',
    )(func)

    # extended options
    func = click.option(
        '--extend_mode',
        type=click.Choice(['none', 'bp', 'gene']),
        default='none',
        show_default=True,
        help="After CGC detection, extend cluster bounds: 'bp' uses --extend_bp; 'gene' uses --extend_gene_count; 'none' disables.",
    )(func)
    func = click.option(
        '--extend_bp',
        type=int,
        default=0,
        show_default=True,
        help='With --extend_mode=bp, extend each side by this many base pairs.',
    )(func)
    func = click.option(
        '--extend_gene_count',
        type=int,
        default=0,
        show_default=True,
        help='With --extend_mode=gene, extend each side by this many flanking genes.',
    )(func)

    # newly added parameters
    func = click.option(
        '--min_core_cazyme',
        type=int,
        default=1,
        show_default=True,
        help='Minimum core CAZyme count required to retain a CGC.',
    )(func)
    func = click.option(
        '--min_cluster_genes',
        type=int,
        default=2,
        show_default=True,
        help='Minimum total genes in a CGC locus.',
    )(func)
    func = click.option(
        '--feature_type',
        'feature_types',
        multiple=True,
        default=["CDS"],
        show_default=True,
        help='GFF feature types to parse (repeat option for multiple values, e.g. --feature_type CDS --feature_type gene).',
    )(func)
    return func

def cgc_substrate_base_options(func):
    """Options shared by substrate_prediction and easy_substrate (step 4)."""
    func = general_options(func)
    func = click.option(
        '--pul',
        required=False,
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help='Path to dbCAN-PUL PUL.faa when you need an explicit PUL database file.',
    )(func)
    func = click.option(
        '-o',
        '--out',
        default='substrate.out',
        show_default=True,
        help='Legacy filename hint for standalone substrate tools (results still go under --output_dir).',
    )(func)
    func = click.option(
        '-w',
        '--workdir',
        default='.',
        show_default=True,
        type=str,
        help='Working directory for legacy substrate scripts (prefer --output_dir for run_dbcan).',
    )(func)
    func = click.option(
        '-rerun',
        '--rerun',
        default=False,
        type=bool,
        show_default=True,
        help='Re-run substrate prediction (pass true/false explicitly; rarely needed).',
    )(func)
    func = click.option(
        '-env',
        '--env',
        default='local',
        show_default=True,
        type=str,
        help='Execution environment label for external wrappers (usually keep local).',
    )(func)
    func = click.option(
        '-odbcan_sub',
        '--odbcan_sub',
        type=bool,
        default=None,
        help='If set to true/false, force exporting extra dbCAN-sub tables; omit to use package default.',
    )(func)
    func = click.option(
        '-odbcanpul',
        '--odbcanpul',
        default=True,
        type=bool,
        show_default=True,
        help='Whether to export dbCAN-PUL homology tables (pass true/false).',
    )(func)
    func = click.option(
        '--db_dir',
        default='./dbCAN_databases',
        required=True,
        show_default=True,
        type=click.Path(file_okay=False, dir_okay=True),
        help='Database directory containing PUL/DIAMOND assets for substrate prediction.',
    )(func)
    return func

def cgc_substrate_homology_params_options(func):
    """Homology (dbCAN-PUL) filters for substrate_prediction."""
    func = click.option(
        '-upghn',
        '--uniq_pul_gene_hit_num',
        default=2,
        type=int,
        show_default=True,
        help='Minimum unique PUL genes hit by BLAST for a valid PUL–CGC link.',
    )(func)
    func = click.option(
        '-uqcgn',
        '--uniq_query_cgc_gene_num',
        default=2,
        type=int,
        show_default=True,
        help='Minimum unique CGC genes participating in a PUL–CGC link.',
    )(func)
    func = click.option(
        '-cpn',
        '--CAZyme_pair_num',
        default=1,
        type=int,
        show_default=True,
        help='Minimum CAZyme–CAZyme pairs required inside the CGC for homology scoring.',
    )(func)
    func = click.option(
        '-tpn',
        '--total_pair_num',
        default=2,
        type=int,
        show_default=True,
        help='Minimum total informative gene pairs (CAZyme + accessory) for a link.',
    )(func)
    func = click.option(
        '-ept',
        '--extra_pair_type',
        default=None,
        type=str,
        help='Optional comma-separated accessory gene pair types (advanced PUL matching).',
    )(func)
    func = click.option(
        '-eptn',
        '--extra_pair_type_num',
        default='0',
        show_default=True,
        type=str,
        help='Comma-separated counts matching --extra_pair_type entries (same order).',
    )(func)
    func = click.option(
        '-iden',
        '--identity_cutoff',
        default=0.0,
        type=float,
        show_default=True,
        help='Minimum BLAST identity (0.0–1.0) for PUL–CGC homology hits.',
    )(func)
    func = click.option(
        '-cov',
        '--coverage_cutoff',
        default=0.0,
        type=float,
        show_default=True,
        help='Minimum BLAST query coverage (0.0–1.0) for PUL–CGC homology hits.',
    )(func)
    func = click.option(
        '-bsc',
        '--bitscore_cutoff',
        default=50.0,
        type=float,
        show_default=True,
        help='Minimum BLAST bit score for PUL–CGC homology hits.',
    )(func)
    func = click.option(
        '-evalue',
        '--evalue_cutoff',
        default=0.01,
        type=float,
        show_default=True,
        help='Maximum BLAST E-value for PUL–CGC homology hits.',
    )(func)
    return func

def cgc_substrate_dbcan_sub_param_options(func):
    """dbCAN-sub evidence filters inside substrate_prediction."""
    func = click.option(
        '-hmmcov',
        '--hmmcov',
        default=0.0,
        type=float,
        show_default=True,
        help='Minimum dbCAN-sub HMM coverage (0.0–1.0) when scoring substrate evidence.',
    )(func)
    func = click.option(
        '-hmmevalue',
        '--hmmevalue',
        default=0.01,
        type=float,
        show_default=True,
        help='Maximum dbCAN-sub HMM E-value allowed for substrate evidence.',
    )(func)
    func = click.option(
        '-ndsc',
        '--num_of_domains_substrate_cutoff',
        default=2,
        type=int,
        show_default=True,
        help='Minimum distinct substrate-associated domains required per prediction.',
    )(func)
    func = click.option(
        '-npsc',
        '--num_of_protein_substrate_cutoff',
        default=2,
        type=int,
        show_default=True,
        help='Minimum proteins supporting a substrate call within a CGC.',
    )(func)
    func = click.option(
        '-subs',
        '--substrate_scors',
        default=2.0,
        type=float,
        show_default=True,
        help='Minimum aggregated dbCAN-sub score (field name substrate_scors) to accept a substrate assignment.',
    )(func)
    return func

def cgc_sub_options(func):
    """total option for cgc substrate prediction"""
    func = cgc_substrate_base_options(func)
    func = cgc_substrate_homology_params_options(func)
    func = cgc_substrate_dbcan_sub_param_options(func)
    return func

def syn_plot_options(func):
    func = click.option(
        '--db_dir',
        required=True,
        type=click.Path(file_okay=False, dir_okay=True),
        help='Database directory for synteny / PUL reference data (syntenic_plot).',
    )(func)
    return func

def cgc_circle_plot_options(func):
    func = output_dir_option(func)
    return func

def topology_annotation_options(func):
    func = click.option(
        '--run_signalp/--no-run_signalp',
        default=False,
        show_default=True,
        help='Run SignalP 6 (BioLib) on translated proteins and append peptide signal columns to overview.',
    )(func)
    func = click.option(
        '--run_deeptmhmm/--no-run_deeptmhmm',
        default=False,
        show_default=True,
        help='Run a user-installed DeepTMHMM predict.py and append transmembrane predictions to overview.',
    )(func)
    func = click.option(
        '--deeptmhmm_dir',
        type=click.Path(exists=False, file_okay=False, dir_okay=True),
        default=None,
        help='Directory that contains DeepTMHMM predict.py (only used with --run_deeptmhmm).',
    )(func)
    func = click.option(
        '--deeptmhmm_python',
        default='python',
        show_default=True,
        help='Python interpreter used to launch DeepTMHMM predict.py.',
    )(func)
    func = click.option(
        '--signalp_org',
        default='other',
        type=click.Choice(['other', 'euk']),
        show_default=True,
        help='SignalP organism class: other (bacteria/archaea) or euk (eukaryotes).',
    )(func)
    func = click.option(
        '--force_topology/--no-force_topology',
        default=False,
        show_default=True,
        help='Recompute topology columns even when overview already contains predictions.',
    )(func)
    return func

def expression_options(func):
    """Options for run_dbcan expression."""
    func = click.option(
        '-i', '--input_dir',
        required=True,
        type=click.Path(exists=True, file_okay=False, dir_okay=True),
        help='run_dbcan output directory (contains cgc_standard_out.tsv, overview.tsv).',
    )(func)
    func = click.option(
        '--samplesheet',
        required=True,
        type=click.Path(exists=True, dir_okay=False),
        help='TSV samplesheet: sample_id, condition, bam and/or r1/r2 columns.',
    )(func)
    func = click.option(
        '--reference-fasta',
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help='Reference assembly FASTA for BWA (required when aligning from FASTQ).',
    )(func)
    func = click.option(
        '--gff',
        required=False,
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help=(
            'Gene GFF for read counting. Default: {input_dir}/cgc.gff if present, '
            'else *.fix.gff or uniInput.gff. Supports protein_id= (cgc.gff) and ID= (fix.gff).'
        ),
    )(func)
    func = click.option(
        '--output-dir',
        'output_dir',
        default=None,
        show_default='{input_dir}.expression',
        help='Expression analysis output directory.',
    )(func)
    func = click.option(
        '--design',
        default='~condition',
        show_default=True,
        help='PyDESeq2 design formula.',
    )(func)
    func = click.option(
        '--abundance',
        default='TPM',
        type=click.Choice(['FPKM', 'RPM', 'TPM']),
        show_default=True,
        help='Normalization for abundance tables.',
    )(func)
    func = click.option(
        '--alpha',
        default=0.05,
        type=float,
        show_default=True,
        help='Adjusted p-value threshold for DEG calling.',
    )(func)
    func = click.option(
        '--lfc-threshold',
        default=1.0,
        type=float,
        show_default=True,
        help='Minimum |log2 fold change| for DEG calling.',
    )(func)
    func = click.option(
        '--cgc-de-rule',
        default='any',
        type=click.Choice(['any', 'majority', 'all']),
        show_default=True,
        help='Rule to label CGC as DE based on gene-level DEGs.',
    )(func)
    func = click.option(
        '--all-genes',
        is_flag=True,
        default=False,
        help='Include all GFF genes in count matrix (default: CAZyme/CGC genes only).',
    )(func)
    func = click.option(
        '--skip-alignment',
        is_flag=True,
        default=False,
        help='Use BAM paths from samplesheet only (no BWA).',
    )(func)
    func = click.option(
        '--skip-deseq2',
        is_flag=True,
        default=False,
        help='Skip differential expression; only count and abundance.',
    )(func)
    func = click.option(
        '--run-plots',
        is_flag=True,
        default=False,
        help='Generate expression and DEG plots after analysis.',
    )(func)
    func = click.option(
        '--run-plots-only-de',
        is_flag=True,
        default=False,
        help='With --run-plots, only plot DE CGCs (not all CGCs).',
    )(func)
    func = click.option(
        '--also-heatmap',
        is_flag=True,
        default=False,
        help='With --run-plots, also output CGC expression heatmap.',
    )(func)
    func = click.option(
        '--max-cgc',
        default=500,
        type=int,
        show_default=True,
        help='Max CGCs to plot without --force.',
    )(func)
    func = click.option(
        '--force',
        is_flag=True,
        default=False,
        help='Allow plotting more than --max-cgc CGCs.',
    )(func)
    func = click.option(
        '--overlap-base-ratio',
        default=0.2,
        type=float,
        show_default=True,
        help='Minimum read overlap fraction for cal_coverage.',
    )(func)
    func = click.option(
        '--mapping-quality',
        default=30,
        type=int,
        show_default=True,
        help='Minimum MAPQ for cal_coverage.',
    )(func)
    func = click.option(
        '--identity',
        default=0.98,
        type=float,
        show_default=True,
        help='Minimum alignment identity for cal_coverage.',
    )(func)
    func = click.option(
        '--hifi',
        is_flag=True,
        default=False,
        help='HiFi read mode for cal_coverage overlap calculation.',
    )(func)
    func = threads_option(func)
    return func


def logging_options(func):
    """Global logging options for all commands."""
    func = click.option(
        '--log-level',
        type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']),
        default='WARNING',
        show_default=True,
        help=(
            'Python logging level for run_dbcan and nested subcommands (including easy_CGC / easy_substrate steps). '
            'DEBUG shows diagnostic messages where implemented.'
        ),
    )(func)
    func = click.option(
        '--log-file',
        type=click.Path(),
        default=None,
        help='Log file path (truncates each run); mirrors console when set.',
    )(func)
    func = click.option(
        '--verbose',
        '-v',
        is_flag=True,
        default=False,
        show_default=True,
        help='Shortcut for --log-level DEBUG (overrides --log-level if both are passed).',
    )(func)
    return func






