#!/usr/bin/env bash
#SBATCH --job-name=chipseq_nf
#SBATCH --partition=batch
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G
#SBATCH --time=120:00:00
#SBATCH --qos=long_jobs
#SBATCH --output=jobs/%x_%j.log
#SBATCH --exclude=cnode-07-[5-29],cnode-08-[1-35]
#
# Nextflow DRIVER for nf-core/chipseq 2.1.0 on ARC. Each pipeline process is
# its own SLURM job (conf/arc.config), so the driver needs only 2 CPU / 10 GB.
# Pattern copied from the lab's working ScaleTagToolkit driver: module
# nextflow/25.10.0 with the engine pinned, --qos=long_jobs (lifts the 36 h
# batch cap for the driver only), no --account, work/ on MaxsonBraunScratch.
#
# Run from the root of this folder:
#   mkdir -p jobs
#   PREVIEW=1 bash run_chipseq_arc.sbatch        # login node: checks + nextflow -preview
#   sbatch run_chipseq_arc.sbatch                # real run
#
# Overrides (environment variables, e.g. SAMPLESHEET=x.csv sbatch --export=ALL run_chipseq_arc.sbatch):
#   SAMPLESHEET  default samplesheet.csv
#   PARAMS       default params/hg38.yaml
#   OUTDIR       default results
#   WORKDIR      default /home/groups/MaxsonBraunScratch/$USER/nxf-work/<this folder name>
#   NXF_SINGULARITY_CACHEDIR  default /home/groups/MaxsonLab/software/singularity-containers/nf-core
#   MIN_FREE_GB  default 500 (free space required under WORKDIR)
# Extra arguments are passed to nextflow run (e.g. --skip_preseq).

set -u

# PIPELINE may be overridden with a local copy (e.g. from `nf-core pipelines
# download`) when compute nodes cannot reach GitHub; -r is then not used.
PIPELINE="${PIPELINE:-nf-core/chipseq}"
PIPELINE_REV="2.1.0"
NXF_MODULE="nextflow/25.10.0"
NXF_PIN="24.10.0"

# sbatch runs a COPY of this script from the spool dir: locate files via the
# submit directory, never via the script's own path.
KIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
cd "$KIT_DIR" || { echo "FATAL: cannot cd to $KIT_DIR" >&2; exit 1; }

SAMPLESHEET="${SAMPLESHEET:-samplesheet.csv}"
PARAMS="${PARAMS:-params/hg38.yaml}"
OUTDIR="${OUTDIR:-results}"
RUN_NAME="${RUN_NAME:-$(basename "$KIT_DIR")}"
WORKDIR="${WORKDIR:-/home/groups/MaxsonBraunScratch/$USER/nxf-work/$RUN_NAME}"
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-/home/groups/MaxsonLab/software/singularity-containers/nf-core}"
export NXF_APPTAINER_CACHEDIR="$NXF_SINGULARITY_CACHEDIR"
MIN_FREE_GB="${MIN_FREE_GB:-500}"
PREVIEW="${PREVIEW:-0}"

if [ -z "${SLURM_JOB_ID:-}" ] && [ "$PREVIEW" != "1" ]; then
    echo "Refusing to run the driver on a login node. Use: sbatch run_chipseq_arc.sbatch   (or PREVIEW=1 bash run_chipseq_arc.sbatch)" >&2
    exit 2
fi
mkdir -p jobs

fail=0
err() { echo "ERROR: $*" >&2; fail=1; }

echo "kit dir     : $KIT_DIR"
echo "samplesheet : $SAMPLESHEET"
echo "params      : $PARAMS"
echo "outdir      : $OUTDIR"
echo "work dir    : $WORKDIR"
echo "image cache : $NXF_SINGULARITY_CACHEDIR"

for f in conf/arc.config "$PARAMS" "$SAMPLESHEET"; do
    [ -r "$f" ] || err "cannot read $f"
done

# --- samplesheet --------------------------------------------------------------
EXPECT="sample,fastq_1,fastq_2,replicate,antibody,control,control_replicate"
if [ -r "$SAMPLESHEET" ]; then
    hdr="$(head -1 "$SAMPLESHEET" | tr -d '\r')"
    [ "$hdr" = "$EXPECT" ] || err "samplesheet header is '$hdr', expected '$EXPECT'"
    n_rows=$(awk 'NR > 1 && NF' "$SAMPLESHEET" | wc -l | tr -d ' ')
    echo "rows        : $n_rows"
    awk -F, 'NR > 1 && NF { print $2; if ($3 != "") print $3 }' "$SAMPLESHEET" | tr -d '\r' > "jobs/.fastqs.$$"
    n_missing=0
    while IFS= read -r fq; do
        case "$fq" in /*) ;; *) err "FASTQ path is not absolute: $fq"; continue ;; esac
        [ -r "$fq" ] || { err "FASTQ not readable: $fq"; n_missing=$((n_missing + 1)); }
    done < "jobs/.fastqs.$$"
    echo "fastq files : $(wc -l < "jobs/.fastqs.$$" | tr -d ' ') listed, $n_missing unreadable"
    rm -f "jobs/.fastqs.$$"
    # every control named must exist as a sample with the named replicate
    bad_ctrl=$(awk -F, 'NR > 1 && NF { gsub(/\r/, ""); have[$1 "," $4] = 1; if ($6 != "") want[$6 "," $7] = NR }
        END { for (k in want) if (!(k in have)) print "line " want[k] ": control " k }' "$SAMPLESHEET")
    [ -z "$bad_ctrl" ] || err "control not present as sample,replicate: $bad_ctrl"
fi

# --- reference files named in the params file ---------------------------------
if [ -r "$PARAMS" ]; then
    for key in fasta gtf blacklist; do
        p=$(awk -v k="$key:" '$1 == k { print $2 }' "$PARAMS")
        [ -n "$p" ] || { err "$key not set in $PARAMS"; continue; }
        [ -r "$p" ] || err "$key not readable: $p"
    done
    bt2=$(awk '$1 == "bowtie2_index:" { print $2 }' "$PARAMS")
    if [ -d "$bt2" ]; then
        n_idx=$(ls "$bt2" | grep -c '\.rev\.1\.bt2l\{0,1\}$')
        [ "$n_idx" -eq 1 ] || err "bowtie2_index dir $bt2 must hold exactly one index (found $n_idx); run scripts/prepare_reference.sbatch"
    else
        err "bowtie2_index is not a directory: $bt2 (run scripts/prepare_reference.sbatch)"
    fi
fi

# --- work dir space, image cache ---------------------------------------------
mkdir -p "$WORKDIR" 2>/dev/null || err "cannot create $WORKDIR"
if [ -d "$WORKDIR" ]; then
    free_gb=$(df -Pk "$WORKDIR" | awk 'NR == 2 { print int($4 / 1048576) }')
    echo "free space  : ${free_gb} GB under work dir (need >= $MIN_FREE_GB)"
    [ "${free_gb:-0}" -ge "$MIN_FREE_GB" ] || err "not enough space for work/ (do NOT move it to /home/exacloud/gscratch, which is campus-shared)"
fi
mkdir -p "$NXF_SINGULARITY_CACHEDIR" 2>/dev/null
touch "$NXF_SINGULARITY_CACHEDIR/.write_test" 2>/dev/null && rm -f "$NXF_SINGULARITY_CACHEDIR/.write_test" \
    || err "image cache not writable: $NXF_SINGULARITY_CACHEDIR (set NXF_SINGULARITY_CACHEDIR)"

# --- Nextflow ----------------------------------------------------------------
# `module` is a shell function; a batch shell sources no rc files.
if ! command -v module >/dev/null 2>&1; then
    for f in /etc/profile.d/modules.sh /etc/profile.d/lmod.sh /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash; do
        [ -r "$f" ] && . "$f" && break
    done
fi
if command -v module >/dev/null 2>&1; then
    module load "$NXF_MODULE" || err "cannot load $NXF_MODULE"
else
    err "module system unavailable in this shell"
fi
if command -v nextflow >/dev/null 2>&1; then
    export NXF_VER="$NXF_PIN"
    nextflow -version 2>&1 | grep -q "version $NXF_PIN" || err "nextflow engine is not the pinned $NXF_PIN"
else
    err "nextflow not on PATH after module load"
fi

if [ "$fail" -ne 0 ]; then
    echo "preflight failed; nothing launched" >&2
    exit 1
fi

MODE=""
[ "$PREVIEW" = "1" ] && MODE="-preview"
REV_ARGS="-r $PIPELINE_REV"
if [ -d "$PIPELINE" ]; then
    REV_ARGS=""
    grep -q "version *= *'$PIPELINE_REV'" "$PIPELINE/nextflow.config" 2>/dev/null \
        || { echo "FATAL: local pipeline $PIPELINE is not version $PIPELINE_REV" >&2; exit 1; }
fi
echo "[ $(date) ] nextflow -c conf/arc.config run $PIPELINE $REV_ARGS -profile singularity $MODE"
# -c is a CORE option and belongs BEFORE `run`.
nextflow -c conf/arc.config run "$PIPELINE" $REV_ARGS \
    -profile singularity \
    -params-file "$PARAMS" \
    --input "$SAMPLESHEET" \
    --outdir "$OUTDIR" \
    -work-dir "$WORKDIR" \
    -resume $MODE "$@"
status=$?
echo "[ $(date) ] nextflow exit status $status"
if [ "$status" -eq 0 ] && [ "$PREVIEW" != "1" ]; then
    echo "QC report : $OUTDIR/multiqc/"
    echo "peaks     : $OUTDIR/bowtie2/merged_library/macs3/"
    echo "work dir can be removed once results are checked: $WORKDIR"
fi
exit "$status"
