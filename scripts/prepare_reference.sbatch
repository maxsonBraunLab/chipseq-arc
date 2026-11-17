#!/usr/bin/env bash
#SBATCH --job-name=chipseq_ref_prep
#SBATCH --partition=batch
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=jobs/%x_%j.log
#
# ONE-TIME: prepare the GRCh38 reference folder that params/hg38.yaml points at.
#   1. bowtie2/ : symlinks to the lab's existing Bowtie2 index files. nf-core
#      finds the index by searching a DIRECTORY for *.rev.1.bt2, so the folder
#      must hold exactly one index.
#   2. <gtf>.chr.gtf : the lab's Ensembl 104 GTF with contigs renamed
#      1..22,X,Y -> chr1..chr22,chrX,chrY and MT -> chrM so they match the
#      no-alt FASTA. Other contigs (unplaced scaffolds) are dropped and counted.
#   3. Checks that GTF and blacklist contig names are all in the FASTA.
# Nothing under the source paths is modified. Re-running is safe; set FORCE=1
# to rebuild the GTF.
#
#   mkdir -p jobs
#   sbatch scripts/prepare_reference.sbatch
#
# Reads the FASTA headers (~3 GB) and the GTF (~1.5 GB) once; a few minutes.

set -u

BT2_PREFIX="${BT2_PREFIX:-/home/groups/MaxsonLab/indices/GRch38/bowtie2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index}"
FASTA="${FASTA:-/home/groups/MaxsonLab/indices/GRch38/bowtie2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta}"
SRC_GTF="${SRC_GTF:-/home/groups/MaxsonLab/indices/GRch38/STAR_index/Homo_sapiens.GRCh38.104.gtf}"
BLACKLIST="${BLACKLIST:-/home/groups/MaxsonLab/indices/GRch38/hg38.blacklist.v2.bed}"
REF_DIR="${REF_DIR:-/home/groups/MaxsonLab/indices/GRch38/nfcore_chipseq}"
FORCE="${FORCE:-0}"

fail=0
die() { echo "FATAL: $*" >&2; exit 1; }
echo "BT2_PREFIX = $BT2_PREFIX"
echo "FASTA      = $FASTA"
echo "SRC_GTF    = $SRC_GTF"
echo "BLACKLIST  = $BLACKLIST"
echo "REF_DIR    = $REF_DIR"

for f in "$FASTA" "$SRC_GTF" "$BLACKLIST"; do
    [ -r "$f" ] || { echo "ERROR: cannot read $f" >&2; fail=1; }
done
n_bt2=0
for ext in 1 2 3 4 rev.1 rev.2; do
    if [ -r "$BT2_PREFIX.$ext.bt2" ] || [ -r "$BT2_PREFIX.$ext.bt2l" ]; then
        n_bt2=$((n_bt2 + 1))
    else
        echo "ERROR: missing $BT2_PREFIX.$ext.bt2(l)" >&2
    fi
done
[ "$n_bt2" -eq 6 ] || fail=1
[ "$fail" -eq 0 ] || die "inputs missing; nothing written"

mkdir -p "$REF_DIR/bowtie2" 2>/dev/null
touch "$REF_DIR/.write_test" 2>/dev/null || die "cannot write to $REF_DIR (set REF_DIR=<dir you own>, then update params/hg38.yaml)"
rm -f "$REF_DIR/.write_test"

# --- 1. Bowtie2 index folder ------------------------------------------------
for f in "$BT2_PREFIX".*bt2 "$BT2_PREFIX".*bt2l; do
    [ -e "$f" ] || continue
    ln -sfn "$f" "$REF_DIR/bowtie2/$(basename "$f")"
done
n_rev=$(ls "$REF_DIR/bowtie2/" | grep -c '\.rev\.1\.bt2l\{0,1\}$')
echo "bowtie2/ : $(ls "$REF_DIR/bowtie2" | wc -l | tr -d ' ') files, $n_rev index"
[ "$n_rev" -eq 1 ] || die "$REF_DIR/bowtie2 must hold exactly one index, found $n_rev"

TMP="$REF_DIR/.prep_tmp.$$"
mkdir -p "$TMP" || die "cannot create $TMP"
trap 'rm -rf "$TMP"' EXIT

# --- FASTA contig names -------------------------------------------------------
if [ -r "$FASTA.fai" ]; then
    cut -f1 "$FASTA.fai" | sort -u > "$TMP/fasta.txt"
else
    grep '^>' "$FASTA" | sed -e 's/^>//' -e 's/[[:space:]].*//' | sort -u > "$TMP/fasta.txt"
fi
echo "FASTA contigs: $(wc -l < "$TMP/fasta.txt" | tr -d ' ')"

# --- 2. GTF with UCSC-style names -----------------------------------------
base="$(basename "$SRC_GTF")"
base="${base%.gz}"
OUT_GTF="$REF_DIR/${base%.gtf}.chr.gtf"
if [ -s "$OUT_GTF" ] && [ "$FORCE" != "1" ]; then
    echo "GTF exists, not rebuilt (FORCE=1 to rebuild): $OUT_GTF"
else
    case "$SRC_GTF" in *.gz) cat_cmd="gzip -dc" ;; *) cat_cmd="cat" ;; esac
    $cat_cmd "$SRC_GTF" | awk 'BEGIN { FS = OFS = "\t" }
        /^#/ { print; next }
        {
            c = $1
            if (c ~ /^chr/)                                   n = c
            else if (c ~ /^([1-9]|1[0-9]|2[0-2]|X|Y)$/)       n = "chr" c
            else if (c == "MT")                               n = "chrM"
            else { dropped++; dc[c] = 1; next }
            $1 = n; kept++; print
        }
        END {
            nd = 0; for (k in dc) nd++
            printf "GTF lines kept %d, dropped %d (on %d unplaced contigs)\n", kept, dropped, nd > "/dev/stderr"
        }' > "$TMP/out.gtf" || die "GTF conversion failed"
    mv "$TMP/out.gtf" "$OUT_GTF" || die "cannot write $OUT_GTF"
fi

# --- 3. Checks ---------------------------------------------------------------
grep -v '^#' "$OUT_GTF" | cut -f1 | sort -u > "$TMP/gtf.txt"
cut -f1 "$BLACKLIST" | sort -u > "$TMP/bl.txt"
n_gtf=$(wc -l < "$TMP/gtf.txt" | tr -d ' ')
miss_gtf=$(comm -23 "$TMP/gtf.txt" "$TMP/fasta.txt" | wc -l | tr -d ' ')
miss_bl=$(comm -23 "$TMP/bl.txt" "$TMP/fasta.txt" | wc -l | tr -d ' ')
n_gene=$(awk -F'\t' '$3 == "gene"' "$OUT_GTF" | wc -l | tr -d ' ')
echo "GTF contigs: $n_gtf, not in FASTA: $miss_gtf; genes: $n_gene"
echo "blacklist contigs not in FASTA: $miss_bl"
[ "$n_gtf" -ge 24 ] || { echo "ERROR: fewer than 24 GTF contigs after renaming" >&2; fail=1; }
[ "$miss_gtf" -eq 0 ] || { echo "ERROR: GTF contigs absent from FASTA:" >&2; comm -23 "$TMP/gtf.txt" "$TMP/fasta.txt" | head >&2; fail=1; }
[ "$miss_bl" -eq 0 ] || { echo "WARNING: blacklist contigs absent from FASTA:" >&2; comm -23 "$TMP/bl.txt" "$TMP/fasta.txt" | head >&2; }

if [ "$fail" -eq 0 ]; then
    echo "OK. params/hg38.yaml should read:"
    echo "  bowtie2_index: $REF_DIR/bowtie2"
    echo "  gtf: $OUT_GTF"
fi
exit "$fail"
