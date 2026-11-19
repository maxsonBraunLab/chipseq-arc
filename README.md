# nf-core/chipseq on ARC (MaxsonBraun lab)

Lab launcher for [nf-core/chipseq](https://nf-co.re/chipseq) **2.1.0** on OHSU ARC:
Bowtie2 on the lab's GRCh38 index, Singularity containers, SLURM `batch`.
The pipeline itself is not copied here; Nextflow fetches the pinned release.
Paste blocks are one command per line.

| file | purpose |
|---|---|
| `conf/arc.config` | SLURM `batch`, no `--account`, Singularity, resource caps, lab Bowtie2 flags |
| `params/hg38.yaml` | GRCh38 FASTA / Bowtie2 index / GTF / blacklist, MACS3 settings |
| `run_chipseq_arc.sbatch` | driver job: preflight checks, then `nextflow run nf-core/chipseq -r 2.1.0` |
| `scripts/prepare_reference.sbatch` | one-time: Bowtie2 index folder + chr-named GTF |
| `scripts/make_samplesheet.py` | drafts `samplesheet.csv` from lab FASTQ names |
| `samplesheet_template.csv` | the 7-column nf-core format |

## One-time setup

1. Can compute nodes reach the internet? The driver downloads containers on
   first use. Expect a line starting `HTTP/`:

```
srun -p batch -t 00:05:00 curl -sI --max-time 15 https://depot.galaxyproject.org/singularity/
```

   If nothing prints, follow *Offline compute nodes* below instead of step 2.

2. Fetch the pipeline and its plugin once, on the login node:

```
module load nextflow/25.10.0
export NXF_VER=24.10.0
nextflow pull nf-core/chipseq -r 2.1.0
nextflow plugin install nf-validation@1.1.3
```

3. Prepare the reference (a few minutes):

```
mkdir -p jobs
sbatch scripts/prepare_reference.sbatch
```

   It writes `/home/groups/MaxsonLab/indices/GRch38/nfcore_chipseq/` and ends
   with `OK.` plus the two paths `params/hg38.yaml` expects. Nothing in the
   existing index folders is changed. It refuses to start if it cannot write
   there; then run it with `REF_DIR=<a folder you own>` and edit the two
   paths in `params/hg38.yaml`.

4. Container images are cached in
   `/home/groups/MaxsonLab/software/singularity-containers/nf-core`
   (override with `NXF_SINGULARITY_CACHEDIR`). The first run fills it
   (a few GB), so it is slower than later runs.

## Per project

```
cd /home/groups/MaxsonLab/braun/analysis
git clone https://github.com/maxsonBraunLab/chipseq-arc.git my_chip_project
cd my_chip_project
python3 scripts/make_samplesheet.py /path/to/fastqs --control IgG -o samplesheet.csv
```

Read `samplesheet.csv` before running: it is a draft. FASTQs must follow
`{condition}_{replicate}_{mark}_R1.fastq.gz` (Illumina `_S#_L00#_R1_001`
suffixes are fine; lanes become extra rows of the same library). A mark whose
replicate has no matching control falls back to the only control replicate,
with a NOTE; if there is no unambiguous control the field is left blank with a
WARNING. One run uses one peak mode, so put broad marks (H3K27me3, H3K9me3)
in a separate samplesheet and run with `narrow_peak: false`.

```
mkdir -p jobs
PREVIEW=1 bash run_chipseq_arc.sbatch
sbatch run_chipseq_arc.sbatch
```

`PREVIEW=1` runs every preflight check plus `nextflow -preview` (validates
parameters, samplesheet and DAG) on the login node, without submitting work.
The driver job (2 CPU, 10 GB, 120 h with `--qos=long_jobs`) submits one SLURM
job per process; `work/` goes to
`/home/groups/MaxsonBraunScratch/$USER/nxf-work/<folder name>`.

Monitoring and restarting:

```
squeue -u $USER
sacct -u $USER -S today --format=JobID%-12,JobName%-40,State,Elapsed,MaxRSS
tail -n 40 jobs/chipseq_nf_*.log
```

To continue after a failure or cancellation, submit again; `-resume` is built in.

Results: `results/multiqc/` (QC report), `results/bowtie2/merged_library/`
(filtered BAMs, bigWigs, `macs3/` peaks and consensus peak sets),
`results/pipeline_info/`. Delete the work directory once results are checked.

## Lab choices in this configuration (review)

- **Bowtie2 flags** match `maxsonBraunLab/chip_seq`: paired-end
  `--local --very-sensitive-local --no-unal --no-mixed -I 10 -X 700`,
  single-end `--local`. nf-core's own default is end-to-end with no fragment
  cap; `-X 700` drops pairs with fragments longer than 700 bp.
- **Peaks**: `narrow_peak: true` (nf-core default is broad). MACS3 uses `BAMPE`
  for paired-end data automatically. `macs_gsize` 2.7e9 as in `chip_seq`.
- **Annotation**: the lab's Ensembl 104 GTF with chromosomes renamed to the
  FASTA's `chr` names; unplaced-scaffold genes are dropped (count printed).
- **Consensus**: `min_reps_consensus: 1`.
- **Limits**: tasks capped at 16 CPU / 128 GB / 36 h (the `batch` limit
  without `long_jobs`); nodes `cnode-07-[5-29]` and `cnode-08-[1-35]` are
  excluded, copied from the lab's working ScaleTagToolkit driver (reason not
  recorded; remove the `clusterOptions` line in `conf/arc.config` to allow them).

## Offline compute nodes

If step 1 printed nothing, download the pipeline and all images on the login
node with nf-core tools, then point the driver at the local copy:

```
conda create -y -n nf-core -c conda-forge -c bioconda nf-core
conda activate nf-core
export NXF_SINGULARITY_CACHEDIR=/home/groups/MaxsonLab/software/singularity-containers/nf-core
nf-core pipelines download nf-core/chipseq --revision 2.1.0 --container-system singularity --container-cache-utilisation amend --compress none --outdir /home/groups/MaxsonLab/software/nf-core/chipseq
PIPELINE=/home/groups/MaxsonLab/software/nf-core/chipseq/2_1_0 sbatch --export=ALL run_chipseq_arc.sbatch
```

The driver checks that a local `PIPELINE` folder really is version 2.1.0.

## Verified so far (off-cluster, 2026-11-17)

With Nextflow 24.10.0 against the 2.1.0 source: the merged configuration
(SLURM `batch`, caps, Singularity, the Bowtie2 flags replacing the pipeline's
own); `-preview` passing parameter and samplesheet validation and building the
Bowtie2 DAG without an index-build step; the driver's preflight catching a
missing FASTQ, an unknown control and a bad reference path;
`prepare_reference.sbatch` and `make_samplesheet.py` on synthetic files.
**Not yet verified on ARC**: internet access from compute nodes, write access
to the index and container folders, the real GTF/FASTA contig names (the prep
script checks them), and any run on real reads. Start with one small
paired-end sample.
