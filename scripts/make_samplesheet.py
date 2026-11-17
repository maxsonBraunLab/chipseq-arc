#!/usr/bin/env python3
"""Draft an nf-core/chipseq 2.1.0 samplesheet from FASTQ file names.

Expects the maxsonBraunLab/chip_seq naming convention:
    {condition}_{replicate}_{mark}_R1.fastq.gz   (+ _R2 for paired-end)
    {condition}_{replicate}_{mark}.fastq.gz      (single-end)
Illumina suffixes are tolerated: ..._{mark}_S12_L001_R1_001.fastq.gz.
{condition} may contain underscores; {replicate} must be an integer;
{mark} must not contain underscores.

Each mark becomes sample "{condition}_{mark}" with antibody "{mark}".
Its control is "{condition}_{CONTROL}" with the same replicate number if that
exists, else the only replicate of that control, else blank (the run then
proceeds without a control for it, and the script warns).
Several lanes of one library become several rows with the same sample and
replicate, which nf-core merges.

The output is a DRAFT: read it before running the pipeline.

    python3 scripts/make_samplesheet.py /path/to/fastqs --control Input -o samplesheet.csv
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("fastq_dir")
ap.add_argument("--control", default="Input", help="mark name used for controls (default Input; e.g. IgG)")
ap.add_argument("-o", "--output", default="samplesheet.csv")
ap.add_argument("--force", action="store_true", help="overwrite an existing output file")
args = ap.parse_args()

if os.path.exists(args.output) and not args.force:
    sys.exit("{} exists; use --force to overwrite".format(args.output))

EXT = re.compile(r"\.(fastq|fq)\.gz$")
ILLUMINA = re.compile(r"_S\d+(_L\d{3})?$|_L\d{3}$")
READ = re.compile(r"_(R?[12])(_001)?$")

libs = defaultdict(lambda: {"R1": [], "R2": [], "SE": []})   # (cond, rep, mark, lane_key) -> reads
problems = []
for root, _, files in os.walk(args.fastq_dir, followlinks=True):
    for f in sorted(files):
        if not EXT.search(f):
            continue
        path = os.path.abspath(os.path.join(root, f))
        stem = EXT.sub("", f)
        m = READ.search(stem)
        read = "SE"
        if m:
            read = "R" + m.group(1)[-1]
            stem = stem[:m.start()]
        lane = ""
        m = ILLUMINA.search(stem)
        if m:
            lane = m.group(0)
            stem = stem[:m.start()]
        parts = stem.split("_")
        if len(parts) < 3 or not parts[-2].isdigit():
            problems.append("cannot parse {} as condition_replicate_mark".format(f))
            continue
        cond, rep, mark = "_".join(parts[:-2]), int(parts[-2]), parts[-1]
        libs[(cond, rep, mark, lane)][read].append(path)

rows = []
reps_of = defaultdict(set)
for (cond, rep, mark, lane), r in libs.items():
    reps_of[(cond, mark)].add(rep)
for (cond, rep, mark, lane), r in sorted(libs.items()):
    if r["SE"] and (r["R1"] or r["R2"]):
        problems.append("{}_{}_{}{}: mixes single-end and paired files".format(cond, rep, mark, lane))
        continue
    if r["R1"] or r["R2"]:
        if len(r["R1"]) != 1 or len(r["R2"]) != 1:
            problems.append("{}_{}_{}{}: expected one R1 and one R2, found {} and {}".format(
                cond, rep, mark, lane, len(r["R1"]), len(r["R2"])))
            continue
        fq1, fq2 = r["R1"][0], r["R2"][0]
    else:
        if len(r["SE"]) != 1:
            problems.append("{}_{}_{}{}: {} single-end files".format(cond, rep, mark, lane, len(r["SE"])))
            continue
        fq1, fq2 = r["SE"][0], ""
    if mark == args.control:
        rows.append([cond + "_" + mark, fq1, fq2, rep, "", "", ""])
        continue
    ctrl_reps = reps_of.get((cond, args.control), set())
    if rep in ctrl_reps:
        crep = rep
    elif len(ctrl_reps) == 1:
        crep = next(iter(ctrl_reps))
        sys.stderr.write("NOTE: {}_{} replicate {} uses {}_{} replicate {} (the only one)\n".format(
            cond, mark, rep, cond, args.control, crep))
    else:
        crep = None
    if crep is None:
        sys.stderr.write("WARNING: no unambiguous {} control for {}_{} replicate {}; control left blank\n".format(
            args.control, cond, mark, rep))
        rows.append([cond + "_" + mark, fq1, fq2, rep, mark, "", ""])
    else:
        rows.append([cond + "_" + mark, fq1, fq2, rep, mark, cond + "_" + args.control, crep])

if problems:
    sys.stderr.write("make_samplesheet: {} problem(s); nothing written\n".format(len(problems)))
    for p in problems:
        sys.stderr.write("  " + p + "\n")
    sys.exit(1)
if not rows:
    sys.exit("no FASTQ files found under " + args.fastq_dir)

with open(args.output, "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["sample", "fastq_1", "fastq_2", "replicate", "antibody", "control", "control_replicate"])
    w.writerows(rows)
n_ctrl = sum(1 for r in rows if r[4] == "")
print("wrote {}: {} rows, {} samples ({} control rows)".format(
    args.output, len(rows), len({r[0] for r in rows}), n_ctrl))
