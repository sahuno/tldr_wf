## Prepare TLDR regions for IGV-Reports

```python
import pandas as pd


from IPython.display import IFrame

# # 1. Read your TLDR breakpoint report
# tldr = pd.read_csv("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/tldr_wf/outputs/postTLDR/data/Insertions_pass.tsv", sep="\t")
# tldr["SampleReads"]
# tldr["EmptyReads"]
# NumSamples	SampleReads	EmptyReads
# # 2. Build the regions DataFrame
# regions = pd.DataFrame({
#     "seq":        tldr["Chrom"],
#     "start":      tldr["Start"],
#     "end":        tldr["End"],
#     "name":       tldr.apply(
#                      lambda r: f"{r['new_samples_name']}_{r['Chrom']}_{r['Start']}", 
#                      axis=1
#                  ),
#     "Family":     tldr["Family"],
#     "Subfamily":  tldr["Subfamily"],
#     "sample":     tldr["new_samples_name"]
# })

# # 3. Save in the IGV-Reports recommended format
# regions.to_csv("tldr_regions.tsv", sep="\t", index=False)

# print("Wrote", len(regions), "regions to tldr_regions.tsv")


# refmm10="/data1/greenbab/database/mm10/mm10.fa"
# # read the break point files
# tldr = pd.read_csv("/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/tldr_wf/outputs/postTLDR/data/Insertions_pass.tsv", sep="\t")

# python -m igv_reports.report \
#     tldr_regions.tsv \
#     --genome mm10 \
#     --fasta $refmm10 \
#     --sequence 1 --begin 2 --end 3 --idlink 4 \
#     --info-columns Family Subfamily sample \
#     --flanking 200 \
#     --tracks /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/mark_duplicates/D-A-1_4000/D-A-1_4000_modBaseCalls_sorted_dup.bam /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/mark_duplicates/D-0-1_4000/D-0-1_4000_modBaseCalls_sorted_dup.bam \
#     --output tldr_igv_report.html



# python -m igv_reports.report \
#     tldr_regions.tsv \
#     --genome mm10 \
#     --fasta /data1/greenbab/database/mm10/mm10.fa \
#     --sequence {Chrom} --begin {Start} --end {End} --idlink {keyword} \
#     --info-columns Family Subfamily sample \
#     --flanking 200 \
#     --tracks ${tumorBamPaths} /data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/mark_duplicates/D-0-1_4000/D-0-1_4000_modBaseCalls_sorted_dup.bam \
#     --output tldr_{keyword}igv_report.html




import pandas as pd
import yaml
import os
import ace_tools as tools

# Placeholder: replace with your actual TLDR file path and samples.yaml path
tldr_file = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/tldr_wf/outputs/postTLDR/data/Insertions_pass.tsv"
samples_yaml = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/tldr_wf/configs/samples_bams_all.yaml"
NormalBam="/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/mark_duplicates/D-0-1_5000_4000/D-0-1_5000_4000_modBaseCalls_sorted_dup.bam"
# 1. Load TLDR table
tldr = pd.read_csv(tldr_file, sep="\t")

# 2. Extract base BAM filenames and add .bam suffix
tldr['tumorBamFileName'] = tldr['SampleReads'].str.split('|').str[0] + '.bam'
tldr['tumorBamSampleName'] = tldr['SampleReads'].str.split('|').str[0]

# 3. Load samples.yaml
with open(samples_yaml, 'r') as fh:
    cfg = yaml.safe_load(fh)

# cfg['samples'] is a dict from sample key to full BAM path
# Build a mapping from bam filename -> full path
bam_mapping = { os.path.basename(path): path
                for path in cfg['samples'].values() }

# 4. Map full paths to a new column
tldr['tumorBamPaths'] = tldr['tumorBamFileName'].map(bam_mapping)


tldr['Chrom_Start_End_Strand_Family_Subfamily'] = tldr.apply(
    lambda r: f"{r.Chrom}_{r.Start}_{r.End}_{r.Strand}_{r.Family}_{r.Subfamily}",
    axis=1
)

strand_map = {'+':'plus', '-':'minus'}
tldr['StrandWord'] = tldr['Strand'].map(strand_map)

# 2) Build the composite
tldr['keyword'] = (
    tldr['Chrom'].astype(str) + '_' +
    tldr['Start'].astype(str) + '_' +
    tldr['End'].astype(str) + '_' +
    tldr['StrandWord'] + '_' +
    tldr['Family'] + '_' +
    tldr['Subfamily']
)



# import subprocess
# import shlex

# GENOME   = "mm10"
# FASTA    = "/data1/greenbab/database/mm10/mm10.fa"
# FLANKING = "200"

# for _, row in tldr[:3].iterrows():
#     chrom   = row.Chrom
#     start   = int(row.Start)
#     end     = int(row.End)
#     name    = row.keyword
#     bam     = row.tumorBamPaths
#     sampleName = row.tumorBamSampleName
#     outname = f"tldr_{name}_igv_report.html"
#     cmd = [
#         "python", "-m", "igv_reports.report",
#         "--genome", GENOME,
#         "--fasta", FASTA,
#         "--sequence", chrom, "--begin", str(start), "--end", str(end), "--idlink", name,
#         "--info-columns", "Family", "Subfamily", "new_samples_name",
#         "--flanking", FLANKING,
#         "--tracks", bam,
#         "--output", outname
#     ]
#     print("Running:", " ".join(shlex.quote(x) for x in cmd))
#     subprocess.run(cmd, check=True)

    
# # 5. Preview the new columns
# tools.display_dataframe_to_user("TLDR with BAM Paths", 
#                                 tldr[['SampleReads', 'tumorBamFileName', 'tumorBamPaths']].head())
import tempfile, subprocess, shlex

GENOME   = "mm10"
FASTA    = "/data1/greenbab/database/mm10/mm10.fa"
FLANKING = "200"

# NormalBam="/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/mark_duplicates/D-0-1_5000_4000/D-0-1_5000_4000_modBaseCalls_sorted_dup.bam"
NormalBam = ["/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/mark_duplicates/D-0-1_5000_4000/D-0-1_5000_4000_modBaseCalls_sorted_dup.bam"]

for _, row in tldr.iterrows():
    name    = row.keyword
    tracks     = [row.tumorBamPaths] + NormalBam
    sampleName = row.tumorBamSampleName
    outname = f"tldr_{name}_{sampleName}_igv_report.html"
    # 1) write one‐row TSV
    with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as tf:
        tf.write("sequence\tbegin\tend\tidlink\n")
        tf.write(f"{row.Chrom}\t{int(row.Start)}\t{int(row.End)}\t{name}\n")
        sites = tf.name
    # 2) call igv_reports with column‐indices
    cmd = [
        "python", "-m", "igv_reports.report",
         sites,
        "--genome", GENOME,
        "--fasta",  FASTA,
        "--sequence", "1",
        "--begin",    "2",
        "--end",      "3",
        "--idlink",   "4",
        "--info-columns", "Family", "Subfamily", "new_samples_name",
        "--flanking", FLANKING,
        "--tracks",   *tracks,
        "--output",   outname
    ]
    print("Running:", " ".join(shlex.quote(x) for x in cmd))
    subprocess.run(cmd, check=True)
