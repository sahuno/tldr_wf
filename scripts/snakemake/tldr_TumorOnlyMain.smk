import os
import glob
from datetime import datetime

#purpose: phase snps and 5mc

# rm -r .snakemake logs results benchmarks

# Assign the current date to a variable
current_datetime = datetime.now()
formatted_datetime = current_datetime.strftime('%Y%m%d_%H_%M_%S')

parent_dir = "/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/tldr_wf/"
#set species
set_species = "mouse"
configfile: parent_dir + "configs/config.yaml"
configfile: parent_dir + "configs/samples_bams_all.yaml" #mouse samples

# configfile: parent_dir + "configs/samples_bams_5000sampleRate.yaml" #mouse samples; test run


# teref.mouse.fa
print(config["samples"]) #sanity check


rule all:
    input:
        expand('results/tldr/{samples}/done.{samples}.txt', samples=config["samples"]),
        expand('results/tldr/{samples}/plot_scripts.{samples}.txt', samples=config["samples"]),
        expand('results/runMethylArtist/{samples}/done.{samples}.txt', samples=config["samples"])

rule tldr:
    input:
        lambda wildcards: config["samples"][wildcards.samples]
    params:
        # RefNormalBam="/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/results/mark_duplicates/D-0-1_5000/D-0-1_5000_modBaseCalls_sorted_dup.bam",
        RefNormalBam="/data1/greenbab/projects/triplicates_epigenetics_diyva/DNA/preprocessed/mergedbams_modkit/results/mark_duplicates/D-0-1_5000_4000/D-0-1_5000_4000_modBaseCalls_sorted_dup.bam",
        reference_genome=lambda wildcards: config["mm10"] if set_species == "mouse" else config["hg38"],
        te_re=lambda wildcards: config["consRef_mm10"] if set_species == "mouse" else config["consRef_human"],
        software_dir=config["software_dir"],
        out_dir='results/tldr/{samples}/',
        threads=12,
        chromFile="/data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/tldr_wf/configs/chromFile.txt"
    singularity: "/data1/greenbab/projects/janjigiy_su2c_sasha/te.sif"
    output:
        done='results/tldr/{samples}/done.{samples}.txt',
        plot_scripts='results/tldr/{samples}/plot_scripts.{samples}.txt'
    log:
      "logs/tldr/{samples}/{samples}.log"
    shell:
        """
        tldr --bams {input},{params.RefNormalBam} \
        --outbase "results/tldr/{wildcards.samples}" --procs 10 --elts /tldr/ref/teref.human.fa \
        --ref {params.reference_genome} \
        --chroms {params.chromFile} \
        --min_te_len 100 --max_cluster_size 100 \
        --trdcol --detail_output --methylartist --keep_pickles 2> {log} && \
        cat results/tldr/{wildcards.samples}/*.methylartist.cmd > {output.plot_scripts} \
        touch {output.done}
        """
#   --chroms {params.chromFile} \
        # f"results/tldr/{samples}/{wildcards.samples}.*.methylartist.cmd" >> {output.plot_scripts}

rule runMethylArtist:
    input:
        plot_scripts="results/tldr/{samples}/plot_scripts.{samples}.txt"
    params:
        threads=12,
        tldr_dir=lambda wildcards, input: "results/tldr/{}".format(wildcards.samples),
        outdir = lambda wildcards: f"results/runMethylArtist/{wildcards.samples}"
    output:
        done='results/runMethylArtist/{samples}/done.{samples}.txt',
        dir=directory('results/runMethylArtist/{samples}/')
    log:
      "logs/runMethylArtist/{samples}/{samples}.log"
    shell:
        r"""
        # 1) activate your methylartist env
        source /home/ahunos/miniforge3/etc/profile.d/conda.sh && conda activate methylartist

        # 2) make & cd into the sample-specific output dir
        mkdir -p {params.outdir}
        pushd {params.outdir}

        # 3) for each line in the plot_scripts file, run it and append logs
        while IFS='' read -r cmd || [ -n "$cmd" ]; do
            echo ">>> $cmd"                             >> ../../{log}
            $cmd                                      >> ../../{log} 2>&1
        done < ../../{input.plot_scripts}

        # 4) come back and “touch” the done file
        popd
        touch {output.done}
        """


# rule runMethylArtist:
#     input:
#         tldr_done="results/tldr/{samples}/done.{samples}.txt"
#     params:
#         threads=12,
#         tldr_dir=lambda wildcards, input: "results/tldr/{}".format(wildcards.samples)
#     output:
#         done='results/runMethylArtist/{samples}/done.{samples}.txt',
#         dir=directory('results/runMethylArtist/{samples}/')
#     log:
#       "logs/runMethylArtist/{samples}/{samples}.log"
#     shell:
#         """
#         # Activate the methylartist environment.
#         source /home/ahunos/miniforge3/etc/profile.d/conda.sh && conda activate methylartist
#         mkdir -p results/runMethylArtist/{wildcards.samples}
#         # Use the directory from params.tldr_dir to search for all command files.
#         for cmd in $(find {params.tldr_dir} -maxdepth 1 -name '*.methylartist.cmd'); do
#             echo "Running $cmd" >> {log}
#             bash "$cmd" 
#             # Optionally move any generated PNG or TSV files into the run directory.
#             mv *.png results/runMethylArtist/{wildcards.samples}/ || true
#             mv *.tsv results/runMethylArtist/{wildcards.samples}/ || true
#         done 2>> {log}
#         touch {output.done}
#         """


#run on slurm
# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/tldr_wf/scripts/snakemake/tldr_TumorOnlyMain.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurmMinimal --jobs 10 --cores all --use-singularity --singularity-args "'-B "/data1/greenbab"'" --keep-going --forceall -np

# methylartist             /home/ahunos/miniforge3/envs/methylartist
# bash /data1/greenbab/projects/methyl_benchmark_spectrum/data/preprocessed/L1_insertions_output/044T_v14_modBaseCalls_sorted_dup_044N_v14_modBaseCalls_sorted_dup/044T_v14_modBaseCalls_sorted_dup.83606efe-63c6-424d-853d-1d9c8c356d45.methylartist.cmd
# 044N_v14_modBaseCalls_sorted_dup.83606efe-63c6-424d-853d-1d9c8c356d45.te.83606efe_0_1592.mh.cohort.ms1.smw18.locus.meth.png


# snakemake -s /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/tldr_wf/scripts/snakemake/tldr_TumorOnlyMain.smk --workflow-profile /data1/greenbab/users/ahunos/apps/workflows/methylation_workflows/DNAme_Ref_LINE1/config/slurmMinimal --jobs 10 --cores all --use-singularity --singularity-args "'-B "/data1/greenbab"'" --keep-going --forceall -np
