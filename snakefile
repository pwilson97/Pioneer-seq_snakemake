import pandas as pd
import os

configfile: "config.json"

localrules: all, unzip, filter, merge, trim, restrict, convert, matrix

df = pd.read_csv(config["metafile"], sep='\t', header=0, index_col=0, dtype={'Ratio': str})
samples = list(df.index)
df.index = samples
replicates = set(df["Replicate"])

def get_gz_pair(sample):
    dir = config["raw_fastq_gz_dir"]
    return tuple(os.path.join(dir + df.loc[str(sample), x]) for x in ('ForwardFastqGZ', 'ReverseFastqGZ'))

def get_ref_lib(sample):
    dir = config["ref_lib_dir"]
    return os.path.join(dir, df.loc[sample]["ReferenceLibrary"])

def get_userouts(replicate):
    dir = config["dir_names"]["align_dir"]
    samples = df.index[df["Replicate"] == replicate].tolist()
    return [f"{dir}/{sample}.userout" for sample in samples]

def get_db(replicate):
    dir = config["ref_lib_dir"]
    db = df.loc[df["Replicate"] == replicate, "ReferenceLibraryDetails"].values[0]
    return dir + db

def get_gelband_ratios(replicate):
    gel_bands = df.loc[df["Replicate"] == replicate, "GelBand"].tolist()
    ratios = df.loc[df["Replicate"] == replicate, "Ratio"].tolist()
    return [gel_band + "_" + ratio for gel_band, ratio in zip(gel_bands, ratios)]

def get_protein(replicate):
    protein = df.loc[df["Replicate"] == replicate, "Protein"].iloc[0]
    return protein

rule all:
    input:
        expand("{dir}/{sample}.userout", dir=config["dir_names"]["align_dir"],sample=samples),
        expand("{dir}/{replicate}.txt", dir=config["dir_names"]["matrix_dir"], replicate=replicates)

rule unzip:
    input:
        lambda wildcards: get_gz_pair(wildcards.sample)
    output:
        config["dir_names"]["unzip_dir"] + "/{sample}_R1.fastq",
        config["dir_names"]["unzip_dir"] + "/{sample}_R2.fastq"        		
    shell:
        "gunzip -c {input[0]} > {output[0]} && gunzip -c {input[1]} > {output[1]}"		

rule filter:
    input:
        rules.unzip.output
    output:
        config["dir_names"]["filter_dir"] + "/{sample}_R1.fq",
        config["dir_names"]["filter_dir"] + "/{sample}_R2.fq"
#    conda:
#        config["envs"]["cutadapt"]
    params:
        qscore = config["params"]["qscore"]
    log:
        config["dir_names"]["filter_log_dir"] + "/{sample}_R1.txt",
        config["dir_names"]["filter_log_dir"] + "/{sample}_R2.txt"
    shell:
       "module load gcccore/11.2.0 cutadapt/3.5 && cutadapt -q {params.qscore} -o {output[0]} {input[0]} | sed '1,6d' > {log[0]} && cutadapt -q {params.qscore} -o {output[1]} {input[1]} | sed '1,6d' > {log[1]} "

rule merge:
    input:
        rules.filter.output
    output:
        config["dir_names"]["merge_dir"] + "/{sample}.fq",
        config["dir_names"]["notmerge_dir"] + "/{sample}_R1.fq",
        config["dir_names"]["notmerge_dir"] + "/{sample}._R2.fq"
    version:
        config["tool_versions"]["vsearch"]
    params:
        maxdiffs = config["params"]["maxdiffs"],
        minovlen = config["params"]["minovlen"]
    log:
        config["dir_names"]["merge_log_dir"] + "/{sample}.txt"
    shell:
        "{version} --fastq_mergepairs {input[0]} --reverse {input[1]} --fastq_maxdiffs {params.maxdiffs} --fastq_minovlen {params.minovlen} --fastqout {output[0]} --fastaout_notmerged_fwd {output[1]} --fastaout_notmerged_rev {output[2]} &> {log}"

rule trim:
    input:
        rules.merge.output
    output:
        config["dir_names"]["trim_dir"] + "/{sample}.fq"
    conda:
        config["envs"]["cutadapt"]
    params:
        pair1 = config["params"]["pair1"],
        pair2 = config["params"]["pair2"],
        pair3 = config["params"]["pair3"],
        pair4 = config["params"]["pair4"]
    log:
        config["dir_names"]["trim_log_dir"] + "/{sample}.txt"
    shell: #not  done yet; might go from four to one param
        "module load gcccore/11.2.0 cutadapt/3.5 && cutadapt -a Pair_1={params.pair1} -a Pair_2={params.pair2} -a Pair_3={params.pair3} -a Pair_4={params.pair4} --discard-untrimmed -o {output} {input[0]} | sed '1,4d' > {log}"

rule restrict:
    input:
        rules.trim.output
    output:
        config["dir_names"]["restrict_dir"] + "/{sample}.fq"
    conda:
        config["envs"]["cutadapt"]
    params:
        min = config["params"]["minimum_length"],
        max = config["params"]["maximum_length"]
    log:
        config["dir_names"]["restrict_log_dir"] + "/{sample}.txt"
    shell:
        "module load gcccore/11.2.0 cutadapt/3.5 && cutadapt -m {params.min} -M {params.max} -o {output} {input} | sed '1,4d' > {log}"

rule convert:
    input:
        rules.restrict.output
    output:
        config["dir_names"]["convert_dir"] + "/{sample}.fa"
    version:
        config["tool_versions"]["fastx_toolkit"]
    shell:
        "{version} -n -i {input} -o {output}"

rule align:
    input:
        rules.convert.output,
#        "reference_libraries/Nucs_v2_DB_no_primers.fa"
        lambda wildcards: get_ref_lib(wildcards.sample)
    output:
        config["dir_names"]["align_dir"] + "/{sample}.userout",
        config["dir_names"]["align_dbout_dir"] + "/{sample}.dbout"
    version:
        config["tool_versions"]["vsearch"]
    params:
        id = config["params"]["global_clustering_threshold"],
        threads = config["params"]["threads"],
        maxrejects = config["params"]["maxrejects"],
        maxaccepts = config["params"]["maxaccepts"],
        qmask = config["params"]["qmask"],
        dbmask = config["params"]["dbmask"],
        userfields = config["params"]["userfields"],
        mincols = config["params"]["mincols"]
    log:
        config["dir_names"]["align_log_dir"] + "/{sample}.txt"
    shell:
        "{version} --usearch_global {input[0]} --db {input[1]} --dbmatched {output[1]} --sizeout --id {params.id} --threads {params.threads} --maxrejects {params.maxrejects} --maxaccepts {params.maxaccepts} --qmask {params.qmask} --dbmask {params.dbmask} --userfields {params.userfields} --userout {output[0]} --top_hits_only --mincols {params.mincols} &> {log}"

rule matrix:
    input:
        userouts = lambda wildcards: get_userouts(wildcards.replicate),
        metafile = config["metafile"],
        db = lambda wildcards: get_db(wildcards.replicate)
    output:
        config["dir_names"]["matrix_dir"] + "/{replicate}.txt"
    params:
        gelband_ratios = lambda wildcards: get_gelband_ratios(wildcards.replicate),
        protein = lambda wildcards: get_protein(wildcards.replicate)
    script:
        config["matrix_script"]
