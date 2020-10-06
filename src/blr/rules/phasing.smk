"""
Rules related to phasing of variants (called or reference set)
"""
import os
import pandas as pd
import numpy as np
from collections import Counter, defaultdict

rule hapcut2_extracthairs:
    """Extract heterozygous variants covered by alignments in BAM"""
    output:
        unlinked = temp("{base}.calling.unlinked.txt")
    input:
        bam = "{base}.calling.bam",
        vcf = "{base}.phaseinput.vcf",
    log: "{base}.hapcut2_extracthairs.log"
    params:
        indels = "1" if config["phase_indels"] else "0"
    shell:
        "extractHAIRS"
        " --10X 1"
        " --indels {params.indels}"
        " --realign_variants 1"  # Improves overall error-rate
        " --ref {config[genome_reference]}"
        " --bam {input.bam}"
        " --VCF {input.vcf}"
        " --out {output.unlinked} 2> {log}"


rule hapcut2_linkfragments:
    """Link heterozygous variants together using barcode information"""
    output:
        linked = temp("{base}.calling.linked.txt")
    input:
        bam = "{base}.calling.bam",
        bai = "{base}.calling.bam.bai",
        vcf = "{base}.phaseinput.vcf",
        unlinked = "{base}.calling.unlinked.txt"
    log: "{base}.hapcut2_linkfragments.log"
    shell:
        "LinkFragments.py"
        " --bam {input.bam}"
        " -v {input.vcf}"
        " --fragments {input.unlinked}"
        " --out {output.linked}"
        " --distance {config[window_size]} &> {log}"


rule hapcut2_phasing:
    """Phase heterozygous variants using HapCUT2. Output phased VCF"""
    output:
        phase = "{base}.calling.phase",
        phased_vcf = "{base}.calling.phased.vcf"
    input:
        linked = "{base}.calling.linked.txt",
        vcf = "{base}.phaseinput.vcf",
    log: "{base}.hapcut2_phasing.log"
    shell:
        "hapcut2"
        " --nf 1"
        " --fragments {input.linked}"
        " --vcf {input.vcf}"
        " --out {output.phase}"
        " --error_analysis_mode 1"
        " --outvcf 1 &> {log}"
        " && "
        "mv {output.phase}.phased.VCF {output.phased_vcf}"


rule hapcut2_stats:
    """Get phasing statistics relative the ground truth. See https://github.com/vibansal/HapCUT2/tree/master/utilities
    for details. """
    output:
        stats = "final.phasing_stats.txt"
    input:
        vcf1 = "final.phased.vcf",
    params:
        vcf2 = f" -v2 {config['phasing_ground_truth']}" if config['phasing_ground_truth'] else "",
        indels = " --indels" if config["phase_indels"] else ""
    shell:
        "blr calculate_haplotype_statistics"
        " -v1 {input.vcf1}"
        " {params.vcf2}"
        " {params.indels}"
        " -o {output.stats}"


def get_haplotag_input(wildcards):
    inputfiles = {"bam": f"{wildcards.base}.calling.bam"}
    if config["haplotag_tool"] == "blr":
        inputfiles.update({
            "hapcut2_phase_file": f"{wildcards.base}.calling.phase"
        })
    elif config["haplotag_tool"] == "whatshap":
        inputfiles.update({
            "vcf": f"{wildcards.base}.calling.phased.vcf.gz",
            "vcf_index": f"{wildcards.base}.calling.phased.vcf.gz.tbi",
        })
    return inputfiles


rule haplotag:
    """
    Transfer haplotype information from the phased VCF file to the bam file.
    Adds HP tag with haplotype (1 or 2) and PS tag with phase set information.
    """
    output:
        bam = "{base}.calling.phased.bam"
    input:
        unpack(get_haplotag_input)
    log: "{base}.haplotag.log"
    run:
        ignore_readgroups = "--ignore-read-groups" if config["reference_variants"] else ""

        commands = {
            "whatshap":
                "whatshap haplotag"
                " {input.vcf}"
                " {input.bam}"
                " --linked-read-distance-cutoff {config[window_size]}"
                " --reference {config[genome_reference]}"
                " -o {output.bam}"
                " {ignore_readgroups}",
            "blr":
                "blr phasebam"
                " --molecule-tag {config[molecule_tag]}"
                " --phase-set-tag {config[phase_set_tag]}"
                " --haplotype-tag {config[haplotype_tag]}"
                " --min-mapq {config[min_mapq]}"
                " -o {output.bam}"
                " {input.bam}"
                " {input.hapcut2_phase_file}"
        }
        command = commands[config["haplotag_tool"]]
        shell(command + " 2> {log}")



rule build_config:
    """
    Builds a config file required for running NAIBR.
    """
    output:
        config = "{base}.naibr.config"
    input:
        bam = "{base}.calling.phased.bam",
        index = "{base}.calling.phased.bam.bai"
    log: "{base}.build_config.log"
    params:
        cwd = os.getcwd()
    shell:
        "blr naibrconfig"
        " --bam-file {input.bam}"
        " --outdir {params.cwd}"
        " --distance {config[window_size]}"
        " --min-mapq {config[naibr_min_mapq]}"
        " --min-sv 1000"
        " --threads 1"
        " --min-overlaps 3"
        " --output {output.config}"
        " 2> {log}"


rule get_naibr_path:
    """Symlink existing NAIBR repo if specified or clone the latest version"""
    output:
        naibr_path = temp(directory("NAIBR"))
    run:
        if config["naibr_path"]:
            shell("ln -s {config[naibr_path]} {output.naibr_path}")
        else:
            shell("git clone https://github.com/raphael-group/NAIBR.git {output.naibr_path}")


rule lsv_calling:
    """
    Runs NAIBR for LSV calling. This involves activating a python2 env, changing wd, running and changing back wd and
    env.
    """
    output:
        results = "{base}.naibr_sv_calls.tsv"
    input:
        config = "{base}.naibr.config",
        naibr_path = "NAIBR"
    log: "{base}.lsv_calling.log"
    threads: 20
    conda: "../naibr-environment.yml"
    params:
        cwd = os.getcwd()
    shell:
        "cd {input.naibr_path}"
        " &&"
        " python"
        " NAIBR.py"
        " {params.cwd}/{input.config}"
        " > {params.cwd}/{log}"
        " &&"
        " cd -"
        " &&"
        " mv NAIBR_SVs.bedpe {output.results}"


rule format_naibr_bedpe:
    output:
        bedpe = "final.naibr_sv_calls.bedpe"
    input:
        tsv = "final.naibr_sv_calls.tsv"
    run:
        # Translation based on this: https://github.com/raphael-group/NAIBR/issues/11
        # Note that this is not entirely accurate as the more complex variants are possible.
        sv_types = {"+-": "DEL", "++": "INV", "--": "INV", "-+": "DUP"}
        zygosity = {"1,1":"HOM" ,"2,2":"HOM" , "1,2":"HET" ,"2,1":"HET"}  # From https://github.com/raphael-group/NAIBR/issues/10
        names = ["Chr1",	"Break1",	"Chr2",	"Break2", "SplitMolecules", "DiscordantReads", "Orientation",
                 "Haplotype", "Score", "PassFilter"]
        data = pd.read_csv(input.tsv, sep="\t", header=0, names=names)
        with open(output.bedpe, "w") as file:
            for nr, row in enumerate(data.itertuples(index=False)):
                # Use BEDPE format according to 10x: https://support.10xgenomics.com/genome-exome/software/pipelines/latest/output/bedpe
                # Supported by IGV (see https://github.com/igvteam/igv/wiki/BedPE-Support)
                print(
                    f"chr{row.Chr1}",
                    row.Break1,
                    row.Break1,
                    f"chr{row.Chr2}",
                    row.Break2,
                    row.Break2,
                    f"call_{nr}",
                    row.Score,
                    "+",
                    "+",
                    "." if row.PassFilter == "PASS" else "Filtered",
                    "Type={};Zygosity={};Split_molecules={};Discordant_reads={};Orientation={};Haplotype={}".format(
                        sv_types[row.Orientation],
                        zygosity.get(row.Haplotype, "UNK"),
                        row.SplitMolecules,
                        row.DiscordantReads,
                        row.Orientation,
                        row.Haplotype,
                    ),
                    sep="\t",
                    file=file
                )

rule aggregate_sv_sizes:
    output:
        tsv = "final.sv_sizes.tsv"
    input:
        tsv = "final.naibr_sv_calls.tsv"
    run:
        sv_trans = {"+-": "DEL", "++": "INV", "--": "INV", "-+": "DUP"}
        sv_types = ["DEL", "INV", "DUP"]
        bins = [0, 1000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000, 1_000_000_000]
        bin_labels = ["0-1k", "1k-10k", "10k-100k", "100k-1M", "1M-10M", "10M-100M", "100M-1G"]
        names = ["Chr1", "Break1", "Chr2", "Break2", "SplitMolecules", "DiscordantReads", "Orientation", "Haplotype",
                 "Score", "PassFilter"]
        data = pd.read_csv(input.tsv, sep="\t", header=0, names=names)
        counts = Counter()
        lengths = defaultdict(list)
        with open(output.tsv, "w") as file:
            # Remove filtered SVs
            data = data[data["PassFilter"] == "PASS"]

            for nr, row in enumerate(data.itertuples(index=False)):
                sv_type = sv_trans[row.Orientation]
                if row.Chr1 != row.Chr2:
                    counts[("Interchrom",sv_type)] += 1
                else:
                    lengths[sv_type].append(abs(row.Break1 - row.Break2))


            for sv_type in sv_types:
                bin_counts, _ = np.histogram(lengths[sv_type], bins=bins)

                for count, label in zip(bin_counts, bin_labels):
                    counts[(label, sv_type)] += count

            print("Size", *sv_types, sep="\t", file=file)
            for label in ["Interchrom"] + bin_labels:
                print(label, *[counts[(label, sv_type)] for sv_type in sv_types], sep="\t", file=file)
