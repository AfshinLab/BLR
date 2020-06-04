"""
Rules related to phasing of variants (called or reference set)
"""


rule hapcut2_extracthairs:
    """Extract heterozygous variants covered by alignments in BAM"""
    output:
        unlinked = "{base}.calling.unlinked.txt"
    input:
        bam = "{base}.calling.bam",
        vcf = "{base}.phaseinput.vcf",
    log: "{base}.hapcut2_extracthairs.log"
    shell:
        "extractHAIRS"
        " --10X 1"
        " --bam {input.bam}"
        " --VCF {input.vcf}"
        " --out {output.unlinked} 2> {log}"


rule hapcut2_linkfragments:
    """Link heterozygous variants together using barcode information"""
    output:
        linked = "{base}.calling.linked.txt"
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
        " --out {output.linked} &> {log}"


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
        " --outvcf 1 2> {log}"
        " && "
        "mv {output.phase}.phased.VCF {output.phased_vcf}"


rule symlink_reference_phased:
    output: "ground_truth.phased.vcf"
    shell:
        "ln -s {config[phasing_ground_truth]} {output}"


rule hapcut2_stats:
    """Get phasing statistics relative the ground truth. See https://github.com/vibansal/HapCUT2/tree/master/utilities
    for details. """
    output:
        stats = expand("{{base}}.phasing_stats.{ext}", ext=["txt", "tsv"])
    input:
        vcf1 = "{base}.calling.phased.vcf",
        vcf2 = "ground_truth.phased.vcf"
    params:
        stats_prefix = "{base}.phasing_stats"
    shell:
        "blr calculate_haplotype_statistics"
        " -v1 {input.vcf1}"
        " -v2 {input.vcf2}"
        " -o {params.stats_prefix}"


def get_haplotag_input(wildcards):
    inputfiles = {"bam": f"{wildcards.base}.calling.bam"}
    if config["reference_variants"] or config["haplotag_tool"] == "blr":
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
    Adds HP tag with haplotype (0 or 1) and PS tag with phase set information.
    """
    output:
        bam = "{base}.calling.phased.bam"
    input:
        unpack(get_haplotag_input)
    log: "{base}.haplotag.log"
    run:
        commands = {
            "whatshap":
                "whatshap haplotag"
                " {input.vcf}"
                " {input.bam}"
                " --linked-read-distance-cutoff 30000"
                " --reference {config[genome_reference]}",
            "blr":
                "blr phasebam"
                " --molecule-tag {config[molecule_tag]}"
                " --phase-set-tag {config[phase_set_tag]}"
                " --haplotype-tag {config[haplotype_tag]}"
                " {input.bam}"
                " {input.hapcut2_phase_file}"
        }
        command = commands[config["haplotag_tool"]]
        shell(command + " 2> {log} |"
            " tee {output.bam}"
            " |"
            " samtools index  - {output.bai}")


rule build_config
    """
    Builds a config file required for running NAIBR.
    """
    output:
        config = "naibr.config"
    log: "build_config.log"
    shell:
        "blr naibrconfig"
        " -o naibr.config"
        " 2> {log}"


rule lsv_calling:
    """
    Runs NAIBR for LSV calling.
    """
    output:
        results = "lsv_results"
    input:
        config = "naibr.config"
        bam = "{base}.calling.phased.bam"
        index = "{base}.calling.phased.bam.bai"
    log: "{base}.lsv_calling.log"
    shell:
        "python"
        " NAIBR.py"
        " {input.config}"
        " 2> {log}"
