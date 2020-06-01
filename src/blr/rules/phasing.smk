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
    log: "{base}.calling.unlinked.txt.extracthairs.log"
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
    log: "{base}.calling.linked.txt.linkfragments.log"
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
        phased_vcf = "{base}.calling.phase.phased.VCF"
    input:
        linked = "{base}.calling.linked.txt",
        vcf = "{base}.phaseinput.vcf",
    log: "{base}.calling.phase.hapcut2.log"
    shell:
         "hapcut2"
         " --nf 1"
         " --fragments {input.linked}"
         " --vcf {input.vcf}"
         " --out {output.phase}"
         " --error_analysis_mode 1"
         " --outvcf 1 2> {log}"


rule symlink_reference_phased:
    output: "ground_truth.phased.vcf"
    shell: "ln -s {config[phasing_ground_truth]} {output}"


rule hapcut2_stats:
    """Get phasing statistics relative the ground truth. See https://github.com/vibansal/HapCUT2/tree/master/utilities
    for details. """
    output:
        stats = expand("{{base}}.phasing_stats.{ext}", ext=["txt", "tsv"])
    input:
         vcf1 = "{base}.calling.phase.phased.VCF",
         vcf2 = "ground_truth.phased.vcf"
    params:
        base = "{base}"
    shell:
         "blr calculate_haplotype_statistics"
         " -v1 {input.vcf1}"
         " -v2 {input.vcf2}"
         " -o {params.base}.phasing_stats"


rule compress_and_index_phased_vcf:
    "Compress and index VCF files."
    output:
        vcf = "{base}.calling.phase.phased.vcf.gz",
        index = "{base}.calling.phase.phased.vcf.gz.tbi"
    input:
        vcf = "{base}.calling.phase.phased.VCF"
    shell:
         "bgzip -c {input.vcf} > {output.vcf} && tabix -p vcf {output.vcf}"


def get_haplotag_input(wildcards):
    inputfiles = {"bam": f"{wildcards.base}.calling.bam"}
    if config["reference_variants"] or config["haplotag_tool"] == "blr":
        inputfiles.update({
            "hapcut2_phase_file": f"{wildcards.base}.calling.phase"
        })
    elif config["haplotag_tool"] == "whatshap":
        inputfiles.update({
            "vcf": f"{wildcards.base}.calling.phase.phased.vcf.gz",
            "vcf_index": f"{wildcards.base}.calling.phase.phased.vcf.gz.tbi",
        })
    return inputfiles


rule haplotag:
    """
    Transfer haplotype information from the phased VCF file to the bam file.
    Adds HP tag with haplotype (0 or 1) and PS tag with phase set information.
    """
    output:
        bam = "{base}.calling.phase.bam"
    input:
        unpack(get_haplotag_input)
    log: "{base}.haplotag.log"
    run:
        commands = {
            "whatshap":
                "whatshap haplotag"
                " {input.vcf}"
                " {input.bam}"
                " -o {output.bam}"
                " --linked-read-distance-cutoff 30000"
                " --reference {config[genome_reference]}",
            "blr":
                "blr phasebam"
                " --molecule-tag {config[molecule_tag]}"
                " --phase-set-tag {config[phase_set_tag]}"
                " --haplotype-tag {config[haplotype_tag]}"
                " {input.bam}"
                " {input.hapcut2_phase_file}"
                " -o {output.bam}"
        }

        command = commands[config["haplotag_tool"]].format(**locals(), **globals())

        shell("{command} 2> {log}")
