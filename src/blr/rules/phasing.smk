"""
Rules related to phasing of variants (called or reference set)
"""


rule hapcut2_extracthairs:
    """Extract heterozygous variants covered by alignments in BAM"""
    output:
        unlinked = temp("{base}.calling.unlinked.txt")
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
        " --outvcf 1 2> {log}"
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
        vcf2 = f" -v2 {config['phasing_ground_truth']}" if config['phasing_ground_truth'] else ""
    shell:
        "blr calculate_haplotype_statistics"
        " -v1 {input.vcf1}"
        " {params.vcf2}"
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
    Adds HP tag with haplotype (0 or 1) and PS tag with phase set information.
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
                " -o {output.bam}"
                " --linked-read-distance-cutoff {config[window_size]}"
                " --reference {config[genome_reference]}"
                " {ignore_readgroups}",
            "blr":
                "blr phasebam"
                " --molecule-tag {config[molecule_tag]}"
                " --phase-set-tag {config[phase_set_tag]}"
                " --haplotype-tag {config[haplotype_tag]}"
                " --min-mapq {config[min_mapq]}"
                " {input.bam}"
                " {input.hapcut2_phase_file}"
                " -o {output.bam}"
        }
        command = commands[config["haplotag_tool"]]
        shell(command + " 2> {log}")
