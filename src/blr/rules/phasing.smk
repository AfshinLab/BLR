"""
Rules related to phasing of variants (called or reference set)
"""

variants = "variants.reference.vcf" if config["reference_variants"] else "variants.called.vcf"
if config["variant_caller"] == "gatk" and config["BQSR"]:
    bamfile_basename = "mapped.sorted.tag.bcmerge.mkdup.mol.filt.BQSR"
else:
    bamfile_basename = "mapped.sorted.tag.bcmerge.mkdup.mol.filt"


rule hapcut2_extracthairs:
    """Extract heterozygous variants covered by alignments in BAM"""
    output:
        unlinked = bamfile_basename + ".unlinked.txt"
    input:
        bam = bamfile_basename + ".bam",
        vcf = variants
    log: "hapcut2_extracthairs.log"
    shell:
         "extractHAIRS"
         " --10X 1"
         " --bam {input.bam}"
         " --VCF {input.vcf}"
         " --out {output.unlinked} 2> {log}"


rule hapcut2_linkfragments:
    """Link heterozygous variants together using barcode information"""
    output:
        linked = bamfile_basename + ".linked.txt"
    input:
        bam = bamfile_basename + ".bam",
        bai = bamfile_basename + ".bam.bai",
        vcf = variants,
        unlinked = bamfile_basename + ".unlinked.txt"
    log: "hapcut2_linkfragments.log"
    shell:
         "LinkFragments.py"
         " --bam {input.bam}"
         " -v {input.vcf}"
         " --fragments {input.unlinked}"
         " --out {output.linked} &> {log}"


rule hapcut2_phasing:
    """Phase heterozygous varinats using HapCUT2. Output phased VCF"""
    output:
        phase = bamfile_basename + ".phase",
        phased_vcf = bamfile_basename + ".phase.phased.VCF"
    input:
        linked = bamfile_basename + ".linked.txt",
        vcf = variants
    log: "hapcut2_phasing.log"
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


phasing_stats_prefix = "phasing_stats"

rule hapcut2_stats:
    """Get phasing statistics relative the ground truth. See https://github.com/vibansal/HapCUT2/tree/master/utilities
    for details. """
    output:
        stats = expand(f"{phasing_stats_prefix}.{{ext}}", ext=["txt", "tsv"])
    input:
         vcf1 = bamfile_basename + ".phase.phased.VCF",
         vcf2 = "ground_truth.phased.vcf"
    shell:
         "blr calculate_haplotype_statistics"
         " -v1 {input.vcf1}"
         " -v2 {input.vcf2}"
         " -o {phasing_stats_prefix}"


rule compress_and_index_phased_vcf:
    "Compress and index VCF files."
    output:
        vcf = bamfile_basename + ".phase.phased.vcf.gz",
        index = bamfile_basename + ".phase.phased.vcf.gz.tbi"
    input:
        vcf = bamfile_basename + ".phase.phased.VCF"
    shell:
         "bgzip -c {input.vcf} > {output.vcf} && tabix -p vcf {output.vcf}"


def get_haplotag_input(wildcards):
    inputfiles = {"bam": bamfile_basename + ".bam"}
    if config["reference_variants"] or config["haplotag_tool"] == "blr":
        inputfiles.update({
            "hapcut2_phase_file": bamfile_basename + ".phase"
        })
    elif config["haplotag_tool"] == "whatshap":
        inputfiles.update({
            "vcf": bamfile_basename + ".phase.phased.vcf.gz",
            "vcf_index": bamfile_basename + ".phase.phased.vcf.gz.tbi",
        })
    return inputfiles


rule haplotag:
    """Transfer haplotype information fron the phased VCF file to the bam file. Adds HP tag with haplotype (0 or 1) and 
    PS tag with phase set information. 
    """
    output:
        bam = bamfile_basename + ".phase.bam"
    input:
        unpack(get_haplotag_input)
    log: "haplotag.log"
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
