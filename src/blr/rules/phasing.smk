
variants = "variants.reference.vcf" if config["reference_variants"] else "variants.called.vcf"
if config["variant_caller"] == "gatk" and config["BQSR"]:
    bamfile_basename = "mapped.sorted.tag.bcmerge.mkdup.mol.filt.BQSR"
else:
    bamfile_basename = "mapped.sorted.tag.bcmerge.mkdup.mol.filt"

rule hapcut2_extracthairs:
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
         " --outvcf 1 2> {log}"


rule symlink_reference_phased:
    output: "ground_truth.phased.vcf"
    shell: "ln -s {config[phasing_ground_truth]} {output}"


rule hapcut2_stats:
    output:
        stats = "phasing_stats.txt"
    input:
         vcf1 = bamfile_basename + ".phase.phased.VCF",
         vcf2 = "ground_truth.phased.vcf"
    shell:
         "blr calculate_haplotype_statistics"
         " -v1 {input.vcf1}"
         " -v2 {input.vcf2}"
         " > {output.stats}"


rule compress_and_index_phased_vcf:
    "Compress and index VCF files"
    output:
        vcf = bamfile_basename + ".phase.phased.vcf.gz",
        index = bamfile_basename + ".phase.phased.vcf.gz.tbi"
    input:
        vcf = bamfile_basename + ".phase.phased.VCF"
    shell:
         "bgzip -c {input.vcf} > {output.vcf} && tabix -p vcf {output.vcf}"

rule haplotag:
    "Transfer haplotype information fron the phased VCF file to the bam file. "
    output:
        bam = bamfile_basename + ".phase.bam"
    input:
        bam = bamfile_basename + ".bam",
        vcf = bamfile_basename + ".phase.phased.vcf.gz",
        vcf_index = bamfile_basename + ".phase.phased.vcf.gz.tbi",
    log: "whatshap_haplotag.log"
    shell:
        "whatshap haplotag"
        " {input.vcf}"
        " {input.bam}"
        " -o {output.bam}"
        " --reference {config[genome_reference]}"
        " 2> {log}"