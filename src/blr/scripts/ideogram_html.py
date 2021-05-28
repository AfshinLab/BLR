"""
Generate ideogram HTML code for MultiQC report integration with the Ideogram.js library
See: https://github.com/eweitz/ideogram
"""
from collections import namedtuple, defaultdict
from itertools import cycle
from pysam import VariantFile
import random


def parse_vcf_phase(vcf_file):
    """Generate phaseblocks from phased VCF"""
    blocks = defaultdict(list)
    chrom_blocks = defaultdict(list)
    with VariantFile(vcf_file) as open_vcf:
        sample_name = open_vcf.header.samples[0]

        prev_chrom = None
        for rec in open_vcf.fetch():
            sample = rec.samples[sample_name]

            if not sample.phased:
                continue

            chrom = rec.chrom
            ps = sample.get("PS")
            pos = rec.start + 1
            if not prev_chrom:
                prev_chrom = chrom

            # If new chromosome, add blocks to chrom_blocks and reset
            if chrom != prev_chrom:
                # Get start/stop for each block
                chrom_blocks[prev_chrom] = [(min(v), max(v)) for k, v in blocks.items()]
                blocks = defaultdict(list)
                prev_chrom = chrom

            if ps:
                # Aggregate all positions in phaseblock
                blocks[ps].append(pos)

    # Final
    chrom_blocks[prev_chrom] = [(min(v), max(v)) for k, v in blocks.items()]
    blocks.clear()
    return chrom_blocks


def parse_blocks(phased_vcf_file):
    """Parse over phaseblocks in VCF"""
    chrom_blocks = parse_vcf_phase(phased_vcf_file)
    Phaseblock = namedtuple("Phaseblock", ("chr", "start", "stop"))
    for chrom, blocks in chrom_blocks.items():
        blocks.sort(key=lambda x: x[0])
        for b in blocks:
            yield Phaseblock(chrom, *b)


def fromat_size(length: int) -> str:
    """Reformat phaseblock size to Mbp or kbp based on size"""
    if len(str(length)) > 5:
        return str(round(length/1_000_000, 1)) + " Mbp"
    elif len(str(length)) > 2:
        return str(round(length / 1_000, 1)) + " kbp"
    else:
        return str(length) + " bp"


#
# Extract phaseblocks from phased VCF
#
chr_stacks = defaultdict(list)
prev = 0
chrom = None
longest = (None, 0, 0)
for phaseblock in parse_blocks(snakemake.input.phased_vcf):  # noqa: F821
    if phaseblock:
        if phaseblock.chr != chrom:
            chrom = phaseblock.chr
            prev = 0

        pos, length = 0, 0
        if phaseblock.stop < prev:  # Phaseblock is completely overlapping previous
            continue
        elif phaseblock.start < prev:  # Phaseblock is partly overlapping previous
            pos, length = prev, int(phaseblock.stop - prev)
        else:
            pos, length = phaseblock.start, int(phaseblock.stop-phaseblock.start)

        chr_stacks[chrom].append((pos, length))

        # Find longest phaseblock
        if length > longest[2]:
            longest = (chrom, pos, length)

        prev = phaseblock.stop

#
# Generate HTML for ideogram
#
letters = "abcdefghijklmnopqrstuvwxyz"
unique_id = "".join(random.sample(letters, 10))

with open(snakemake.output.html, "w") as out:  # noqa: F821
    # Based on https://eweitz.github.io/ideogram/annotations-overlaid
    html = f"""
<!--
id: 'phaseblock-overview-{unique_id}'
section_name: 'Phaseblock overview'
description: "Ideogram of chromsomes with phaseblocks overlayed in alternating green and blue. \
Plot was generated using <a href='https://eweitz.github.io/ideogram/'>Ideogram.js</a>. Click on a \
chromsome to enlarge it. Click on the chromosome once again to return to the overview. Hover over block \
to get the phaseblock location along with it size. The longest phaseblock is highlighted in red."
-->

<div class="ideogram-{unique_id}">
</div>
<script type="text/javascript">
  var config = {{
    container: '.ideogram-{unique_id}',
    organism: 'human',
    assembly: '{snakemake.params.assembly}',
    chrHeight: 500,
    chrMargin: 1,
    annotations: {{
        "keys":["chr","name","start","length","color"],
        "annots":[
"""
    colors = [
        "rgba(63, 191, 63, 0.65)",  # Green
        "rgba(63, 127, 191, 0.65)"  # Blue
    ]
    longest_color = "rgba(191, 63, 63, 0.65)"  # Red
    # Separate entries for each chromosome
    for nr, (chrom, phaseblocks) in enumerate(chr_stacks.items()):
        html += "            "
        chr_id = chrom.replace("chr", "")
        html += f"{{'chr': '{chr_id}', 'annots':[\n"
        html += "                "
        blocks = []
        for (start, length), color in zip(phaseblocks, cycle(colors)):
            size = fromat_size(length)

            # Highlight longest phaseblock in different color
            if longest == (chrom, start, length):
                color = longest_color
                size += " (TOP)"

            blocks.append(f"['{chr_id}_{start}', '{size}', {start}, {length}, '{color}']")

        html += ",\n                ".join(blocks)
        html += "\n            ]},"
        if nr < len(chr_stacks) - 1:
            html += "\n"

    html += """
        ]},
    annotationsLayout: 'overlay',
    orientation: 'horizontal'
  };

  var ideogram = new Ideogram(config);
</script>
"""
    out.write(html)
