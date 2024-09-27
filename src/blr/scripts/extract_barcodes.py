"""
Get BED with non-overlapping regions covered by SVs from VCF
"""
import argparse
from collections import defaultdict
import os
import re
import subprocess
import sys

from Bio import SeqIO
import gfapy
from gfapy.sequence import rc
import pysam


# source: https://github.com/anne-gcd/MTG-Link/blob/851e4c14034a8a11be371c8e87561c222b4702df/helpers.py#L41
class Gap:
    """
    Class defining a gap/target characterized by:
    - its ID
    - its length
    - its left flanking sequence's name
    - its right flanking sequence's name
    """

    def __init__(self, gap):
        self._identity = gap.gid
        self._length = gap.disp
        self._left = gap.sid1
        self._right = gap.sid2

    @property
    def identity(self):
        return self._identity

    @property
    def length(self):
        return self._length

    @property
    def left(self):
        return self._left

    @property
    def right(self):
        return self._right

    def __getattr__(self, attr):
        """If Python doesn't find the attribute "attr", it calls this method and print an alert"""
        print("WARNING: There is no attribute {} here !".format(attr))

    def __delattr_(self, attr):
        """We can't delete an attribute, we raise the exception AttributeError"""
        raise AttributeError("You can't delete attributes from this class")

    def label(self):
        """Method to label the gap/target"""
        if self._identity == "*":
            return str(self.left) + "_" + str(self.right)
        else:
            return str(self.identity)

    def info(self):
        """Method to get some information on the gap/target"""
        if self.identity == "*":
            print(
                "WORKING ON TARGET: between contigs {} & {}; length {}".format(
                    self.left, self.right, self.length
                )
            )
        else:
            print("WORKING ON TARGET: {}; length {}".format(self.identity, self.length))

    def __repr__(self):
        return "Target: id ({}), length ({}), left flanking seq ({}), right flanking seq ({})".format(
            self.identity, self.length, self.left, self.right
        )


# source: https://github.com/anne-gcd/MTG-Link/blob/851e4c14034a8a11be371c8e87561c222b4702df/helpers.py#L111
class Scaffold(Gap):
    """
    Class defining a scaffold characterized by:
    - the gap/target it is linked to
    - its name
    - its orientation
    - its length
    - the path of its sequence
    """

    def __init__(self, gap, scaffold, gfa_file):
        super().__init__(gap)
        self.gap = gap
        self.scaffold = scaffold
        self._name = scaffold.name
        self._orient = scaffold.orient
        self._slen = scaffold.line.slen
        self._seq_path = scaffold.line.UR
        self.gfa_file = gfa_file

    @property
    def name(self):
        return self._name

    @property
    def orient(self):
        return self._orient

    @property
    def slen(self):
        return self._slen

    @property
    def seq_path(self):
        return self._seq_path

    def __getattr__(self, attr):
        """If Python doesn't find the attribute "attr", it calls this method and print an alert"""
        print("WARNING: There is no attribute {} here !".format(attr))

    def __delattr_(self, attr):
        """We can't delete an attribute, we raise the exception AttributeError"""
        raise AttributeError("You can't delete attributes from this class")

    def sequence(self):
        """Method to get the sequence of the scaffold"""
        # if relative path
        if not str(self.seq_path).startswith("/"):
            seq_link = (
                str("/".join(str(self.gfa_file).split("/")[:-1]))
                + "/"
                + str(self.seq_path)
            )
        # if absolute path
        else:
            seq_link = self.seq_path

        # get the sequence of the scaffold
        for record in SeqIO.parse(seq_link, "fasta"):
            if re.match(self.name, record.id):
                return record.seq

    def chunk(self, c):
        """Method to get the region of the chunk/flank"""
        # ----------------------------------------------------
        # For gaps/targets into scaffolds' sequences
        # ----------------------------------------------------
        if ("-L" in self.name) or ("-R" in self.name):
            coordsOnScaffold = re.findall(r"[0-9]+-[0-9]+", str(self.name))[0]

            # if left scaffold
            if self.scaffold == self.left:
                start = int(str(coordsOnScaffold).split("-")[1]) - c
                end = int(str(coordsOnScaffold).split("-")[1])

            # if right scaffold
            elif self.scaffold == self.right:
                start = int(str(coordsOnScaffold).split("-")[0])
                end = int(str(coordsOnScaffold).split("-")[0]) + c

            # NB: The 'contig_name' should match the contig name on the BAM file
            contig_name = re.split(r"_[0-9]+-[0-9]+", str(self.name))[0]
            return str(contig_name) + ":" + str(max(0, start)) + "-" + str(end)

        # ----------------------------------------------------
        # For gaps/targets between scaffolds' sequences
        # ----------------------------------------------------
        else:
            # if left_fwd or right_rev
            if (self.orient == "+" and self.scaffold == self.left) or (
                self.orient == "-" and self.scaffold == self.right
            ):
                start = self.slen - c
                end = self.slen
            # if right_fwd or left_rev
            elif (self.orient == "+" and self.scaffold == self.right) or (
                self.orient == "-" and self.scaffold == self.left
            ):
                start = 0
                end = c
            return str(self.name) + ":" + str(start) + "-" + str(end)

    def __repr__(self):
        return "Scaffold: name ({}), orientation ({}), length ({}), sequence's file ({})".format(
            self.name, self.orient, self.slen, self.seq_path
        )


def extractBarcodesWithPysam(bam, gapLabel, region, barcodesOccurrencesDict):
    with pysam.AlignmentFile(bam, 'rb') as reader:
        for read in reader.fetch(region=region):
            try:
                barcode = read.get_tag("BX")
            except KeyError:
                continue
            barcodesOccurrencesDict[barcode] += 1

    return barcodesOccurrencesDict


# source: https://github.com/anne-gcd/MTG-Link/blob/851e4c14034a8a11be371c8e87561c222b4702df/barcodesExtraction.py#L102
def extractBarcodesFromChunkRegions(
    current_gap, gfaFile, bamFile, chunkSize, barcodesMinOcc
):
    """
    To extract the barcodes of reads mapping on chunk/flank regions.

    Args:
        - current_gap: str
            current gap/target identification
        - gfaFile: file
            GFA file containing the gaps' coordinates
        - bamFile: file
            indexed BAM file obtained after mapping the linked reads onto the draft assembly
        - chunkSize: int
            size of the chunk/flank region
        - barcodesMinOcc: int
            minimal occurrence of barcodes observed in the union set from the two flanking gap/target sequences

    Return:
        - unionBarcodesFile: file
            file containing the extracted barcodes of the union of both left and right gap/target flanking sequences
    """
    # ----------------------------------------------------
    # Pre-Processing
    # ----------------------------------------------------
    ## Create the object 'gap' from the class 'Gap'
    gap = Gap(current_gap)
    if not gap:
        print(
            "Unable to create the object 'gap' from the class 'Gap'.",
            file=sys.stderr,
        )
        sys.exit(1)

    # Get some information on the current gap/target we are working on.
    gap.info()
    gapLabel = gap.label()

    # Create two objects ('leftScaffold' and 'rightScaffold') from the class 'Scaffold'.
    leftScaffold = Scaffold(current_gap, gap.left, gfaFile)
    if not leftScaffold:
        print(
            "Unable to create the object 'leftScaffold' from the class 'Scaffold'.",
            file=sys.stderr,
        )
        sys.exit(1)

    rightScaffold = Scaffold(current_gap, gap.right, gfaFile)
    if not rightScaffold:
        print(
            "Unable to create the object 'rightScaffold' from the class 'Scaffold'.",
            file=sys.stderr,
        )
        sys.exit(1)

    chunk_L = chunkSize
    chunk_R = chunkSize

    # ----------------------------------------------------
    # Extract Barcodes
    # ----------------------------------------------------
    # Initiate a dictionary to count the occurences of each barcode extracted on the chunk/flank regions.
    barcodesOccurrencesDict = defaultdict(int)

    # Obtain the left barcodes extracted on the left region (left chunk/flank) and store the barcodes and their occurences in the dict 'barcodes_occ'.
    leftRegion = leftScaffold.chunk(chunk_L)
    if not leftRegion:
        print("Unable to obtain the left region (left flank).", file=sys.stderr)
        sys.exit(1)
    extractBarcodesWithPysam(
        bamFile, gapLabel, leftRegion, barcodesOccurrencesDict
    )

    # Obtain the right barcodes extracted on the right region (right chunk/flank) and store the barcodes and their occurences in the dict 'barcodes_occ'.
    rightRegion = rightScaffold.chunk(chunk_R)
    if not rightRegion:
        print("Unable to obtain the right region (right flank).", file=sys.stderr)
        sys.exit(1)
    extractBarcodesWithPysam(
        bamFile, gapLabel, rightRegion, barcodesOccurrencesDict
    )

    if len(barcodesOccurrencesDict) == 0:
        print(
            "Error while extracting the barcodes.", file=sys.stderr
        )
        sys.exit(1)

    # Do the union of the barcodes on both left and right regions (e.g. both left and right chunks/flanks).
    gfa_name = gfaFile.split("/")[-1]
    unionBarcodesFile = (
        f"{gfa_name}.{gapLabel}.g{gap.length}.flank{chunkSize}.occ{barcodesMinOcc}.bxu"
    )
    try:
        with open(unionBarcodesFile, "w") as unionBarcFile:
            ## Filter barcodes by the minimal occurrence of barcodes observed in the union set from the two flanking gap/target sequences ('barcodesMinOcc')
            for (barcode, occurrence) in barcodesOccurrencesDict.items():
                if occurrence >= barcodesMinOcc:
                    unionBarcFile.write(barcode + "\n")
    except IOError as err:
        print(
            f"Unable to open or write to the output 'unionBarcodesFile' {unionBarcodesFile}",
            file=sys.stderr,
        )
        raise err

def get_region(current_gap, gfaFile, flanksize, chunkSize):
    gap = Gap(current_gap)

    # Get some information on the current gap/target we are working on.
    gap.info()
    gapLabel = gap.label()

    # Create two objects ('leftScaffold' and 'rightScaffold') from the class 'Scaffold'.
    leftScaffold = Scaffold(current_gap, gap.left, gfaFile)
    rightScaffold = Scaffold(current_gap, gap.right, gfaFile)

    chunk_L = chunkSize
    chunk_R = chunkSize

    leftRegion = leftScaffold.chunk(chunk_L)
    rightRegion = rightScaffold.chunk(chunk_R)

    gfa_name = gfaFile.split("/")[-1]
    unionBarcodesFile = (
        f"{gfa_name}.{gapLabel}.g{gap.length}.flank{chunkSize}.occ{barcodesMinOcc}.bxu"
    )
    

def main(outdir, gfa, bam, flanksize, minbarcocc):
    gfa_path = os.path.abspath(gfa)
    bam_path = os.path.abspath(bam)
    print("GFA:", gfa_path)
    os.makedirs(outdir, exist_ok=True)
    print("OUTDIR:", outdir)
    os.chdir(outdir)
    print("CWD:", os.getcwd())
    open_gfa = gfapy.Gfa.from_file(gfa_path)
    for gap in open_gfa.gaps:
        extractBarcodesFromChunkRegions(
            gap, gfa_path, bam_path, flanksize, minbarcocc
        )


if __name__ == "__main__":
    if len(sys.argv) == 1:
        outdir = snakemake.output.dir  # noqa: F821
        bam = snakemake.input.bam  # noqa: F821
        gfa = snakemake.input.gfa  # noqa: F821
        flanksize = snakemake.params.flanksize  # noqa: F821
        minbarcocc = snakemake.params.minbarcocc  # noqa: F821
        log = snakemake.log[0]  # noqa: F821
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument("outdir", help="Output directory")
        parser.add_argument("--gfa", help="Input GFA file", required=True)
        parser.add_argument("--bam", help="Input BAM file", required=True)
        parser.add_argument(
            "--flanksize",
            type=int,
            default=10_000,
            help="Flank size for extracting barcodes",
        )
        parser.add_argument(
            "--minbarcocc",
            type=int,
            default=2,
            help="Minimum barcode occurence for keeping in region",
        )
        args = parser.parse_args()

        outdir = args.outdir
        gfa = args.gfa
        bam = args.bam
        flanksize = args.flanksize
        minbarcocc = args.minbarcocc
        log = f"{args.outdir}.log"

    # Write stdout to log file
    with open(log, "w") as sys.stdout:
        main(outdir, gfa, bam, flanksize, minbarcocc)
