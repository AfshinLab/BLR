#!/usr/bin/env python
""" BLR MultiQC plugin module for general stats"""

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

from multiqc_blr.utils import bin_sum

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')


class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Halt execution if we've disabled the plugin
        if config.kwargs.get('disable_plugin', True):
            return None

        # Initialise the parent module Class object
        super(MultiqcModule, self).__init__(
            name="HapCUT2",
            target="HapCUT2",
            href="https://github.com/vibansal/HapCUT2",
            anchor="hapcut2",
            info=" is a package for haplotype assembly from sequencing data."
        )
        n_phasing_stats_reports = self.gather_phasing_stats()
        if n_phasing_stats_reports > 0:
            log.info("Found {} phasing stats reports".format(n_phasing_stats_reports))

        n_phaseblock_reports = self.gather_phaseblocks()
        if n_phaseblock_reports > 0:
            log.info("Found {} phaseblock reports".format(n_phaseblock_reports))

    def gather_phasing_stats(self):
        # Create headers
        headers = OrderedDict()
        headers['switch rate'] = {
            'title': 'Switch rate',
            'description': 'switch errors as a fraction of possible positions for switch errors',
            'format': '{:,.7f}',
            'placement': 1
            }

        headers['mismatch rate'] = {
            'title': 'Mismatch rate',
            'description': 'mismatch errors as a fraction of possible positions for mismatch errors',
            'format': '{:,.7f}',
            'placement': 2
        }

        headers['flat rate'] = {
            'title': 'Flat rate',
            'description': 'flat errors as a fraction of possible positions for flat errors',
            'format': '{:,.7f}',
            'hidden': True,
        }

        headers['phased count'] = {
            'title': 'Phased count',
            'description': 'count of total SNVs phased in the test haplotype',
            'format': '{:,.0f}',
            'placement': 3
        }

        headers['AN50'] = {
            'title': 'AN50 (Mbp)',
            'description': 'the AN50 metric of haplotype completeness',
            'format': '{:,.3f}',
            'hidden': True
        }

        headers['N50'] = {
            'title': 'N50 (Mbp)',
            'description': 'the N50 metric of haplotype completeness',
            'format': '{:,.3f}',
            'placement': 4
        }

        headers['num snps max blk'] = {
            'title': 'SNPs in max blk',
            'description': 'the fraction of SNVs in the largest (most variants phased) block',
            'format': '{:,.0f}',
            'placement': 5
        }

        # Find and load any input files for this module
        phasing_data = dict()
        for f in self.find_log_files('hapcut2/phasing_stats', filehandles=True):
            sample_name = self.clean_s_name(f["fn"], f["root"])
            phasing_data[sample_name] = dict()

            for parameter, value in self.parse_phasing_stats(f["f"]):
                phasing_data[sample_name][parameter] = value

        if len(phasing_data) > 0:
            # Write parsed report data to a file
            self.write_data_file(phasing_data, "hapcut2_phasing_stats")

            pconfig = {
                'id': 'hapcut2_phasing_stats_table',
                'title': "HapCUT2 phasing stats",
                'scale': False,
                'share_key': False
            }
            table_html = table.plot(phasing_data, headers, pconfig)

            # Add a report section with table
            self.add_section(
                name="HapCUT2 phasing stats",
                description="Statistics table",
                helptext='''
                Description of statistics (taken from https://github.com/vibansal/HapCUT2/tree/master/utilities):
                ''',
                plot=table_html
            )

            # Add N50 to general stats table
            general_stats_data = {
                sample: {"N50_phaseblock": data["N50"]} for sample, data in phasing_data.items()
            }

            general_stats_header = OrderedDict({
                "N50_phaseblock": {
                    'title': 'N50 phaseblock',
                    'description': 'N50 statistic for phaseblock lengths',
                    'scale': 'Blues',
                    'suffix': 'Mbp',
                    'format': '{:,.3f}'
                }})

            self.general_stats_addcols(general_stats_data, general_stats_header)

        return len(phasing_data)

    def gather_phaseblocks(self):
        # Collect rawdata of lengths from file
        rawdata = dict()
        for f in self.find_log_files('hapcut2/phaseblocks', filehandles=True):
            sample_name = self.clean_s_name(f["fn"], f["root"])
            rawdata[sample_name] = list()
            for phaseblock in self.parse_phaseblocks(f["f"]):
                rawdata[sample_name].append(phaseblock["phaseblock_length"])

        # Generate bins relative to max phaseblock length and sum for each bin and sample to get plot data
        phaseblock_lengths = dict()
        if rawdata:
            binsize = 50000
            max_length = max([max(v) for v in rawdata.values()])
            bins = range(0, max_length + binsize, binsize)
            for sample, data in rawdata.items():
                _, weights = bin_sum(data, binsize=binsize, normalize=True)
                phaseblock_lengths[sample] = {
                    int(b / 1000): w for b, w in zip(bins, weights)  # bin per kbp
                }

        if phaseblock_lengths:
            pconfig = {
                'id': 'hapcut2_phasingblock_lengths',
                'title': "HapCUT2 phaseblock lengths",
                'xlab': "Phaseblock length (kbp)",
                'ylab': 'Total DNA density',
                'yCeiling': 1,
                'tt_label': '{point.x} kbp: {point.y:.4f}',
            }

            plot_html = linegraph.plot(phaseblock_lengths, pconfig)

            # Add a report section with plot
            self.add_section(
                name="HapCUT2 phaseblock lengths",
                description="Phaseblock lengths as reported by HapCUT2",
                plot=plot_html
            )

            # Make new dict with keys as strings for writable output.
            phaseblock_lengths_writable = dict()
            for sample, data in phaseblock_lengths.items():
                phaseblock_lengths_writable[sample] = {str(k): v for k, v in data.items()}

            # Write parsed report data to a file
            self.write_data_file(phaseblock_lengths_writable, "hapcut2_phaseblock_lengths")

        return len(phaseblock_lengths)

    @staticmethod
    def parse_phasing_stats(file):
        """
        This generator yields key-value pairs for the data from the line following `---` until the next line
        staring with `===`.
        """
        for line in file:
            # Collect parameter and value
            parameter, value = line.strip().split(":", maxsplit=1)
            value = float(value.strip())

            # Make N50 and AN50 stats per Mbp instead of bp.
            if parameter in ["N50", "AN50"]:
                value /= 1_000_000

            yield parameter, value

    @staticmethod
    def parse_phaseblocks(file):
        """
        Format  description from https://github.com/vibansal/HapCUT2/blob/master/outputformat.md

        Example entry.
        ```
        BLOCK: offset: 6 len: 4 phased: 2 SPAN: 23514 fragments 1
        ```

        offset: <SNV offset>
        len: <SNV span of block>
        phased: <# SNVs phased>
        SPAN: <base pair span of block>
        fragments <# of fragments in block>

        """
        for line in file:
            if line.startswith("BLOCK:"):
                contents = line.split()
                yield {
                    "snv_span": int(contents[4]),
                    "phased_snvs": int(contents[6]),
                    "phaseblock_length": int(contents[8]),
                    "fragments": int(contents[10])
                }
            else:
                continue
