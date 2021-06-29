#!/usr/bin/env python
""" BLR MultiQC plugin module for general stats"""

from collections import OrderedDict
import logging

from multiqc import config
from multiqc.plots import table, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

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
        n_phasing_stats_reports, n_phasing_stats_reports_per_chrom = self.gather_phasing_stats()
        if n_phasing_stats_reports > 0:
            log.info("Found {} phasing stats reports".format(n_phasing_stats_reports))

        if n_phasing_stats_reports_per_chrom > 0:
            log.info("Found {} phasing stats reports with chromsome breakdown".format(
                n_phasing_stats_reports_per_chrom
                ))

    def gather_phasing_stats(self):
        # Create headers
        headers = OrderedDict()
        headers['switch rate'] = {
            'title': 'Switch rate',
            'description': 'switch errors as a fraction of possible positions for switch errors',
            'format': '{:,.7f}',
            'scale': 'Blues',
            'placement': 1
            }

        headers['mismatch rate'] = {
            'title': 'Mismatch rate',
            'description': 'mismatch errors as a fraction of possible positions for mismatch errors',
            'format': '{:,.7f}',
            'scale': 'Blues',
            'placement': 2
        }

        headers['flat rate'] = {
            'title': 'Flat rate',
            'description': 'flat errors as a fraction of possible positions for flat errors',
            'format': '{:,.7f}',
            'scale': 'Blues',
            'hidden': True,
        }

        headers['phased count'] = {
            'title': 'Phased count',
            'description': 'count of total SNVs phased in the test haplotype',
            'format': '{:,.0f}',
            'scale': 'Blues',
            'placement': 3
        }

        headers['AN50'] = {
            'title': 'AN50 (Mbp)',
            'description': 'the AN50 metric of haplotype completeness',
            'format': '{:,.3f}',
            'scale': 'Blues',
            'hidden': True
        }

        headers['N50'] = {
            'title': 'N50 (Mbp)',
            'description': 'the N50 metric of haplotype completeness',
            'format': '{:,.3f}',
            'scale': 'Blues',
            'placement': 4
        }

        headers['num snps max blk'] = {
            'title': 'SNPs in max blk',
            'description': 'the fraction of SNVs in the largest (most variants phased) block',
            'format': '{:,.0f}',
            'scale': 'Blues',
            'placement': 5
        }

        # Find and load any input files for this module
        phasing_data = dict()
        phasing_data_per_chrom = [{} for h in headers]
        param_to_index = {param: index for index, param in enumerate(headers)}
        for f in self.find_log_files('hapcut2/phasing_stats', filehandles=True):
            sample_name = self.clean_s_name(f["fn"], f["root"]).replace(".phasing_stats", "")

            if sample_name in phasing_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(sample_name))

            self.add_data_source(f)

            phasing_data[sample_name] = dict()
            for d in phasing_data_per_chrom:
                d[sample_name] = {}

            for chrom, parameter, value in self.parse_phasing_stats(f["f"]):
                if chrom == "All":
                    phasing_data[sample_name][parameter] = value
                else:
                    index = param_to_index[parameter]
                    phasing_data_per_chrom[index][sample_name][chrom] = value

        # Filter out samples to ignore
        phasing_data = self.ignore_samples(phasing_data)

        # Skip if no data
        if not phasing_data.values():
            return 0, 0

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
                'title': 'N50 block',
                'description': 'N50 statistic for phaseblock lengths',
                'scale': 'Blues',
                'suffix': ' Mbp',
                'format': '{:,.3f}'
            }})

        self.general_stats_addcols(general_stats_data, general_stats_header)

        # Return if not any per-chrom statistics
        nr_stats_per_chrom = sum(data != {} for sample, data in phasing_data_per_chrom[0].items())
        if nr_stats_per_chrom == 0:
            return len(phasing_data), 0

        # TODO Write data per chrom to file

        pconfig_per_chrom = {
            'id': 'hapcut2_phasing_stats_per_chrom_graph',
            'title': "HapCUT2: Phasing stats per chromosome",
            'xlab': "Chromosome",
            'categories': True,
            'tt_decimals': 4,
            "ymin": 0,
            'data_labels': [
                {'name': label["title"], 'ylab': label["title"]} for label in headers.values()
            ]
        }
        plot_html = linegraph.plot(phasing_data_per_chrom, pconfig_per_chrom)

        # Add a report section with plot
        self.add_section(
            name="Phasing stats per chromsomes",
            description="Breakdown of phasing stats for each chromosome.",
            plot=plot_html
        )

        return len(phasing_data), nr_stats_per_chrom

    @staticmethod
    def parse_phasing_stats(file):
        """
        This generator yields chromsome-key-value tuples for the data. "All" is used for the total.
        """
        chrom = "All"
        for line in file:
            # Search for lines of fromat `---- chrA -----` to get chromsome, else use default.
            if line.startswith("-"):
                chrom = line.split(" ", maxsplit=2)[1]
                continue

            # Collect parameter and value
            parameter, value = line.strip().split(":", maxsplit=1)

            if value.strip() == "n/a":
                continue

            value = float(value.strip())

            # Make N50 and AN50 stats per Mbp instead of bp.
            if parameter in ["N50", "AN50"]:
                value /= 1_000_000

            yield chrom, parameter, value
