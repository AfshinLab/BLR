#!/usr/bin/env python
""" BLR MultiQC plugin module for whatshap stats"""

from __future__ import print_function
import logging
import pandas as pd
from collections import OrderedDict, defaultdict

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
            name="Whatshap",
            target="WhatsHap",
            anchor="whatshap",
            info="""
             is a software for phasing genomic variants using DNA sequencing reads, also called read-based phasing or
             haplotype assembly. It is especially suitable for long reads, but works also well with short reads.
             """
        )

        table_data, snvs_phased_data = self.parse_stats()

        if table_data:
            table_headers = self.get_stats_table_headers()

            # Add a report section with table
            self.add_section(
                name="Phasing stats",
                description="Phasing statistics table",
                helptext='''
                Collection of phasing statistics related to variants and phaseblocks
                ''',
                plot=table.plot(table_data, table_headers, {
                    'id': 'whatshap_phasing_stats_table',
                    'title': "WhatsHap phasing stats",
                    'scale': False,
                    'share_key': False
                })
            )

            # Add N50 to general stats table
            general_stats_data = {
                sample: {"percent_SNVs_phased": data["percent_SNVs_phased"]} for sample, data in table_data.items()
            }

            general_stats_header = OrderedDict({
                "percent_SNVs_phased": {
                    'title': 'SNVs phased',
                    'description': 'Percentage of heterozygous SNVs that are phased',
                    'scale': 'Blues',
                    'suffix': '%',
                    'format': '{:,.3f}'
                }})

            self.general_stats_addcols(general_stats_data, general_stats_header)

        if snvs_phased_data:
            # Plot percent SNVs phased per chromosome
            self.add_section(
                name='Percent SNVs phased',
                description='Percent of heterozygous SNVs phased per chromosome. ALL shows the values for all '
                            'chromosomes combined',
                plot=linegraph.plot(snvs_phased_data, {
                    'title': 'Percent of heterozygous SNVs phased',
                    'xlab': 'Chromosome',
                    'ylab': '% SNVs phased',
                    'ymax': 100,
                    'ymin': 80,
                    'categories': True,
                    'tt_label': '<b>{point.x}</b>: {point.y:.3f}%',
                })
            )

        haplotag_data = self.parse_haplotag()
        if haplotag_data:
            table_headers = self.get_haplotag_table_headers()

            # Add a report section with table
            self.add_section(
                name="Haplotag info",
                description="Information about haplotype assingment",
                plot=table.plot(haplotag_data, table_headers, {
                    'id': 'whatshap_haplotag_table',
                    'title': "WhatsHap haplotag",
                    'scale': False,
                    'share_key': False
                })
            )

    def parse_stats(self):
        table_data = dict()
        snvs_phased_data = dict()
        for f in self.find_log_files('whatshap/stats', filehandles=True):
            s_name = self.clean_s_name(f["fn"], f["root"]).replace(".whatshap_stats", "")
            s_data = pd.read_csv(f["f"], sep="\t")

            # Add custom columns
            s_data['percent_SNVs_phased'] = 100 * s_data["phased_snvs"] / s_data["heterozygous_snvs"]
            s_data['percent_variants_phased'] = 100 * s_data["phased"] / s_data["heterozygous_variants"]

            # For multiple chromosome only keep chromsome ALL which is the aggregate data.
            if len(s_data) > 1:
                all_data = s_data[s_data["chromosome"] == "ALL"].drop(["#sample", "file_name", "chromosome"], axis=1)
            elif len(s_data) == 1:
                all_data = s_data.drop(["#sample", "file_name", "chromosome"], axis=1)
            else:
                continue

            table_data[s_name] = dict()
            for parameter, value in all_data.to_dict("records")[0].items():
                table_data[s_name][parameter] = value

            snvs_phased_data[s_name] = dict()
            for row in s_data.itertuples():
                snvs_phased_data[s_name][row.chromosome] = row.percent_SNVs_phased

        # Filter out samples to ignore
        table_data = self.ignore_samples(table_data)
        snvs_phased_data = self.ignore_samples(snvs_phased_data)

        if len(table_data) == 0:
            log.debug("Could not find any whatshap stats in {}".format(config.analysis_dir))
            return table_data, snvs_phased_data

        # Write parsed report data to a file
        self.write_data_file(table_data, "whatshap_stats")
        self.write_data_file(snvs_phased_data, "whatshap_stats_snvs_phased")

        return table_data, snvs_phased_data

    def parse_haplotag(self):
        data = defaultdict(dict)
        for f in self.find_log_files("whatshap/haplotag", filehandles=True):
            s_name = self.clean_s_name(f["fn"], f["root"]).replace(".haplotag", "")

            s_data = dict()
            collect_data = False
            for line in f["f"]:
                if line.strip() == "== SUMMARY ==":
                    collect_data = True
                    continue

                if collect_data and not line.startswith("haplotag"):
                    param, value = line.strip().split(":", maxsplit=1)
                    s_data[param] = int(value.strip())

            # Calculate percent of reads that were tagged.
            s_data["Percentage tagged"] = 100 * s_data['Alignments that could be tagged'] / s_data['Total alignments processed']
            data[s_name] = s_data

        data = self.ignore_samples(data)
        if len(data) == 0:
            log.debug("Could not find any whatshap haplotag data in {}".format(config.analysis_dir))
            return data

        self.write_data_file(data, "whatshap_haplotag")

        return data

    @staticmethod
    def get_stats_table_headers():
        headers = OrderedDict()
        headers['variants'] = {
            'title': 'Variants',
            'description': 'The total number of variants.',
            'format': '{:,}',
            'placement': 1
            }
        headers['phased'] = {
            'title': 'Phased variants',
            'description': 'The number of variants that are phased.',
            'format': '{:,}',
            'placement': 4
        }
        headers['unphased'] = {
            'title': 'Unphased variants',
            'description': 'The number of variants that are not phased',
            'format': '{:,}',
            'placement': 6,
            'hidden': True,
        }
        headers['singletons'] = {
            'title': 'Singletons',
            'description': 'The number of phaseblocks covering only a single variant.',
            'format': '{:,}',
            'placement': 7,
            'hidden': True,
        }
        headers['blocks'] = {
            'title': 'Phaseblocks',
            'description': 'The total number of phaseblocks',
            'format': '{:,}',
            'placement': 8,
        }
        headers['variant_per_block_median'] = {
            'title': 'Variants per block (median)',
            'description': 'The median number of variants covered by phaseblocks.',
            'format': '{:,.3f}',
            'placement': 9,
            'hidden': True,
        }
        headers['variant_per_block_avg'] = {
            'title': 'Variants per block (average)',
            'description': 'The average number of basepairs covered by phaseblocks.',
            'format': '{:,.3f}',
            'placement': 10,
            'hidden': True,
        }
        headers['variant_per_block_min'] = {
            'title': 'Variants per block (min)',
            'description': 'The minimum of variants covered by a phaseblock.',
            'format': '{:,}',
            'placement': 11,
            'hidden': True,
        }
        headers['variant_per_block_max'] = {
            'title': 'Variants per block (max)',
            'description': 'The maximum of variants covered by a phaseblock i.e. the shortest phaseblock.',
            'format': '{:,}',
            'placement': 12,
            'hidden': True,
        }
        headers['variant_per_block_sum'] = {
            'title': 'Variants per block (sum)',
            'description': 'The total sum of variants covered by phaseblocks.',
            'format': '{:,}',
            'placement': 13,
            'hidden': True,
        }
        headers['bp_per_block_median'] = {
            'title': 'Bp per block (median)',
            'description': 'The median number of basepairs covered by phaseblocks.',
            'format': '{:,.3f}',
            'placement': 14,
            'hidden': True,
        }
        headers['bp_per_block_avg'] = {
            'title': 'Bp per block (average)',
            'description': 'The average number of basepairs covered by phaseblocks.',
            'format': '{:,.3f}',
            'placement': 15,
            'hidden': True,
        }
        headers['bp_per_block_min'] = {
            'title': 'Bp per block (min)',
            'description': 'The minimum of basepairs covered by a phaseblock i.e. the shortest phaseblock.',
            'format': '{:,}',
            'placement': 16,
            'hidden': True,
        }
        headers['bp_per_block_max'] = {
            'title': 'Bp per block (max)',
            'description': 'The maximum of basepairs covered by a phaseblock i.e. the longest phaseblock.',
            'format': '{:,}',
            'placement': 17,
            'hidden': True,
        }
        headers['bp_per_block_sum'] = {
            'title': 'Bp per block (sum)',
            'description': 'The total sum of basepairs covered by phaseblocks',
            'format': '{:,}',
            'placement': 18,
            'hidden': True,
        }
        headers['heterozygous_variants'] = {
            'title': 'Heterozygous variants',
            'description': 'Number of heterozygous variants',
            'format': '{:,}',
            'placement': 2,
            'hidden': False,
        }
        headers['heterozygous_snvs'] = {
            'title': 'Heterozygous SNVs',
            'description': 'Number of heterozygous SNVs',
            'format': '{:,}',
            'placement': 3,
            'hidden': False,
        }
        headers['phased_snvs'] = {
            'title': 'Phased SNVs',
            'description': 'Number of phased SNVs',
            'format': '{:,}',
            'placement': 5,
            'hidden': False,
        }
        headers['block_n50'] = {
            'title': 'Phaseblock N50',
            'description': 'Phaseblock N50 related to genome length.',
            'format': '{:,.3f}',
            'hidden': True,
            'placement': 19,
        }

        # Custom headers added below
        headers['percent_SNVs_phased'] = {
            'title': 'SNVs phased',
            'description': 'Percentage of heterozygous SNVs that are phased.',
            'format': '{:,.3f}%',
            'hidden': False,
            'placement': 0,
        }
        headers['percent_variants_phased'] = {
            'title': 'Variants phased',
            'description': 'Percentage of heterozygous variants that are phased.',
            'format': '{:,.3f}%',
            'hidden': False,
            'placement': 1,
        }
        return headers

    @staticmethod
    def get_haplotag_table_headers():
        headers = OrderedDict()
        headers['Total alignments processed'] = {
            'title': 'Total alignments',
            'description': 'The total number of alignments processed.',
            'format': '{:,}',
            'placement': 2
            }
        headers['Alignments that could be tagged'] = {
            'title': 'Alignments tagged',
            'description': 'The number of alignments that could be assigned to a haplotype',
            'format': '{:,}',
            'placement': 3
            }
        headers['Alignments spanning multiple phase sets'] = {
            'title': 'Multiple PS',
            'description': 'The number of alignments spanning multiple phase sets.',
            'format': '{:,}',
            'placement': 4
            }
        headers["Percentage tagged"] = {
            'title': 'Tagged',
            'description': 'The percentage of alignments that were assigned to a haplotype.',
            'format': '{:.2f}',
            'suffix': "%",
            'placement': 1
            }
        return headers
