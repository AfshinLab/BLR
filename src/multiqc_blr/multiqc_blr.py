#!/usr/bin/env python
""" MultiQC example plugin functions

We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""

import logging

from multiqc.utils import config
from blr import __version__

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_blr_version = __version__


# Add default config options for the things that are used in multiqc_blr
def before_config():
    # Use blr template by default
    config.template = "blr"


def execution_start():
    """ Code to execute after the config files and
    command line flags have been parsedself.

    This setuptools hook is the earliest that will be able
    to use custom command line flags.
    """

    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_plugin', True):
        return None

    log.info("Running MultiQC BLR Plugin v{}".format(config.multiqc_blr_version))

    # Add to the main MultiQC config object.
    # User config files have already been loaded at this point
    #   so we check whether the value is already set. This is to avoid
    #   clobbering values that have been customised by users.

    # Add to the search patterns used by modules
    if 'stats' not in config.sp:
        # Current looking for file ending with ".log" and having the content "SETTINGS FOR:" on the first line
        config.update_dict(config.sp,
                           {'stats': {'fn': '*.log',
                                      'contents': 'SETTINGS FOR:',
                                      'num_lines': 1,
                                      'max_filesize': 16384}})

    if 'stats/phaseblock_data' not in config.sp:
        config.update_dict(config.sp,
                           {'stats/phaseblock_data': {'fn': '*.phaseblock_data.tsv'}})

    if 'stats/molecule_lengths' not in config.sp:
        config.update_dict(config.sp,
                           {'stats/molecule_lengths': {'fn': '*.molecule_lengths.tsv'}})

    if 'stats/sv_sizes' not in config.sp:
        config.update_dict(config.sp,
                           {'stats/sv_sizes': {'fn': '*.sv_sizes.tsv',
                                               'contents': 'Size	DEL	INV	DUP',
                                               'num_lines': 1}})

    if 'stats/molecule_stats' not in config.sp:
        config.update_dict(config.sp,
                           {'stats/molecule_stats': {'fn': '*.molecule_stats.txt',
                                                     'contents': '# Stats compiled from molecule_stats.py',
                                                     'num_lines': 1}})

    if 'stats/barcode_stats' not in config.sp:
        config.update_dict(config.sp,
                           {'stats/barcode_stats': {'fn': '*.barcode_stats.txt',
                                                    'contents': '# Stats compiled from barcode_stats.py',
                                                    'num_lines': 1}})

    if 'hapcut2/phasing_stats' not in config.sp:
        # Current looking for file containing the string "switch rate:" on the first line.
        config.update_dict(config.sp,
                           {'hapcut2/phasing_stats': {'fn': '*.txt',
                                                      'contents': 'switch rate:',
                                                      'num_lines': 2}})

    if "whatshap/stats" not in config.sp:
        config.update_dict(config.sp,
                           {'whatshap/stats': {'fn': '*.tsv',
                                               'contents': '#sample	chromosome	file_name	variants	phased	unphased	singletons	blocks	variant_per_block_median	variant_per_block_avg	variant_per_block_min	variant_per_block_max	variant_per_block_sum	bp_per_block_median	bp_per_block_avg	bp_per_block_min	bp_per_block_max	bp_per_block_sum	heterozygous_variants	heterozygous_snvs	phased_snvs	block_n50',  # noqa: E501
                                               'num_lines': 1}})

    if "whatshap/haplotag" not in config.sp:
        config.update_dict(config.sp,
                           {'whatshap/haplotag': {'fn': '*haplotag.log',
                                                  'contents': 'haplotag - total processing time:'}})
