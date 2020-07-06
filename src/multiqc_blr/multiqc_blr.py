#!/usr/bin/env python
""" MultiQC example plugin functions

We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""

from __future__ import print_function
import logging

from multiqc.utils import config
from blr import __version__

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.multiqc_blr_version = __version__


# Add default config options for the things that are used in MultiQC_NGI
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
                                      'contents_re': '^SETTINGS FOR:*',
                                      'num_lines': 1}})

    if 'hapcut2/phasing_stats' not in config.sp:
        # Current looking for file containing the string "switch rate:" on the first line.
        config.update_dict(config.sp,
                           {'hapcut2/phasing_stats': {'fn': '*.txt',
                                                      'contents_re': '^switch rate:*',
                                                      'num_lines': 1}})
    if 'hapcut2/phaseblocks' not in config.sp:
        # Currently looking for file containing the string "switch rate:" on the first line.
        config.update_dict(config.sp,
                           {'hapcut2/phaseblocks': {
                               'fn': '*.phase',
                               'contents_re': "^BLOCK:*",
                               'num_lines': 1}})
