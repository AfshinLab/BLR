#!/usr/bin/env python
""" BLR MultiQC plugin module for general stats"""

from __future__ import print_function
from collections import OrderedDict
import logging
import pandas as pd

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
            name="Stats",
            target="stats",
            anchor="stats",
            info=" is collection of statistics from the different BLR commandline tools."
        )

        self.gather_stats_logs()

        n_phaseblock_reports = self.gather_phaseblock_data()
        if n_phaseblock_reports > 0:
            log.info("Found {} phaseblock reports".format(n_phaseblock_reports))

    def gather_stats_logs(self):
        # Find and load any input files for this module
        headers = dict()
        stats_data = dict()
        for f in self.find_log_files('stats', filehandles=True):
            tool_name = self.get_tool_name(f["f"])

            # If tool_name is None then there are no stats in the file --> skip.
            if not tool_name:
                continue

            if tool_name not in stats_data:
                stats_data[tool_name] = dict()
                headers[tool_name] = OrderedDict()

            sample_name = self.clean_s_name(f["fn"], f["root"])

            log.debug(f"Found report for tool {tool_name} with sample {sample_name}")

            if sample_name in stats_data[tool_name]:
                log.debug(f"Duplicate sample name found for tool {tool_name}! Overwriting: {sample_name}")

            stats_data[tool_name][sample_name] = dict()

            for parameter, value in self.parse(f["f"]):
                header_name = parameter.lower().replace(" ", "_")
                stats_data[tool_name][sample_name][header_name] = value

                headers[tool_name][header_name] = {
                    'title': parameter
                }

        # Nothing found - raise a UserWarning to tell MultiQC
        if len(stats_data) == 0:
            log.debug("Could not find any stats logs in {}".format(config.analysis_dir))

        log.info(f"Found {len(stats_data)} tools (Report per tool: "
                 f"{', '.join([tool + '=' + str(len(reps)) for tool, reps in stats_data.items()])})")

        # For each tool generat a separat statistics table for all found samples.
        for tool_name, data in stats_data.items():
            tool_name_title = tool_name.capitalize()

            # Write parsed report data to a file
            self.write_data_file(data, f"{tool_name}_stats")

            pconfig = {
                'id': 'blr_stats_table',
                'title': f"{tool_name_title} stats",
            }
            table_html = table.plot(data, headers[tool_name], pconfig)

            # Add a report section with table
            self.add_section(
                name=tool_name_title,
                description=f"Statistics table for data from BLR tool {tool_name}",
                helptext='''
                This longer description explains what exactly the numbers mean
                and supports markdown formatting. This means that we can do _this_:

                * Something important
                * Something else important
                * Best of all - some `code`

                Doesn't matter if this is copied from documentation - makes it
                easier for people to find quickly.
                ''',
                plot=table_html
            )

    def gather_phaseblock_data(self):
        data_lengths = dict()
        for f in self.find_log_files('stats/phaseblock_data', filehandles=True):
            sample_name = self.clean_s_name(f["fn"], f["root"]).replace(".phaseblock_data", "")
            sample_data = pd.read_csv(f["f"], sep="\t")
            data_lengths[sample_name] = sample_data["Length"].to_list()

        if len(data_lengths) == 0:
            log.debug("Could not find any phaseblock data in {}".format(config.analysis_dir))
            return 0

        data_lengths_binned = dict()
        binsize = 50000
        max_length = max([max(v) for v in data_lengths.values()])
        bins = range(0, max_length + binsize, binsize)
        for sample, values in data_lengths.items():
            _, weights = bin_sum(values, binsize=binsize, normalize=True)
            data_lengths_binned[sample] = {
                int(b / 1000): w for b, w in zip(bins, weights)  # bin per kbp
            }

        # Add longest phaseblock to general stats table
        general_stats_data = {
            name: {"longest_phaseblock": max(data) / 1_000_000} for name, data in data_lengths.items()  # Length in Mbp
        }
        general_stats_header = OrderedDict({
            "longest_phaseblock": {
                'title': 'Longest phaseblock',
                'description': 'Longest phaseblock created',
                'scale': 'Blues',
                'suffix': ' Mbp',
                'format': '{:,.3f}'
            }})

        self.general_stats_addcols(general_stats_data, general_stats_header)

        pconfig = {
            'id': 'phasingblock_lengths',
            'title': "Phaseblock lengths",
            'xlab': "Phaseblock length (kbp)",
            'ylab': 'Total DNA density',
            'yCeiling': 1,
            'tt_label': '{point.x} kbp: {point.y:.4f}',
        }
        plot_html = linegraph.plot(data_lengths_binned, pconfig)

        # Add a report section with plot
        self.add_section(
            name="Phaseblock lengths",
            description="Phaseblock lengths",
            plot=plot_html
        )

        # Make new dict with keys as strings for writable output.
        lengths_writable = dict()
        for sample, data in data_lengths_binned.items():
            lengths_writable[sample] = {str(k): v for k, v in data.items()}

        # Write parsed report data to a file
        self.write_data_file(lengths_writable, "stats_phaseblock_lengths")

        return len(data_lengths)

    @staticmethod
    def get_tool_name(file):
        """ Get the tools name by locating the line starting with 'STATS SUMMARY' which contains the tool name in the
         format e.g 'STATS SUMMARY - blr.cli.tool_name. Return tool name"""
        for line in file:
            if line.startswith("STATS SUMMARY"):
                return line.strip().split(".")[-1]

        return None

    @staticmethod
    def parse(file):
        """
        This generator yields key-value pairs for the data from the line following `---` until the next line
        staring with `===`.
        """
        collect = False
        for line in file:
            # Collect stats after first line starting with '---' and stop when starts with '==='
            if line.startswith("---"):
                collect = True
                continue
            elif line.startswith("===") and collect:
                break

            if collect:
                # Collect parameter and value
                parameter, value = list(filter(None, line.strip().split("  ")))
                value = value.strip().replace(",", "")

                yield parameter, float(value)
