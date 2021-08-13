#!/usr/bin/env python
""" BLR MultiQC plugin module for general stats"""

from collections import OrderedDict, defaultdict
import logging
import pandas as pd

from multiqc import config
from multiqc.plots import table, linegraph
from multiqc.modules.base_module import BaseMultiqcModule

from multiqc_blr.utils import bin_sum, get_tail_x


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

        n_molecule_length_reports = self.gather_molecule_lengths()
        if n_molecule_length_reports > 0:
            log.info("Found {} molecule length reports".format(n_molecule_length_reports))

        n_sv_size_reports = self.gather_sv_sizes()
        if n_sv_size_reports > 0:
            log.info("Found {} SV size reports".format(n_sv_size_reports))

        n_stats = self.gather_stats()
        if n_stats > 0:
            log.info("Found {} stats reports".format(n_stats))

    def gather_stats_logs(self):
        # Find and load any input files for this module
        headers = dict()
        data = dict()
        for f in self.find_log_files('stats', filehandles=True):
            tool_name = self.get_tool_name(f["f"])

            # If tool_name is None then there are no stats in the file --> skip.
            if not tool_name:
                continue

            if tool_name not in data:
                data[tool_name] = dict()
                headers[tool_name] = OrderedDict()

            sample_name = self.clean_s_name(f["fn"], f["root"]).replace(f".{tool_name}", "")

            log.debug(f"Found report for tool {tool_name} with sample {sample_name}")

            if sample_name in data[tool_name]:
                log.debug(f"Duplicate sample name found for tool {tool_name}! Overwriting: {sample_name}")

            self.add_data_source(f)

            data[tool_name][sample_name] = dict()

            for parameter, value in self.parse(f["f"]):
                header_name = parameter.lower().replace(" ", "_")
                data[tool_name][sample_name][header_name] = value

                headers[tool_name][header_name] = {
                    'title': parameter
                }

            # Remove sample if no data
            if not data[tool_name][sample_name]:
                data[tool_name].pop(sample_name)

        # Filter out samples to ignore for each tool
        data = {tool: self.ignore_samples(data) for tool, data in data.items() if self.ignore_samples(data)}

        # Nothing found - raise a UserWarning to tell MultiQC
        if len(data) == 0:
            log.debug("Could not find any stats logs in {}".format(config.analysis_dir))

        log.info(f"Found {len(data)} tools (Report per tool: "
                 f"{', '.join([tool + '=' + str(len(reps)) for tool, reps in data.items()])})")

        # For each tool generat a separat statistics table for all found samples.
        for tool_name, tool_data in data.items():
            tool_name_title = tool_name.capitalize()

            # Write parsed report data to a file
            self.write_data_file(tool_data, f"{tool_name}_stats")

            pconfig = {
                'id': 'blr_stats_table',
                'title': f"{tool_name_title} stats",
            }
            table_html = table.plot(tool_data, headers[tool_name], pconfig)

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

            # Include select stats from "plot" tool in general stats table.
            if tool_name == "plot":
                general_stats_data = dict()
                for name, data in tool_data.items():
                    general_stats_data[name] = {
                        "n50_lpm": data["n50_reads_per_molecule"],
                        "mean_molecule_length_kbp": data["mean_molecule_length"] / 1000,
                        "median_molecule_length_kbp": data["median_molecule_length"] / 1000,
                        "dna_in_molecules_20_kbp_percent": data["dna_in_molecules_>20_kbp_(%)"],
                        "dna_in_molecules_100_kbp_percent": data["dna_in_molecules_>100_kbp_(%)"],
                        "nr_barcodes_final": data["barcodes_final"],
                        "median_molecule_count": data["median_molecule_count"]
                    }

                # Scale number of barcodes
                nr_barcodes_max = max(v.get("nr_barcodes_final", 0) for v in general_stats_data.values())
                barcode_multiplier = 0.000001 if nr_barcodes_max > 1_000_000 else 0.001
                barcode_suffix = " M" if nr_barcodes_max > 1_000_000 else " K"

                general_stats_header = OrderedDict({
                    "n50_lpm": {
                        'title': 'N50 LPM',
                        'description': 'N50 linked-reads per molecule',
                        'scale': 'OrRd',
                        'format': '{:,.0f}'
                    },
                    "mean_molecule_length_kbp": {
                        'title': 'Mean len',
                        'description': 'Mean molecule length in kbp',
                        'scale': 'PuBu',
                        'suffix': ' kbp',
                        'format': '{:,.1f}'
                    },
                    "median_molecule_length_kbp": {
                        'title': 'Median len',
                        'description': 'Median molecule length in kbp',
                        'scale': 'BuPu',
                        'suffix': ' kbp',
                        'format': '{:,.1f}'
                    },
                    "dna_in_molecules_20_kbp_percent": {
                        'title': '>20kbp',
                        'description': 'Percent of DNA in molecules longer than 20 kbp',
                        'scale': 'Oranges',
                        'suffix': '%',
                        'format': '{:.1f}'
                    },
                    "dna_in_molecules_100_kbp_percent": {
                        'title': '>100kbp',
                        'description': 'Percent of DNA in molecules longer than 100 kbp',
                        'scale': 'YlOrBr',
                        'suffix': '%',
                        'format': '{:.1f}'
                    },
                    "nr_barcodes_final": {
                        'title': '# Bc ',
                        'description': 'Number of barcodes in final data.',
                        'modify': lambda x: x * barcode_multiplier,
                        'suffix': barcode_suffix,
                        'scale': 'BuGn',
                        'format': '{:,.1f}'
                    },
                    "median_molecule_count": {
                        'title': ' # Mol',
                        'description': 'Median number of molecules per barcode',
                        'scale': 'OrRd',
                        'format': '{:,.0f}'
                    },
                })

                self.general_stats_addcols(general_stats_data, general_stats_header)

    def gather_phaseblock_data(self):
        data_lengths = dict()
        for f in self.find_log_files('stats/phaseblock_data', filehandles=True):
            sample_name = self.clean_s_name(f["fn"], f["root"]).replace(".phaseblock_data", "")

            if sample_name in data_lengths:
                log.debug("Duplicate sample name found! Overwriting: {}".format(sample_name))

            self.add_data_source(f)

            sample_data = pd.read_csv(f["f"], sep="\t")
            data_lengths[sample_name] = sample_data["Length"].to_list()

        # Filter out samples to ignore
        data_lengths = self.ignore_samples(data_lengths)

        if sum(map(len, data_lengths.values())) == 0:
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
            name: {"longest_phaseblock": max(data)} for name, data in data_lengths.items()
        }
        longest_phaseblock_overall = max(v["longest_phaseblock"] for v in general_stats_data.values())
        multiplier = 0.000001 if longest_phaseblock_overall > 1_000_000 else 0.001
        suffix = " Mbp" if longest_phaseblock_overall > 1_000_000 else " kbp"
        general_stats_header = OrderedDict({
            "longest_phaseblock": {
                'title': 'Top block',
                'description': 'Longest phaseblock created',
                'scale': 'YlGn',
                'modify': lambda x: x * multiplier,
                'format': '{:,.3f}',
                'suffix': suffix,
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

    def gather_molecule_lengths(self):
        data_lengths = dict()
        xmax = 100
        for f in self.find_log_files('stats/molecule_lengths', filehandles=True):
            sample_name = self.clean_s_name(f["fn"], f["root"]).replace(".molecule_lengths", "")

            if sample_name in data_lengths:
                log.debug("Duplicate sample name found! Overwriting: {}".format(sample_name))

            self.add_data_source(f)

            sample_data = pd.read_csv(f["f"], sep="\t")
            sample_data["LengthSumNorm"] = sample_data["LengthSum"] / sample_data["LengthSum"].sum()
            data_lengths[sample_name] = {int(row.Bin/1000): row.LengthSumNorm for row in sample_data.itertuples()}

            xmax = max(get_tail_x(data_lengths[sample_name]), xmax)

        # Filter out samples to ignore
        data_lengths = self.ignore_samples(data_lengths)

        if len(data_lengths) == 0:
            log.debug("Could not find any molecule lengths data in {}".format(config.analysis_dir))
            return 0

        pconfig = {
            'id': 'molecule_lengths',
            'title': "Stats: Molecule lengths",
            'xlab': "Molecule length (kbp)",
            'ylab': 'Total DNA density',
            'xmax': xmax,
            'tt_label': '{point.x} kbp: {point.y:.4f}',
        }
        plot_html = linegraph.plot(data_lengths, pconfig)

        # Add a report section with plot
        self.add_section(
            name="Molecule lengths",
            description="Molecule lengths binned in 1 kbp and their length summed for eached bin to get the density.",
            plot=plot_html
        )

        # Make new dict with keys as strings for writable output.
        lengths_writable = dict()
        for sample, data in data_lengths.items():
            lengths_writable[sample] = {str(k): v for k, v in data.items()}

        # Write parsed report data to a file
        self.write_data_file(lengths_writable, "stats_molecule_lengths")

        return len(data_lengths)

    def gather_sv_sizes(self):
        data = [{}, {}, {}, {}]

        for f in self.find_log_files('stats/sv_sizes', filehandles=True):
            sample_name = self.clean_s_name(f["fn"], f["root"]).replace(".sv_sizes", "")
            sample_df = pd.read_csv(f["f"], sep="\t")
            for d in data:
                d[sample_name] = {}
            for row in sample_df.itertuples():
                data[0][sample_name][row.Size] = row.DEL + row.INV + row.DUP
                data[1][sample_name][row.Size] = row.DEL
                data[2][sample_name][row.Size] = row.INV
                data[3][sample_name][row.Size] = row.DUP

        # Filter out samples to ignore
        data = [self.ignore_samples(d) for d in data]

        if len(data[0]) == 0:
            log.debug("Could not find any molecule lengths data in {}".format(config.analysis_dir))
            return 0

        pconfig = {
            'id': 'sv_sizes',
            'title': "Stats: SV size distribution",
            'xlab': "SV size range",
            'yMinRange': (0, 10),
            'categories': True,
            'data_labels': [
                {'name': 'Total', 'ylab': 'Count'},
                {'name': 'DEL', 'ylab': 'Count'},
                {'name': 'INV', 'ylab': 'Count'},
                {'name': 'DUP', 'ylab': 'Count'},
            ]
        }
        plot_html = linegraph.plot(data, pconfig)

        # Add a report section with plot
        self.add_section(
            name="SV size distribution",
            description="Size distrobution of called structural variants (SV).",
            plot=plot_html
        )

        general_stats_data = {name: {"svs": sum(total_data.values())} for name, total_data in data[0].items()}
        general_stats_header = OrderedDict({
            "svs": {
                'title': 'SVs',
                'description': 'Total number of detected structural variants',
                'scale': 'Blues',
                'format': '{:,}'
            },
        })

        self.general_stats_addcols(general_stats_data, general_stats_header)

        return len(data[0])

    def gather_stats(self):
        names = ["MB", "MC", "RB"]
        xmax = {"MB": 10, "RB": 20}
        data = {name: dict() for name in names}
        for f in self.find_log_files('stats/general_stats', filehandles=True):
            sample_name = self.clean_s_name(f["fn"], f["root"]).replace(".stats", "")

            if sample_name in data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(sample_name))

            self.add_data_source(f)

            sample_data = defaultdict(list)
            for line in f["f"]:
                if line.startswith("MB"):
                    _, mol_per_bc, count = line.split("\t")
                    sample_data["MB"].append((int(mol_per_bc), int(count)))

                if line.startswith("MC"):
                    _, coverage_bin, count = line.split("\t")
                    sample_data["MC"].append((float(coverage_bin), int(count)))

                if line.startswith("RB"):
                    _, reads_bin, count = line.split("\t")
                    sample_data["RB"].append((int(reads_bin), int(count)))

            total_count = sum(s[1] for s in sample_data["MB"])
            data["MB"][sample_name] = {
                mol_per_bc: 100*count/total_count for mol_per_bc, count in sample_data["MB"]
            }
            xmax["MB"] = max(max(get_tail_x(data["MB"][s], threshold=0.99) for s in data["MB"]), xmax["MB"])

            total_count = sum(s[1] for s in sample_data["MC"])
            data["MC"][sample_name] = {
                coverage_bin: 100 * count / total_count for coverage_bin, count in sample_data["MC"]
            }

            total_count = sum(s[1] for s in sample_data["RB"])
            data["RB"][sample_name] = {
                reads_bin: 100 * count / total_count for reads_bin, count in sample_data["RB"]
            }
            xmax["RB"] = max(max(get_tail_x(data["RB"][s], threshold=0.99) for s in data["RB"]), xmax["RB"])

        # Filter out samples to ignore
        data = {name: self.ignore_samples(d) for name, d in data.items()}
        if any(len(d) == 0 for d in data.values()):
            log.debug("Could not find any stats reports in {}".format(config.analysis_dir))
            return 0

        self.add_section(
            name="Molecules per barcode",
            description="Molecule per barcode",
            plot=linegraph.plot(
                data["MB"],
                {
                    'id': 'molecules_per_barcode',
                    'title': "Stats: Molecules per barcode",
                    'xlab': "Molecules per barcode",
                    'ylab': 'Fraction of total',
                    'xmax': xmax["MB"],
                    'yCeiling': 100,
                    'yLabelFormat': '{value}%',
                    'tt_label': '{point.x} molecules: {point.y:.1f}%',
                })
        )

        self.add_section(
            name="Molecule coverage",
            description="Molecule coverage",
            plot=linegraph.plot(
                data["MC"],
                {
                    'id': 'molecule_coverage',
                    'title': "Stats: Molecule coverage",
                    'xlab': "Molecules coverage",
                    'ylab': 'Fraction of total',
                    'yCeiling': 100,
                    'yLabelFormat': '{value}%',
                    'tt_label': 'Coverage {point.x}: {point.y:.1f}%',
                })
        )

        self.add_section(
            name="Reads per barcode",
            description="Graph showing reads per barcode. Reads are counted for each barcode and then binned and the "
                        "total nr of barcodes counted for each bin. The number for each bin relates to the lower bin "
                        "threshold.",
            plot=linegraph.plot(
                data["RB"],
                {
                    'id': 'reads_per_barcode',
                    'title': "Stats: Reads per barcode",
                    'xlab': "Read count bin",
                    'ylab': 'Fraction of total',
                    'xmax': xmax["RB"],
                    'yCeiling': 100,
                    'yLabelFormat': '{value}%',
                    'tt_label': 'Read count {point.x}: {point.y:.1f}%',
                })
        )

        return len(data.popitem()[1])

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
