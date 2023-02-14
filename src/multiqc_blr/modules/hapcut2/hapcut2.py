#!/usr/bin/env python
""" BLR MultiQC plugin module for general stats"""

from collections import OrderedDict
from itertools import cycle
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

        n_phasing_plots = self.gather_phasing_plots()
        if n_phasing_plots > 0:
            log.info("Found {} phasing plot data reports".format(
                n_phasing_plots
                ))

    def gather_phasing_plots(self):
        stat_names = [
            "NX",
            "ANX",
            "NGX",
            "QAN",
            "QNX",
            "QNG",
        ]
        data = {name: dict() for name in stat_names}
        phasing_data_per_plot = [{} for n in stat_names]
        stat_to_index = {name: index for index, name in enumerate(stat_names)}
        for f in self.find_log_files('hapcut2/phasing_plots', filehandles=True):
            sample_name = self.clean_s_name(f["fn"], f["root"]).replace(".phasing_stats.plot", "")

            if sample_name in data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(sample_name))

            self.add_data_source(f)
            for d in phasing_data_per_plot:
                d[sample_name] = {}

            for line in f["f"]:
                if line.startswith("#"):
                    continue

                stat, x, y = line.split("\t")
                assert stat in stat_names
                index = stat_to_index[stat]
                phasing_data_per_plot[index][sample_name][int(x)] = float(y)

        # Filter out samples to ignore
        phasing_data_per_plot = [
            {name: self.ignore_samples(d) for name, d in data.items()} for data in phasing_data_per_plot
        ]

        if all(len(d) == 0 for d in phasing_data_per_plot[0].values()):
            log.debug("Could not find any summary reports in {}".format(config.analysis_dir))
            return 0

        phasing_data_per_plot_final = []
        stat_names_final = []
        for name in stat_names:
            index = stat_to_index[name]
            data = phasing_data_per_plot[index]
            # Skip stats without data
            if not any(len(d) > 0 for d in data.values()):
                continue

            # Scale values and set unit
            max_y = max(max(d.values()) for d in data.values())
            multiplier = 0.000001 if max_y > 1_000_000 else 0.001
            unit = "Mbp" if max_y > 1_000_000 else "kbp"
            data = {sample: {x: y*multiplier for x, y in sample_data.items()} for sample, sample_data in data.items()}
            phasing_data_per_plot_final.append(data)
            stat_names_final.append({
                "name": name,
                "label": f"{name} [{unit}]",
                "unit": unit,
            })

        pconfig_per_plot = {
            'id': 'hapcut2_phasing_plot',
            'title': "HapCUT2: Phasing plots",
            'xlab': "X [%]",
            "ymin": 0,
            'data_labels': [
                {'name': name["name"], 'ylab': name["label"]} for name in stat_names_final
            ]
        }

        plot_html = linegraph.plot(phasing_data_per_plot_final, pconfig_per_plot)

        # Add a report section with plot
        self.add_section(
            name="Phasing plots",
            description="Plots showing countious phasing stats various levels of coverage.",
            helptext='''
            # Description of plots

            The plots display continous phasing statistics for more accurate comparison
            samples. All plot a fashioned on the so call *Nx-curve* for plotting phasing
            contiguity. Here we can define a value Nx for which phase blocks no shorter
            than Nx covers x% of the total phase genome length. We then plot Nx
            as a function of x, with x ranging from 0 to 100.

            Instead of using comparing phase block length to the total phased genome length
            other phasing statistics can be evaluated for differen insights.

            **NX**
            Relates the phase block length to the total length on the phased assembly.
            Compare to the `N50` (value at x=50) and `auN` (area under curve) metrics above.

            **ANX**
            Relates the corrected phase block length to the total number of phased variants.
            Phase block length corrected to proportion of phased variants with in the block
            range. Compare to the `AN50` (value at x=50) metric above.

            **NGX**
            Relates the phase block length to the total genome length. Compare to the `NG50`
            (value at x=50) and `auNG` (area under curve) metrics above.

            **QNX, QAN, QNG**
            Similar to plots detailed above but the phase block length have now been split at
            switch locations.
            ''',
            plot=plot_html
        )

        return sum(len(d) > 0 for d in phasing_data_per_plot_final[0].values())

    def gather_phasing_stats(self):
        # Create headers
        headers = OrderedDict()
        headers['switch rate'] = {
            'title': 'Switch rate',
            'description': 'Switch errors (aka long-switch error) as a fractions of possible switch positions',
            'format': '{:,.7f}',
            'placement': 1
            }

        headers['mismatch rate'] = {
            'title': 'Mismatch rate',
            'description': 'Mismatch errors (aka short-switch errors or point errors) as a fraction of possible '
                           'mismatch positions',
            'format': '{:,.7f}',
            'placement': 2
        }

        headers['flat rate'] = {
            'title': 'Flat rate',
            'description': 'Flat errors (aka Hamming errors) as a fraction of possible flat error positions.',
            'format': '{:,.7f}',
            'hidden': True,
        }

        headers['phased count'] = {
            'title': 'Phased count',
            'description': 'Count of total variants phased',
            'format': '{:,.0f}',
            'placement': 3
        }

        headers['AN50'] = {
            'title': 'AN50 block',
            'description': 'AN50 metric for haplotype contiguity.',
            'format': '{:,.3f}',
            'hidden': True
        }

        headers['N50'] = {
            'title': 'N50 block',
            'description': 'N50 metric for haplotype contiguity.',
            'format': '{:,.3f}',
            'placement': 4
        }

        headers['num snps max blk'] = {
            'title': 'Variants in max blk',
            'description': 'the fraction of variants in the largest (most variants phased) block',
            'format': '{:,.0f}',
            'hidden': True,
        }

        headers['auN'] = {
            'title': 'auN block',
            'description': 'auN metric for haplotype contiguity.',
            'format': '{:,.3f}',
            'placement': 5
        }

        headers['NG50'] = {
            'title': 'NG50 block',
            'description': 'NG50 metric for haplotype contiguity. Similar to N50 but relative the genome length.',
            'format': '{:,.3f}',
            'placement': 6
        }

        headers['auNG'] = {
            'title': 'auNG block',
            'description': 'auNG metric for haplotype contiguity. Similar to auN but for NGx curve.',
            'format': '{:,.3f}',
            'hidden': True
        }

        headers['switch count'] = {
            'title': 'Switch count',
            'description': 'Switch error counts.',
            'format': '{:,.0f}',
            'hidden': True
        }

        headers['switch positions'] = {
            'title': 'Switch positions',
            'description': 'The number of positions where switch errors are assayed.',
            'format': '{:,.0f}',
            'hidden': True
        }

        headers['mismatch count'] = {
            'title': 'Mismatch count',
            'description': 'Mismatch error counts.',
            'format': '{:,.0f}',
            'hidden': True
        }

        headers['mismatch positions'] = {
            'title': 'Mismatch positions',
            'description': 'The number of positions where mismatch errors are assayed.',
            'format': '{:,.0f}',
            'hidden': True
        }

        headers['flat count'] = {
            'title': 'Flat count',
            'description': 'Flat  error counts.',
            'format': '{:,.0f}',
            'hidden': True
        }

        headers['flat positions'] = {
            'title': 'Flat positions',
            'description': 'The number of positions where flat errors are assayed.',
            'format': '{:,.0f}',
            'hidden': True
        }

        headers['QAN50'] = {
            'title': 'QAN50 block',
            'description': 'Similar to AN50 but each phase block has been split at switch and mismatch locations.',
            'format': '{:,.3f}',
            'hidden': True
        }

        headers['QN50'] = {
            'title': 'QN50 block',
            'description': 'Similar to N50 but each phase block has been split at switch and mismatch locations.',
            'format': '{:,.3f}',
            'hidden': True
        }

        headers['auQN'] = {
            'title': 'auQN block',
            'description': 'Similar to auN but each phase block has been split at switch and mismatch locations.',
            'format': '{:,.3f}',
            'hidden': True
        }

        headers['QNG50'] = {
            'title': 'QNG50 block',
            'description': 'Similar to NG50 but each phase block has been split at switch and mismatch locations.',
            'format': '{:,.3f}',
            'hidden': True
        }

        headers['auQNG'] = {
            'title': 'auQNG block',
            'description': 'Similar to auNG but each phase block has been split at switch and mismatch locations.',
            'format': '{:,.3f}',
            'hidden': True
        }

        headers['phased count ref'] = {
            'title': 'Phased count ref',
            'description': 'Count of total variants phased in reference haplotype',
            'format': '{:,.0f}',
            'hidden': True
        }

        headers['phased rate asm'] = {
            'title': '% Phased asm',
            'description': 'Percentage of phased variants in the assembly haplotype that are phased in the reference.',
            'format': '{:,.3%}',
            'hidden': True
        }

        headers['phased rate ref'] = {
            'title': '% Phased ref',
            'description': 'Percentage of phased variant sin the reference haplotype that are phased in the assembly.',
            'format': '{:,.3%}',
            'hidden': True
        }

        # Colorbrewer2 scales from
        # https://github.com/axismaps/colorbrewer/blob/9a37cbbfe7cde61c060c68ecdd1fd3a5095ef4a5/flash/colorbrewer.js
        scales = ["Blues", "Greens", "Greys", "Oranges", "Purples", "Reds", "BuGn", "BuPu", "GnBu", "OrRd", "PuBu",
                  "PuBuGn", "PuRd", "RdPu", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"]

        for header, scale in zip(headers, cycle(scales)):
            headers[header]["scale"] = scale

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
        if sum(map(len, phasing_data.values())) == 0:
            return 0, 0

        # Write parsed report data to a file
        self.write_data_file(phasing_data, "hapcut2_phasing_stats")

        pconfig = {
            'id': 'hapcut2_phasing_stats_table',
            'title': "HapCUT2 phasing stats",
            'scale': False,
            'share_key': False
        }

        # Scale headers automatically
        for metric in ["AN50", "N50", "auN", "NG50", "auNG", "QAN50", "QN50", "auQN", "QNG50", "auQNG"]:
            metric_max = max(v.get(metric, 0) for v in phasing_data.values())
            multiplier = 0.000001 if metric_max > 1_000_000 else 0.001
            suffix = " Mbp" if metric_max > 1_000_000 else " kbp"
            # Include local variable `multiplier` in lambda function.
            headers[metric]["modify"] = lambda x, multiplier=multiplier: x * multiplier
            headers[metric]["suffix"] = suffix

        # TODO - maybe split table into two tables
        #  one for comparative stats (e.g. switch errors)
        #  one for stats that do not require a reference haplotype (e.g N50)
        table_html = table.plot(phasing_data, headers.copy(), pconfig)

        # Add a report section with table
        self.add_section(
            name="Phasing stats",
            description="Table of multiple metrics relevant for phased variants.",
            helptext='''
            # Description of statistics

            **Switch rate**

            A *switch error* is defined as a position where the phase is switched from that of the previous
            heterozygous variant, when compared to the reference haplotype. Two switch errors in a row are instead
            counted as a *mismatch error*, a single position where the phase differs from the reference haplotype.
            The rate here refers to the fraction of switch errors and the the possible positions for switch error
            (the first and last variant in each phaseblock is excluded as wells as phaseblocks with less than 4
            variants). Sometimes referred to as *long switch error rate*. Note that:
            `switch rate` = `switch count` / `switch positions`.

            **Mismatch rate**

            The fraction of mismatch errors (see Switch rate for distiction to switch errors) to all possible
            positions where mismatch errors can occur (any position in a phaseblock longer than one variant is
            considered valid). Sometimes referred to as *short switch error rate*.

            **Flat rate**

            A *flat error* corresponds to the minimum hamming distance  between the two assembled haplotypes (for a
            given block) and the reference haplotype. This is an alternative metric to observing switch/mismatch
            errors in tandem. In general, this metric is thought to penalize switch errors too harshly. It may be
            of interest for a dataset with extremely low incidence of switch errors. Sometimes referred to as
            *Hamming error rate*.

            **AN50**

            AN50 metric for haplotype contiguity. Defined as the span (in base pairs) of a block such that half
            (50%) of all phased variants are in a block of that span or longer. Blocks are adjusted for unphased
            variants by multipling the base-pair span by the fraction of variants spanned by the block that are
            phased. For more info see:
            https://doi.org/10.1101%2Fgr.213462.116 and https://doi.org/10.1186%2F1471-2105-12-S1-S24

            **N50**

            N50 metric for haplotype contiguity. Defined as the span (in base pairs) of a block such that half
            (50%) of the total block length are in a block of that span or longer.

            **auN**

            auN metric for haplotype contiguity. Defined as the area under the Nx-curve where for x in [0, 100]
            A more stable metric for contiguity then the N50. See:
            https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity

            **NG50**

            NG50 metric for haplotype contiguity. Similar to N50 but relative the genome length.

            **auNG**

            auNG metric for haplotype contiguity. Similar to auN but for NGx curve.

            **Switch count**

            Switch error counts for all overlapping blocks between the assembly and reference haplotype.
            A *switch* is defined as a stretch of multiple variants being assigned to the wrong haplotype
            in relation to the reference.

            **Switch positions**

            The number of positions where switch errors are assayed. The number is the total number of variant
            pairs in blocks excluding the ends (ends reported as mismatch errors). Only blocks with at least 4
            variants recovered in both the assembly and reference are eligable for switch errors.

            **Mismatch count**

            Mismatch error counts for all overlapping blocks between the assembly and reference haplotype.
            A *mismatch* is defined as two consequtive switches over a single position in relation to the
            reference.

            **Mismatch positions**

            The number of positions where mismatch errors are assayed. This is every position shared between
            the assebled and referece haplotype in blocks of at least length 2.

            **Flat count**

            Flat error counts for all overlapping blocks between the assembly and reference haplotype. This the
            total hamming distance between the assembled and reference haplotype summed over all overlapping
            blocks.

            **Flat positions**

            The number of positions where flat errors are assayed. This is every position shared between
            the assebled and referece haplotype in blocks of at least length 2. Same as `mismatch positions`.

            **QAN50**

            Similar to AN50 but each phase block has been split at switch and mismatch locations.

            **QN50**

            Similar to N50 but each phase block has been split at switch and mismatch locations.

            **auQN**

            Similar to auN but each phase block has been split at switch and mismatch locations.

            **QNG50**

            Similar to NG50 but each phase block has been split at switch and mismatch locations.

            **auQNG**

            Similar to auNG but each phase block has been split at switch and mismatch locations.
            ''',
            plot=table_html
        )

        # Add N50 to general stats table
        general_stats_data = {
            sample: {"N50_phaseblock": data["N50"]} for sample, data in phasing_data.items() if "N50" in data
        }
        max_n50 = max(v.get("N50_phaseblock", 0) for v in general_stats_data.values())
        multiplier = 0.000001 if max_n50 > 1_000_000 else 0.001
        suffix = " Mbp" if max_n50 > 1_000_000 else " kbp"
        general_stats_header = OrderedDict({
            "N50_phaseblock": {
                'title': 'N50 block',
                'description': 'N50 statistic for phaseblock lengths',
                'scale': 'Blues',
                'modify': lambda x: x * multiplier,
                'format': '{:,.3f}',
                'suffix': suffix,
            }})

        self.general_stats_addcols(general_stats_data, general_stats_header)

        # Return if not any per-chrom statistics
        nr_stats_per_chrom = sum(data != {} for sample, data in phasing_data_per_chrom[0].items())
        has_multiple_chroms = any(len(data) > 0 for sample, data in phasing_data_per_chrom[0].items())
        if nr_stats_per_chrom == 0 or has_multiple_chroms:
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
            yield chrom, parameter, value
