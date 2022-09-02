"""
Parse version from dict to YAML for MultiQC report
"""
from string import ascii_lowercase
from random import sample
import os

from blr.cli.config import load_yaml


with open(snakemake.output.yaml, "w") as out:  # noqa: F821
    unique_id = "".join(sample(ascii_lowercase, 10))
    workdir = os.path.dirname(os.path.abspath(snakemake.input.yaml))  # noqa: F821
    text = "parent_id: 'about'\n"
    text += "parent_name: 'About'\n"
    text += f"id: 'versions-{unique_id}'\n"
    text += "section_name: 'Software versions'\n"
    text += "plot_type: 'html'\n"
    text += f"description: 'Software versions used for run found in `{workdir}`.'\n"
    text += "data: |\n"
    text += "\t<dl class=\"dl-horizontal\">\n".expandtabs(4)
    versions, _ = load_yaml(snakemake.input.yaml)  # noqa: F821
    for software, version in versions.items():
        text += f"\t\t<dt>{software}</dt><dd><samp>{version}</samp></dd>\n".expandtabs(4)
    text += "\t</dl>".expandtabs(4)
    print(text, file=out)
