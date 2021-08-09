"""
Parse configs from dict to YAML for MultiQC report
"""
from string import ascii_lowercase
from random import sample
import os


def parse_nested(d, root=None):
    for key, value in d.items():
        new_root = f"{root}.{key}" if root else str(key)
        if isinstance(value, dict):
            yield from parse_nested(value, root=new_root)
        else:
            yield new_root, value


with open(snakemake.output.yaml, "w") as out:
    unique_id = "".join(sample(ascii_lowercase, 10))
    workdir = os.path.dirname(os.path.abspath(snakemake.input.yaml))
    text = "parent_id: 'configs'\n"
    text += "parent_name: 'Configurations'\n"
    text += f"id: 'run-configs-{unique_id}'\n"
    text += f"section_name: 'Run configs'\n"
    text += "plot_type: 'html'\n"
    text += f"description: 'Configuration parameters used for run found in `{workdir}`.'\n"
    text += "data: |\n"
    text += f"\t<button class='btn btn-default' data-toggle='collapse' data-target='#{unique_id}'>Show configs</button>".expandtabs(4)
    text += f"\t<div id='{unique_id}' class='collapse'>".expandtabs(4)
    text += "\t<dl class=\"dl-horizontal\">\n".expandtabs(4)
    for param, value in parse_nested(snakemake.params.configs):
        if param.startswith("_"):
            continue
        text += f"\t\t<dt>{param}</dt><dd><samp>{value}</samp></dd>\n".expandtabs(4)
    text += "\t</dl>".expandtabs(4)
    text += "\t</div>".expandtabs(4)
    print(text, file=out)
