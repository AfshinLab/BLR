"""
Parse configs from dict to YAML for MultiQC report
"""

def parse_nested(d, root=None):
    for key, value in d.items():
        new_root = f"{root}.{key}" if root else str(key)
        if isinstance(value, dict):
            yield from parse_nested(value, root=new_root)
        else:
            yield new_root, value

with open(snakemake.output.yaml, "w") as out:
    text = "id: 'run_configs'\n"
    text += "section_name: 'Run configs'\n"
    text += "plot_type: 'html'\n"
    text += "description: ' for current analysis. Describes which parameter values was used.'\n"
    text += "data: |\n"
    text += "\t<dl class=\"dl-horizontal\">\n".expandtabs(4)
    for param, value in parse_nested(snakemake.params.configs):
        if param.startswith("_"):
            continue
        text += f"\t\t<dt>{param}</dt><dd><samp>{value}</samp></dd>\n".expandtabs(4)
    text += "\t</dl>".expandtabs(4)
    print(text, file=out)
