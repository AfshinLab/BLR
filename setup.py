from setuptools import setup, find_namespace_packages

with open("README.md") as f:
    long_description = f.read()

setup(
    name="blr",
    author="Tobias Frick",
    url="https://github.com/TobiasFrick/BLR/",
    description="Barcoded long reads pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.6",
    install_requires=[
        "pysam",
        "dnaio",
        "tqdm",
        "pandas",
        "numpy",
        "matplotlib",
        "snakemake",
        "importlib_resources; python_version<'3.7'",
        "dataclasses; python_version<'3.7'",
        "multiqc"
    ],
    extras_require={
        "dev": ["flake8"],
    },
    package_dir={"": "src"},
    packages=find_namespace_packages("src"),
    package_data={
        "blr": [
            "Snakefile",
            "rules/*.smk",
            "config.schema.yaml",
            "blr.yaml",
            "multiqc_config.yaml",
            "naibr.config",
            "naibr-environment.yml",
            "run_anew.smk",
        ]
    },
    entry_points={
        "console_scripts": [
            "blr = blr.__main__:main"
        ],
        "multiqc.modules.v1": [
            "stats = multiqc_blr.modules.stats:MultiqcModule",
            "hapcut2 = multiqc_blr.modules.hapcut2:MultiqcModule",
            "whatshap = multiqc_blr.modules.whatshap:MultiqcModule",
        ],
        "multiqc.cli_options.v1": [
            "disable_plugin = multiqc_blr.cli:disable_plugin"
        ],
        "multiqc.hooks.v1": [
            "before_config = multiqc_blr.multiqc_blr:before_config",
            "execution_start = multiqc_blr.multiqc_blr:execution_start"
        ],
        "multiqc.templates.v1": [
            "blr = multiqc_blr.templates.blr"
        ]
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ]
)
