import os
import sys

from setuptools import setup

# ensure the current directory is on sys.path so versioneer can be imported
# when pip uses PEP 517/518 build rules.
# https://github.com/python-versioneer/python-versioneer/issues/193
sys.path.append(os.path.dirname(__file__))

# import versioneer  # noqa: E402


setup(
    name="Iliad",
    # version=versioneer.get_version(),
    # cmdclass=versioneer.get_cmdclass(),
    author="Noah Herrick",
    author_email="ncherric@iupui.edu",
    description="Snakemake workflow for Genomic data processing "
    "from genotyping by sequence or genotyping by microarray.",
    long_description_content_type="text/markdown",
    zip_safe=False,
    license="MIT",
    url="https://iliad.readthedocs.io",
    project_urls={
        "Source": "https://github.com/iliad",
    },
    entry_points={
    },
    package_data={"": ["*.css", "*.sh", "*.html", "*.jinja2", "*.js", "*.svg"]},
    python_requires=">=3.8",
    install_requires=[
        "pyyaml",
        "configargparse",
        "appdirs",
        "datrie",
        "jsonschema",
        "docutils",
        "gitpython",
        "psutil",
    ],
)