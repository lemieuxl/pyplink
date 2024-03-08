#!/usr/bin/env python
"""Setup."""

# How to build source distribution
#   - python setup.py sdist --format bztar
#   - python setup.py sdist --format gztar
#   - python setup.py sdist --format zip
#   - python setup.py bdist_wheel


import os
from setuptools import setup


MAJOR = 1
MINOR = 3
MICRO = 8
VERSION = f"{MAJOR}.{MINOR}.{MICRO}b1"


def write_version_file(fn=None):
    """Write the version file."""
    if fn is None:
        fn = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            os.path.join("pyplink", "version.py"),
        )

    content = (
        "\n# THIS FILE WAS GENERATED AUTOMATICALLY BY PYPLINK SETUP.PY\n"
        'pyplink_version = "{version}"\n'
    )

    with open(fn, "w") as f:
        f.write(content.format(version=VERSION))


def setup_package():
    """Setup the package."""
    # Saving the version into a file
    write_version_file()

    setup(
        name="pyplink",
        version=VERSION,
        description="Python module to read binary Plink files.",
        author="Louis-Philippe Lemieux Perreault",
        author_email="louis-philippe.lemieux.perreault@statgen.org",
        url="https://github.com/lemieuxl/pyplink",
        license="MIT",
        packages=["pyplink", "pyplink.tests"],
        package_data={"pyplink.tests": ["data/*"], },
        test_suite="pyplink.tests.test_suite",
        install_requires=[
            "numpy >= 1.8.2",
            "pandas >= 0.17.1",
            "importlib_resources >= 5.12.0",
        ],
        classifiers=["Operating System :: POSIX :: Linux",
                     "Operating System :: MacOS :: MacOS X",
                     "Operating System :: Microsoft",
                     "Programming Language :: Python",
                     "Programming Language :: Python :: 3",
                     "Programming Language :: Python :: 3.7",
                     "Programming Language :: Python :: 3.8",
                     "Programming Language :: Python :: 3.9",
                     "Programming Language :: Python :: 3.10",
                     "Programming Language :: Python :: 3.11",
                     "Programming Language :: Python :: 3.12",
                     "License :: OSI Approved :: MIT License",
                     "Topic :: Scientific/Engineering :: Bio-Informatics"],
        keywords="bioinformatics format Plink binary",
    )


if __name__ == "__main__":
    setup_package()
