#!/usr/bin/env python
"""Setup."""

# How to build source distribution
#   - python setup.py sdist --format bztar
#   - python setup.py sdist --format gztar
#   - python setup.py sdist --format zip
#   - python setup.py bdist_wheel --universal

# How to build for conda (do both with 2.7 and 3.4)
#   - cd conda_recipe
#   - conda clean -ytps; conda build purge; conda build --python $VERSION .
#   - cp $FILE ../conda_dist/linux-64
#   - conda convert -p all ../conda_dist/linux-64/$FILE -o ../conda_dist
#   - cd ../conda_dist && conda index *


import os
import sys
from setuptools import setup


MAJOR = 1
MINOR = 3
MICRO = 6
VERSION = "{}.{}.{}".format(MAJOR, MINOR, MICRO)


def write_version_file(fn=None):
    if fn is None:
        fn = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            os.path.join("pyplink", "version.py"),
        )

    content = (
        "\n# THIS FILE WAS GENERATED AUTOMATICALLY BY PYPLINK SETUP.PY\n"
        'pyplink_version = "{version}"\n'
    )

    a = open(fn, "w")
    try:
        a.write(content.format(version=VERSION))
    finally:
        a.close()


def get_requirements():
    # Initial requirements
    requirements = ["numpy >= 1.8.2", "pandas >= 0.17.1", "six >= 1.9.0"]

    # Checking if python 2 (requires mock)
    if sys.version_info[0] == 2:
        requirements.append(["mock >= 2.0.0"])

    return requirements


def setup_package():
    # Saving the version into a file
    write_version_file()

    setup(
        zip_safe=False,
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
        install_requires=get_requirements(),
        classifiers=["Operating System :: POSIX :: Linux",
                     "Operating System :: MacOS :: MacOS X",
                     "Operating System :: Microsoft",
                     "Programming Language :: Python",
                     "Programming Language :: Python :: 2.7",
                     "Programming Language :: Python :: 3",
                     "Programming Language :: Python :: 3.7",
                     "Programming Language :: Python :: 3.8",
                     "Programming Language :: Python :: 3.9",
                     "Programming Language :: Python :: 3.10",
                     "Programming Language :: Python :: 3.11",
                     "License :: OSI Approved :: MIT License",
                     "Topic :: Scientific/Engineering :: Bio-Informatics"],
        keywords="bioinformatics format Plink binary",
    )

    return


if __name__ == "__main__":
    setup_package()
