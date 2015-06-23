#!/usr/bin/env python

# How to build source distribution
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip


import os
from setuptools import setup


MAJOR = 0
MINOR = 3
VERSION = "{}.{}".format(MAJOR, MINOR)


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


def setup_package():
    # Saving the version into a file
    write_version_file()

    setup(
        name="pyplink",
        version=VERSION,
        description="Python module to read binary Plink files.",
        author="Louis-Philippe Lemieux Perreault",
        author_email="louis-philippe.lemieux.perreault@statgen.org",
        url="http://www.statgen.org",
        license="GPL",
        packages=["pyplink"],
        install_requires=["numpy >= 1.8.2", "pandas >= 0.14.1"],
        classifiers=["Operating System :: Linux",
                     "Programming Language :: Python",
                     "Programming Language :: Python :: 2.7",
                     "Programming Language :: Python :: 3.3",
                     "Programming Language :: Python :: 3.4",
                     "License :: OSI Approved :: MIT License"],
    )

    return


if __name__ == "__main__":
    setup_package()
