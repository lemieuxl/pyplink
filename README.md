[![Build Status](https://github.com/lemieuxl/pyplink/actions/workflows/python-tests.yml/badge.svg?branch=master)](https://github.com/lemieuxl/pyplink/actions)
[![PyPI version](https://badge.fury.io/py/pyplink.svg)](http://badge.fury.io/py/pyplink)

# pyplink - Module to process Plink's binary files

`PyPlink` is a Python module to read and write Plink's binary files. Short
documentation available at
[https://lemieuxl.github.io/pyplink/](https://lemieuxl.github.io/pyplink/).

## Dependencies

The tool requires a standard [Python](http://python.org/) installation (3.7 or
higher are supported) with the following modules:

1. [numpy](http://www.numpy.org/)
2. [pandas](http://pandas.pydata.org/)

The tool has been tested on *Linux* only, but should work on *MacOS* and
*Windows* operating systems as well.

## Installation

Using `pip`:

```bash
pip install pyplink
```

Using `conda` (make sure [Bioconda](https://bioconda.github.io/) is configured;
thanks to @nebfield for the recipe):

```bash
conda install pyplink
```

### Updating

To update the module using `pip`:

```bash
pip install -U pyplink
```

To update the module using `conda`:

```bash
# Make sure bioconda is configured
conda update pyplink
```

## Testing

To test the module, just perform the following command:

```console
$ python -m pyplink.tests
.............................................
----------------------------------------------------------------------
Ran 45 tests in 0.334s

OK
```

## Example

The following
[notebook](http://nbviewer.ipython.org/github/lemieuxl/pyplink/blob/master/demo/PyPlink%20Demo.ipynb)
contains a demonstration of the `PyPlink` module.
