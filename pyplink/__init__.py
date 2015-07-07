"""A module to read Plink's binary files."""

# This file is part of pyplink.
#
# The MIT License (MIT)
#
# Copyright (c) 2014 Louis-Philippe Lemieux Perreault
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


from .pyplink import *

try:
    from .version import pyplink_version as __version__
except ImportError:
    __version__ = None


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014 Louis-Philippe Lemieux Perreault"
__credits__ = ["Louis-Philippe Lemieux Perreault"]
__license__ = "MIT"
__maintainer__ = "Louis-Philippe Lemieux Perreault"
__email__ = "louis-philippe.lemieux.perreault@statgen.org"
__status__ = "Development"


def test(verbosity=1):
    """Executes all the tests for pyplink.

    Args:
        verbosity (int): the verbosity level for unittest

    Just set ``verbosity`` to an integer higher than ``1`` to have more
    information about the tests.

    """
    import unittest
    from .tests import test_suite

    # Testing
    unittest.TextTestRunner(verbosity=verbosity).run(test_suite)
