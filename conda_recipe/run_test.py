

import sys
import unittest

from pyplink.tests import test_suite as pyplink_test_suite


result = unittest.TextTestRunner().run(pyplink_test_suite)
sys.exit(0 if result.wasSuccessful() else 1)
