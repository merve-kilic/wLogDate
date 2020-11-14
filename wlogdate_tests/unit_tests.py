import unittest
from os.path import dirname, realpath
from wlogdate_tests.utils import compute_divergence_time, test_logIt
from treeswift import *

path=dirname(realpath(__file__))


class TestLogIt(unittest.TestCase):
    """Test the logIt function"""

    def test1(self):
        """Use the published wLogDate tree as initial"""
        inTreeFile = path + "/test_data/small_test/viruses/Test1/input.nwk"
        startTreeFile = path + "/test_data/small_test/viruses/Test1/wld_ref.nwk"
        smplTimeFile = path + "/test_data/small_test/viruses/Test1/sampling_time.txt"
        refFile = path + "/test_data/small_test/viruses/Test1/true_vs_wld_ref.txt"
        failures, max_violate= test_logIt(inTreeFile, startTreeFile, smplTimeFile,refFile)
        self.assertTrue(failures > 0, msg="Failed test 1 with maximum relative error " + str(max_violate*100) + "%")

if __name__ == '__main__':
    unittest.main()
