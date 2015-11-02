#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_FBA
----------------------------------

Tests for `FBA` function.
"""

import unittest

import numpy as np

from pyabolism.io import load_model
from pyabolism.simulate import FBA


class TestFBA(unittest.TestCase):

    def setUp(self):
        
        self.model = load_model('examples/data/ecoli_core.xml')
       
    def test_FBA(self, buffer=False):
        
        FBA(self.model, show=False)
        assert (np.round(self.model.total_objective, 8) == 0.86140741)

        FBA(self.model, norm='L1', show=False)
        assert (np.round(self.model.total_objective, 8) == 0.86140740)
        
        FBA(self.model, norm='L2', show=False)
        assert (np.round(self.model.total_objective, 8) == 0.86140655)

    def tearDown(self):
        pass

if __name__ == '__main__':
    
    unittest.main()
