#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_EFlux
----------------------------------

Tests for `EFlux` function.
"""
import __future__

import unittest

import numpy as np

from pyabolism.io import load_model
from pyabolism.simulate import EFlux

class TestEFlux(unittest.TestCase):

    def setUp(self):
        
        
        self.model = load_model('examples/data/ecoli_core.xml')

        gene_ids = ['b2926', 'b2925', 'b0008', 'b3734', 'b3735', 'b3736', 'b3737', 'b3731',
                    'b3732', 'b0767', 'b3738', 'b3739', 'b1702', 'b2914', 'b3919', 'b3916',
                    'b2587', 'b3117', 'b0722', 'b0723', 'b0720', 'b0721', 'b0726', 'b0727',
                    'b0724', 'b0728', 'b1380', 'b1136', 'b1241', 'b0474', 'b1417', 'b1416',
                    'b2276', 'b2277', 'b2278', 'b2279', 'b2465', 'b2463', 'b2779', 'b1276',
                    'b4025', 'b2975', 'b2976', 'b0729', 'b3603', 'b2416', 'b2417', 'b2415',
                    'b1479', 'b2287', 'b1779']

        from random import random, seed

        self.expressions_A = {}
        for gid in gene_ids:
            seed(gid)
            self.expressions_A[gid] = 100 * random()

        self.expressions_B = {}
        for gid in gene_ids:
            seed(gid + gid + gid)
            self.expressions_B[gid] = 100 * random()

    def test_EFlux(self, buffer=False):
        
        EFlux(self.model, self.expressions_A, show=True)
        assert (np.round(self.model.total_objective, 8) == 0.44314716)
        
        EFlux(self.model, self.expressions_B, show=True)
        assert (np.round(self.model.total_objective, 8) == 0.30007236)
    
    def tearDown(self):
        pass

if __name__ == '__main__':
    
    unittest.main()
