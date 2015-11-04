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

from pyabolism.simulate.LP import grb, GRB, generate_basic_lp


class TestFBA(unittest.TestCase):

    def setUp(self):
        
        self.model = load_model('examples/data/ecoli_core.xml')
       
    def test_FBA(self, buffer=True):
        
        FBA(self.model, show=False)
        assert (np.round(self.model.total_objective, 8) == 0.86140741)

        FBA(self.model, norm='L1', show=False)
        assert (np.round(self.model.total_objective, 8) == 0.86140740)
        
        FBA(self.model, norm='L2', show=False)
        assert (np.round(self.model.total_objective, 8) == 0.86140655)

    def test_FBA_compound_LP(self, buffer=True):
        
        conditions = [('ModelA_', -5.0),
                      ('ModelB_', -10.0),
                      ('ModelC_', -15.0),
                      ]

        results = [0.4086461041, 0.8614074126, 1.3141687211]

        lp = grb.Model('Multi_LP_testing')

        models = []
        for (model_name, lower_bound) in conditions:
            
            model = load_model('examples/data/ecoli_core.xml')
            
            model.name                                = model_name
            model.reaction['R_EX_glc_e_'].lower_bound = lower_bound
            model.lp                                  = lp
            
            generate_basic_lp(model, add_to_existing=True)
            models.append(model)

        lp.optimize()

        for model, result in zip(models, results):
            for reaction in [r for r in model.reactions() if r.objective_coefficient != 0]:
                assert np.round(reaction.lp_var.X, 8) == np.round(result, 8)

    def tearDown(self):
        pass

if __name__ == '__main__':
    
    unittest.main()
