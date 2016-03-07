#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_FBA
----------------------------------

Tests for `FBA` function.
"""

import unittest

import os
import tempfile

from pyabolism.io import load_model, save_model


class TestIO(unittest.TestCase):

    def setUp(self):

        self.model = load_model('examples/data/ecoli_core.xml')

    def test_sbml_io(self, buffer=True):

        model = load_model('examples/data/ecoli_core.xml', filetype='sbml')

        save_model(model, os.path.sep.join([tempfile.gettempdir(), 'test.xml']), filetype='sbml')

    def test_pickle_io(self, buffer=True):

        filepath = os.path.sep.join([tempfile.gettempdir(), 'test.xml'])

        save_model(self.model, filepath, filetype='pickle')

        load_model(filepath, filetype='pickle')

    def tearDown(self):
        pass

if __name__ == '__main__':

    unittest.main()
