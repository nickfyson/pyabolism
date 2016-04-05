#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_FBA
----------------------------------

Tests for `FBA` function.
"""

import unittest

from pyabolism.io import load_model
from pyabolism.simulate import FVA


class TestFBA(unittest.TestCase):

    def setUp(self):

        self.model = load_model('examples/data/ecoli_core.xml')

    def test_FVA(self, buffer=False):
        # the test of FVA is limited, since we only check whether it reproduces the behaviour seen
        # when first implemented.
        # TODO ideally we would use date from an alternative implementation in future

        FVA(self.model)

        assert round(self.model.total_objective, 4) == 0.8614

        assert round(self.model.reaction['R_ACKr'].flux_range[0], 4)           == 0.0000
        assert round(self.model.reaction['R_ACKr'].flux_range[1], 4)           == 0.0000
        assert round(self.model.reaction['R_ACONT'].flux_range[0], 4)          == 6.2649
        assert round(self.model.reaction['R_ACONT'].flux_range[1], 4)          == 6.2649
        assert round(self.model.reaction['R_ACt2r'].flux_range[0], 4)          == -0.0000
        assert round(self.model.reaction['R_ACt2r'].flux_range[1], 4)          == 0.0000
        assert round(self.model.reaction['R_ADHEr'].flux_range[0], 4)          == 0.0000
        assert round(self.model.reaction['R_ADHEr'].flux_range[1], 4)          == 0.0000
        assert round(self.model.reaction['R_ADK1'].flux_range[0], 4)           == 0.0000
        assert round(self.model.reaction['R_ADK1'].flux_range[1], 4)           == 0.0000
        assert round(self.model.reaction['R_AKGDH'].flux_range[0], 4)          == 5.3355
        assert round(self.model.reaction['R_AKGDH'].flux_range[1], 4)          == 5.3355
        assert round(self.model.reaction['R_AKGt2r'].flux_range[0], 4)         == -0.0000
        assert round(self.model.reaction['R_AKGt2r'].flux_range[1], 4)         == 0.0000
        assert round(self.model.reaction['R_ATPM'].flux_range[0], 4)           == 7.6000
        assert round(self.model.reaction['R_ATPM'].flux_range[1], 4)           == 7.6000
        assert round(self.model.reaction['R_ATPS4r'].flux_range[0], 4)         == 39.7470
        assert round(self.model.reaction['R_ATPS4r'].flux_range[1], 4)         == 39.7470
        assert round(self.model.reaction['R_CO2t'].flux_range[0], 4)           == -23.3424
        assert round(self.model.reaction['R_CO2t'].flux_range[1], 4)           == -23.3424
        assert round(self.model.reaction['R_CS'].flux_range[0], 4)             == 6.2649
        assert round(self.model.reaction['R_CS'].flux_range[1], 4)             == 6.2649
        assert round(self.model.reaction['R_CYTBD'].flux_range[0], 4)          == 44.6930
        assert round(self.model.reaction['R_CYTBD'].flux_range[1], 4)          == 44.6930
        assert round(self.model.reaction['R_D_LACt2'].flux_range[0], 4)        == -0.0000
        assert round(self.model.reaction['R_D_LACt2'].flux_range[1], 4)        == 0.0000
        assert round(self.model.reaction['R_ENO'].flux_range[0], 4)            == 14.8491
        assert round(self.model.reaction['R_ENO'].flux_range[1], 4)            == 14.8491
        assert round(self.model.reaction['R_ETOHt2r'].flux_range[0], 4)        == -0.0000
        assert round(self.model.reaction['R_ETOHt2r'].flux_range[1], 4)        == 0.0000
        assert round(self.model.reaction['R_EX_ac_e_'].flux_range[0], 4)       == 0.0000
        assert round(self.model.reaction['R_EX_ac_e_'].flux_range[1], 4)       == 0.0000
        assert round(self.model.reaction['R_EX_akg_e_'].flux_range[0], 4)      == 0.0000
        assert round(self.model.reaction['R_EX_akg_e_'].flux_range[1], 4)      == 0.0000
        assert round(self.model.reaction['R_EX_co2_e_'].flux_range[0], 4)      == 23.3424
        assert round(self.model.reaction['R_EX_co2_e_'].flux_range[1], 4)      == 23.3424
        assert round(self.model.reaction['R_EX_etoh_e_'].flux_range[0], 4)     == 0.0000
        assert round(self.model.reaction['R_EX_etoh_e_'].flux_range[1], 4)     == 0.0000
        assert round(self.model.reaction['R_EX_for_e_'].flux_range[0], 4)      == 0.0000
        assert round(self.model.reaction['R_EX_for_e_'].flux_range[1], 4)      == -0.0000
        assert round(self.model.reaction['R_EX_fum_e_'].flux_range[0], 4)      == 0.0000
        assert round(self.model.reaction['R_EX_fum_e_'].flux_range[1], 4)      == 0.0000
        assert round(self.model.reaction['R_EX_glc_e_'].flux_range[0], 4)      == -10.0000
        assert round(self.model.reaction['R_EX_glc_e_'].flux_range[1], 4)      == -10.0000
        assert round(self.model.reaction['R_EX_h2o_e_'].flux_range[0], 4)      == 24.9201
        assert round(self.model.reaction['R_EX_h2o_e_'].flux_range[1], 4)      == 24.9201
        assert round(self.model.reaction['R_EX_h_e_'].flux_range[0], 4)        == 9.1129
        assert round(self.model.reaction['R_EX_h_e_'].flux_range[1], 4)        == 9.1129
        assert round(self.model.reaction['R_EX_lac_D_e_'].flux_range[0], 4)    == 0.0000
        assert round(self.model.reaction['R_EX_lac_D_e_'].flux_range[1], 4)    == 0.0000
        assert round(self.model.reaction['R_EX_o2_e_'].flux_range[0], 4)       == -22.3465
        assert round(self.model.reaction['R_EX_o2_e_'].flux_range[1], 4)       == -22.3465
        assert round(self.model.reaction['R_EX_pi_e_'].flux_range[0], 4)       == -3.1689
        assert round(self.model.reaction['R_EX_pi_e_'].flux_range[1], 4)       == -3.1689
        assert round(self.model.reaction['R_EX_pyr_e_'].flux_range[0], 4)      == 0.0000
        assert round(self.model.reaction['R_EX_pyr_e_'].flux_range[1], 4)      == -0.0000
        assert round(self.model.reaction['R_EX_succ_e_'].flux_range[0], 4)     == 0.0000
        assert round(self.model.reaction['R_EX_succ_e_'].flux_range[1], 4)     == 0.0000
        assert round(self.model.reaction['R_FBA'].flux_range[0], 4)            == 7.5708
        assert round(self.model.reaction['R_FBA'].flux_range[1], 4)            == 7.5708
        assert round(self.model.reaction['R_FBP'].flux_range[0], 4)            == 0.0000
        assert round(self.model.reaction['R_FBP'].flux_range[1], 4)            == 0.0000
        assert round(self.model.reaction['R_FORt'].flux_range[0], 4)           == 0.0000
        assert round(self.model.reaction['R_FORt'].flux_range[1], 4)           == 0.0000
        assert round(self.model.reaction['R_FRD'].flux_range[0], 4)            == 0.0000
        assert round(self.model.reaction['R_FRD'].flux_range[1], 4)            == 999993.6645
        assert round(self.model.reaction['R_FUM'].flux_range[0], 4)            == 5.3355
        assert round(self.model.reaction['R_FUM'].flux_range[1], 4)            == 5.3355
        assert round(self.model.reaction['R_FUMt2_2'].flux_range[0], 4)        == 0.0000
        assert round(self.model.reaction['R_FUMt2_2'].flux_range[1], 4)        == 0.0000
        assert round(self.model.reaction['R_G6PDH2r'].flux_range[0], 4)        == 4.7171
        assert round(self.model.reaction['R_G6PDH2r'].flux_range[1], 4)        == 4.7171
        assert round(self.model.reaction['R_GAPD'].flux_range[0], 4)           == 16.1377
        assert round(self.model.reaction['R_GAPD'].flux_range[1], 4)           == 16.1377
        assert round(self.model.reaction['R_GLCpts'].flux_range[0], 4)         == 10.0000
        assert round(self.model.reaction['R_GLCpts'].flux_range[1], 4)         == 10.0000
        assert round(self.model.reaction['R_GND'].flux_range[0], 4)            == 4.7171
        assert round(self.model.reaction['R_GND'].flux_range[1], 4)            == 4.7171
        assert round(self.model.reaction['R_H2Ot'].flux_range[0], 4)           == -24.9201
        assert round(self.model.reaction['R_H2Ot'].flux_range[1], 4)           == -24.9201
        assert round(self.model.reaction['R_ICDHyr'].flux_range[0], 4)         == 6.2649
        assert round(self.model.reaction['R_ICDHyr'].flux_range[1], 4)         == 6.2649
        assert round(self.model.reaction['R_ICL'].flux_range[0], 4)            == 0.0000
        assert round(self.model.reaction['R_ICL'].flux_range[1], 4)            == 0.0000
        assert round(self.model.reaction['R_LDH_D'].flux_range[0], 4)          == -0.0000
        assert round(self.model.reaction['R_LDH_D'].flux_range[1], 4)          == 0.0000
        assert round(self.model.reaction['R_MALS'].flux_range[0], 4)           == 0.0000
        assert round(self.model.reaction['R_MALS'].flux_range[1], 4)           == 0.0000
        assert round(self.model.reaction['R_MDH'].flux_range[0], 4)            == 5.3355
        assert round(self.model.reaction['R_MDH'].flux_range[1], 4)            == 5.3355
        assert round(self.model.reaction['R_ME1'].flux_range[0], 4)            == 0.0000
        assert round(self.model.reaction['R_ME1'].flux_range[1], 4)            == 0.0000
        assert round(self.model.reaction['R_ME2'].flux_range[0], 4)            == 0.0000
        assert round(self.model.reaction['R_ME2'].flux_range[1], 4)            == 0.0000
        assert round(self.model.reaction['R_NADH11'].flux_range[0], 4)         == 39.3575
        assert round(self.model.reaction['R_NADH11'].flux_range[1], 4)         == 39.3575
        assert round(self.model.reaction['R_NADTRHD'].flux_range[0], 4)        == 0.0000
        assert round(self.model.reaction['R_NADTRHD'].flux_range[1], 4)        == 0.0000
        assert round(self.model.reaction['R_O2t'].flux_range[0], 4)            == 22.3465
        assert round(self.model.reaction['R_O2t'].flux_range[1], 4)            == 22.3465
        assert round(self.model.reaction['R_PDH'].flux_range[0], 4)            == 9.4933
        assert round(self.model.reaction['R_PDH'].flux_range[1], 4)            == 9.4933
        assert round(self.model.reaction['R_PFK'].flux_range[0], 4)            == 7.5708
        assert round(self.model.reaction['R_PFK'].flux_range[1], 4)            == 7.5708
        assert round(self.model.reaction['R_PFL'].flux_range[0], 4)            == 0.0000
        assert round(self.model.reaction['R_PFL'].flux_range[1], 4)            == -0.0000
        assert round(self.model.reaction['R_PGI'].flux_range[0], 4)            == 5.1063
        assert round(self.model.reaction['R_PGI'].flux_range[1], 4)            == 5.1063
        assert round(self.model.reaction['R_PGK'].flux_range[0], 4)            == -16.1377
        assert round(self.model.reaction['R_PGK'].flux_range[1], 4)            == -16.1377
        assert round(self.model.reaction['R_PGL'].flux_range[0], 4)            == 4.7171
        assert round(self.model.reaction['R_PGL'].flux_range[1], 4)            == 4.7171
        assert round(self.model.reaction['R_PGM'].flux_range[0], 4)            == -14.8491
        assert round(self.model.reaction['R_PGM'].flux_range[1], 4)            == -14.8491
        assert round(self.model.reaction['R_PIt'].flux_range[0], 4)            == -3.1689
        assert round(self.model.reaction['R_PIt'].flux_range[1], 4)            == -3.1689
        assert round(self.model.reaction['R_PPC'].flux_range[0], 4)            == 2.4684
        assert round(self.model.reaction['R_PPC'].flux_range[1], 4)            == 2.4684
        assert round(self.model.reaction['R_PPCK'].flux_range[0], 4)           == 0.0000
        assert round(self.model.reaction['R_PPCK'].flux_range[1], 4)           == 0.0000
        assert round(self.model.reaction['R_PPS'].flux_range[0], 4)            == 0.0000
        assert round(self.model.reaction['R_PPS'].flux_range[1], 4)            == 0.0000
        assert round(self.model.reaction['R_PTAr'].flux_range[0], 4)           == 0.0000
        assert round(self.model.reaction['R_PTAr'].flux_range[1], 4)           == 0.0000
        assert round(self.model.reaction['R_PYK'].flux_range[0], 4)            == 1.9335
        assert round(self.model.reaction['R_PYK'].flux_range[1], 4)            == 1.9335
        assert round(self.model.reaction['R_PYRt2r'].flux_range[0], 4)         == 0.0000
        assert round(self.model.reaction['R_PYRt2r'].flux_range[1], 4)         == 0.0000
        assert round(self.model.reaction['R_RPE'].flux_range[0], 4)            == 2.5256
        assert round(self.model.reaction['R_RPE'].flux_range[1], 4)            == 2.5256
        assert round(self.model.reaction['R_RPI'].flux_range[0], 4)            == -2.1916
        assert round(self.model.reaction['R_RPI'].flux_range[1], 4)            == -2.1916
        assert round(self.model.reaction['R_SUCCt2_2'].flux_range[0], 4)       == 0.0000
        assert round(self.model.reaction['R_SUCCt2_2'].flux_range[1], 4)       == 0.0000
        assert round(self.model.reaction['R_SUCCt2b'].flux_range[0], 4)        == 0.0000
        assert round(self.model.reaction['R_SUCCt2b'].flux_range[1], 4)        == 0.0000
        assert round(self.model.reaction['R_SUCD1i'].flux_range[0], 4)         == 5.3355
        assert round(self.model.reaction['R_SUCD1i'].flux_range[1], 4)         == 999999.0000
        assert round(self.model.reaction['R_SUCD4'].flux_range[0], 4)          == 5.3355
        assert round(self.model.reaction['R_SUCD4'].flux_range[1], 4)          == 5.3355
        assert round(self.model.reaction['R_SUCOAS'].flux_range[0], 4)         == -5.3355
        assert round(self.model.reaction['R_SUCOAS'].flux_range[1], 4)         == -5.3355
        assert round(self.model.reaction['R_TALA'].flux_range[0], 4)           == 1.4183
        assert round(self.model.reaction['R_TALA'].flux_range[1], 4)           == 1.4183
        assert round(self.model.reaction['R_THD2'].flux_range[0], 4)           == 0.0000
        assert round(self.model.reaction['R_THD2'].flux_range[1], 4)           == 0.0000
        assert round(self.model.reaction['R_TKT1'].flux_range[0], 4)           == 1.4183
        assert round(self.model.reaction['R_TKT1'].flux_range[1], 4)           == 1.4183
        assert round(self.model.reaction['R_TKT2'].flux_range[0], 4)           == 1.1073
        assert round(self.model.reaction['R_TKT2'].flux_range[1], 4)           == 1.1073
        assert round(self.model.reaction['R_TPI'].flux_range[0], 4)            == 7.5708
        assert round(self.model.reaction['R_TPI'].flux_range[1], 4)            == 7.5708

        # and this is the somewhat hacky way to generate the code for the above checks!
        # to be run when confident the code is working. not ideal, but all we have...
        # with open('tempfile.txt', 'w') as f:
        #     for r in self.model.reactions():
        #         f.write('assert ')
        #         f.write(("round(self.model.reaction['%s'].flux_range[0], 4)" % (r.id)).ljust(63))
        #         f.write(" == %.4f" % (r.flux_range[0]))
        #         f.write('\n')
        #         f.write('assert ')
        #         f.write(("round(self.model.reaction['%s'].flux_range[1], 4)" % (r.id)).ljust(63))
        #         f.write(" == %.4f" % (r.flux_range[1]))
        #         f.write('\n')

    def tearDown(self):
        pass

if __name__ == '__main__':

    unittest.main()
