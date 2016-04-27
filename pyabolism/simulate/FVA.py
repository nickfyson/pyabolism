
import numpy as np

from .LP import GRB
from FBA import FBA


def FVA(model, norm='L2', obj_ratio=1.0, show=False):

    assert obj_ratio > 0.0, 'obj_ratio must be strictly positive.'

    # we run standard FBA to find the maximum objective attainable
    FBA(model, show=False, norm=norm)

    # we extract the linear program built for the original FBA problem
    lp = model.lp

    # we add constraints to the linear program such that the objective value is maintained
    # up to the proportion indicated by obj_ratio
    for r in [r for r in model.reactions() if r.objective_coefficient != 0]:
        if r.flux_value > 0:
            r.lp_var.lb = float(obj_ratio) * r.flux_value
            r.lp_var.ub = r.flux_value * np.infty
        else:
            # in the case the we have (for whatever reason) a negative objective flux
            # the bounds must be set such that the objective remains just as *negative*
            # as the original solution
            r.lp_var.ub = float(obj_ratio) * r.flux_value
            r.lp_var.lb = r.flux_value * np.infty

    # we iterate over all the reactions in the model to determine the range of freedom available
    # in each flux value
    # we reuse the existing LP, saving the overheard of building from scratch and allowing
    # the optimization library a 'hot start' on solving each problem
    for r in model.reactions():

        # we can extract the variable tied to each reaction
        variable = r.lp_var
        # and set it as the sole objective for the LP
        lp.setObjective(variable)

        # we first extract the minimum permitted value for the flux...
        lp.modelSense = GRB.MINIMIZE
        lp.optimize()
        minimum = variable.X

        # ...and then likewise obtain the maximum
        lp.modelSense = GRB.MAXIMIZE
        lp.optimize()
        maximum = variable.X

        # the solution to this problem is no longer a single flux value for each reaction
        # but instead we store a tuple with the min and max values
        r.flux_range = (minimum, maximum)

    return

