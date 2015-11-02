
import os
import random
import string

import gurobipy as grb
from gurobipy import GRB

grb.setParam('OutputFlag', 0)
grb.setParam('LogFile', os.getenv("HOME") + '/gurobi.log')
try:
    os.unlink('gurobi.log')
except:
    pass


def generate_basic_lp(model, lp=None):
    """docstring for _generate_basic_lp"""
    # we first initialise the linear program
    if lp:
        prefix = lp.getAttr('ModelName') + '_'
    else:
        lp = grb.Model(''.join(random.choice(string.ascii_uppercase) for _ in range(5)))
        prefix = lp.getAttr('ModelName') + '_'

    lp.modelSense = GRB.MAXIMIZE
    
    # reaction fluxes are the variables in the linear model
    # we take the objective coefficient from the corresponding reaction property,
    # which generally be non-zero for only the biomass reaction
    for reaction in model.reactions():
        
        if reaction.reversible:
            lp.addVar(reaction.lower_bound,
                      reaction.upper_bound,
                      reaction.objective_coefficient,
                      GRB.CONTINUOUS,
                      prefix + reaction.id)
        else:
            lp.addVar(max(reaction.lower_bound, 0.0),
                      reaction.upper_bound,
                      reaction.objective_coefficient,
                      GRB.CONTINUOUS,
                      prefix + reaction.id)
    lp.update()
    
    # each metabolite becomes a constraint,
    # since subject to stoichiometry of reactions the net production must be zero
    for metabolite in model.metabolites():
        
        # each reaction that features the metabolite becomes an element in the constraint
        variables       = []
        stoichiometries = []
        for reaction in model.reaction.get_by_contains(metabolite):
            variables.append(lp.getVarByName(prefix + reaction.id))
            stoichiometries.append(reaction.participants[metabolite])

        if metabolite.boundaryCondition:
            # by the nomenclature of SBML, being a boundary metabolite means
            # the net production is *not* constrained
            pass
        else:
            # the weighted sum of all the reactions featuring the metabolite
            # must have a net total of zero
            lp.addConstr(grb.LinExpr(stoichiometries, variables), GRB.EQUAL, 0.0, prefix + metabolite.id)
    
    lp.update()
    
    return lp
