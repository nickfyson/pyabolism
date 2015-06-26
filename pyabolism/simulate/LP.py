
import os

import gurobipy as grb
from gurobipy import GRB

grb.setParam('OutputFlag', 0)
grb.setParam('LogFile', os.getenv("HOME") + '/gurobi.log')
try:
    os.unlink('gurobi.log')
except:
    pass


def generate_basic_lp(model):
    """docstring for _generate_basic_lp"""
    # we first initialise the linear program
        
    lp = grb.Model('fba')
    
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
                      reaction.id)
        else:
            lp.addVar(max(reaction.lower_bound, 0.0),
                      reaction.upper_bound,
                      reaction.objective_coefficient,
                      GRB.CONTINUOUS,
                      reaction.id)
    lp.update()
    
    # each metabolite becomes a constraint,
    # since subject to stoichiometry of reactions the net production must be zero
    for metabolite in model.metabolites():
        
        # each reaction that features the metabolite becomes an element in the constraint
        variables       = []
        stoichiometries = []
        for reaction in model.reaction.get_by_contains(metabolite):
            variables.append(lp.getVarByName(reaction.id))
            stoichiometries.append(reaction.participants[metabolite])

        if metabolite.boundaryCondition:
            # by the nomenclature of SBML, being a boundary metabolite means
            # the net production is *not* constrained
            pass
        else:
            # the weighted sum of all the reactions featuring the metabolite
            # must have a net total of zero
            lp.addConstr(grb.LinExpr(stoichiometries, variables), GRB.EQUAL, 0.0, metabolite.id)
    
    lp.update()
    
    return lp
