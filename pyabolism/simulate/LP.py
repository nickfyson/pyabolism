
import os

import gurobipy as grb
from gurobipy import GRB

grb.setParam('OutputFlag', 0)
grb.setParam('LogFile', os.getenv("HOME") + '/gurobi.log')
try:
    os.unlink('gurobi.log')
except:
    pass


def generate_basic_lp(model, add_to_existing=False):
    """docstring for generate_basic_lp"""
    # we first initialise the linear program
    if hasattr(model, 'lp') and add_to_existing:
        pass
    else:
        model.lp = grb.Model('LP_' + model.name)

    model.lp.modelSense = GRB.MAXIMIZE

    # reaction fluxes are the variables in the linear model
    # we take the objective coefficient from the corresponding reaction property,
    # which generally be non-zero for only the biomass reaction
    for reaction in model.reactions():

        if reaction.reversible:
            var = model.lp.addVar(reaction.lower_bound,
                                  reaction.upper_bound,
                                  reaction.objective_coefficient,
                                  GRB.CONTINUOUS,
                                  model.name + reaction.id)
        else:
            var = model.lp.addVar(max(reaction.lower_bound, 0.0),
                                  reaction.upper_bound,
                                  reaction.objective_coefficient,
                                  GRB.CONTINUOUS,
                                  model.name + reaction.id)
        reaction.lp_var = var

    model.lp.update()

    # each metabolite becomes a constraint,
    # since subject to stoichiometry of reactions the net production must be zero
    for metabolite in model.metabolites():

        # each reaction that features the metabolite becomes an element in the constraint
        variables       = []
        stoichiometries = []
        for reaction in model.reaction.get_by_contains(metabolite):
            variables.append(reaction.lp_var)
            stoichiometries.append(reaction.participants[metabolite])

        if metabolite.boundaryCondition:
            # by the nomenclature of SBML, being a boundary metabolite means
            # the net production is *not* constrained
            # metabolite.lp_constr = None
            pass
        else:
            # the weighted sum of all the reactions featuring the metabolite
            # must have a net total of zero
            constr = model.lp.addConstr(grb.LinExpr(stoichiometries, variables),
                                        GRB.EQUAL, 0.0, model.name + metabolite.id)
            metabolite.lp_constr = constr
    model.lp.update()
