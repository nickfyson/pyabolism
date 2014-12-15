
import os,sys

import numpy as np

from .LP import grb, GRB, generate_basic_lp


def FBA(model,show=False,norm=''):
    """builds and solves an FBA linear program, updating the flux_value property of each reaction accordingly"""
    
    lp = generate_basic_lp(model)

    lp.optimize()

    if lp.status == GRB.status.INFEASIBLE:
        if show: print 'model infeasible : no LP solution'
        model.growing = False
        model.total_objective = None
        return
    
    # the fluxes are extracted from the LP and stored in the reactions of the MetaModel
        # preserving the fluxes as found in the model before FBA
    for reaction in model.reactions.values():
        reaction.loaded_flux = reaction.flux_value
        reaction.flux_value  = lp.getVarByName(reaction.id).X
    
    # shadow price of a constraint is a property that can be useful for some analyses
    for metabolite in model.metabolites.values():
        if lp.getConstrByName(metabolite.id):
            metabolite.shadow = lp.getConstrByName(metabolite.id).getAttr('Pi')
        else:
            metabolite.shadow = 0
    
    # in general the flux vector returned by FBA is not unique, so optional minimization of the norm is offered
    if norm:

        # the 'taxicab' norm (L1) can only return a minimised vector where all the signs of the components are unchanged
        if norm == 'L1':
            objective = grb.LinExpr()

            for reaction in model.reactions.values():
                
                var        = lp.getVarByName(reaction.id)
                flux_value = reaction.flux_value

                if reaction.objective_coefficient != 0:
                    # to avoid numerical issues in floating point calculations, we slightly loosen bounds on the objective
                    var.lb = flux_value*(1. - 1e-8)
                    var.ub = np.infty
                else:
                    # if not part of the obective, we require only that the flux remain of the same *sign* as found originally
                    var.lb = min(0,flux_value)
                    var.ub = max(0,flux_value)

                # all fluxes are in the objective function, such that we can minimise the *magnitudes*
                objective += var*float(np.sign(flux_value))
        

        # the euclidean norm requires non-linear objective, but allows the sign of each component to vary freely
        elif norm == 'L2':
            objective = grb.QuadExpr()
            for reaction in model.reactions.values():
                var        = lp.getVarByName(reaction.id)
                flux_value = reaction.flux_value
 
                if reaction.objective_coefficient != 0:
                    # to avoid numerical issues in floating point calculations, we slightly loosen bounds on the objective
                    var.lb = flux_value*(1. - 1e-8)
                    var.ub = np.infty
                
                # the square of every flux is in the objective function
                objective += var*var
        else:
            raise Exception('Unknown norm type...')
        
        # we set the objective function
        lp.setObjective(objective)

        # we now wish to minimise the total (squared) flux
        lp.modelSense = GRB.MINIMIZE
        
        lp.optimize()

        # we store the new fluxes found thanks to the minimisation
        for reaction in model.reactions.values():
            reaction.flux_value = lp.getVarByName(reaction.id).X
    
    # we store the total objective achieved as a property of the model
    model.total_objective = 0
    for reaction in [reaction for reaction in model.reactions.values() if reaction.objective_coefficient != 0]:
            model.total_objective += reaction.flux_value * reaction.objective_coefficient
            
            # where required, the flux through targeted reactions is output to the console
            if show: print '%s flux = %18.10f'%(reaction.id,reaction.flux_value)

    return