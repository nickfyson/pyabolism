
import os,sys

import numpy as np

from .FBA import FBA

from ..tools import get_transport_reactions, get_exchange_reactions


def _get_capacity(gene_association,expressions):
    """find the upper bound on a reaction with the given gene_association string"""
    
    gene_association = gene_association.replace('(',' ').replace(')',' ')
    
    # we refer to each set of genes sufficient to catalyse the reaction as a 'clause'
    clauses = gene_association.split(' or ')
    
    capacity = 0
    for clause in clauses:
        # within each clause are a number of genes
        genes = [gene.strip() for gene in clause.split(' and ')]
        # for each of these, the contribution is limited by the lowest capacity
        capacity += np.min([expressions.get(gene,np.infty) for gene in genes])

    return capacity


def EFlux(model,expressions,limit_unpaired='',norm='L2',show=False,unlimited_transports=False):

    # we compile a dictionary of capacities for all reactions in the model
    capacities = {}
    for r in model.reactions():
        gene_association = r.notes.get('GENE_ASSOCIATION','')

        capacities[r.id] = _get_capacity(gene_association,expressions)
    
    # in situations where we limit reactions that are *not* limited by expression
        # we may wish to cap reaction flux in some other way
    if limit_unpaired:
        # limit all fluxes by the maximum of the determined bounds
        if limit_unpaired=='maximum':
            default_bound = np.max( [v for k,v in capacities.items() if v] )
        # limit by median of determined bounds
        elif limit_unpaired=='median':
            default_bound = np.median( [v for k,v in capacities.items() if v] )
        else:
            sys.exit('Unknown value for limit_unpaired argument')
    # if no limit type chosen, the maximum flux for unbounded reactions is infinite
    else:
        default_bound = np.infty
    
    # the bounds on each reaction are set in turn
    for r in model.reactions():

        capacity = capacities[r.id]
        
        # if expression-inspired limits are not available in this case...
        if not capacity:
            # we use the default_bound determined earlier
            capacity = default_bound
        
        # reversible reactions use the same enzyme, and hence get the same bounds
        if r.reversible:
            r.lower_bound = -capacity
        # irreversible reactions cannot have negative flux under any circumstances
        else:
            r.lower_bound = 0.0
        
        r.upper_bound = capacity
    
    # make any activated exchanges and transports unlimited
        # while respecting existing directionality of bounds
    for r in get_exchange_reactions(model):
        if r.lower_bound < 0.0 and r.reversible:
            r.lower_bound = -np.infty
        if r.upper_bound > 0.0:
            r.upper_bound = np.infty
    
    if unlimited_transports:
        for r in get_transport_reactions(model):
            if r.lower_bound < 0.0 and r.reversible:
                r.lower_bound = -np.infty
            if r.upper_bound > 0.0:
                r.upper_bound = np.infty
    # with the reaction bounds appropriately set, the problem is just standard FBA
    FBA(model,show=show,norm=norm)
    
    return model
