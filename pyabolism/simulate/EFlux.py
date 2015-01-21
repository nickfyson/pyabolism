
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


def EFlux(model,expressions,norm='L2',show=False,unlimited_transports=False):
    """implementation of the EF-Flux algorithm (Colijn et al)"""
    
    for r in model.reactions():
        r.reset_bounds()
 
    exchanges  = get_exchange_reactions(model)
    transports = get_transport_reactions(model)
    
    # every reaction gets a maxiumum capacity, dictated by its GPR string
    capacities = {}
    for r in model.reactions():
        gene_association = r.notes.get('GENE_ASSOCIATION','')

        capacities[r.id] = _get_capacity(gene_association,expressions)
    
    for r in model.reactions():

        if r in exchanges:
            continue
        
        if unlimited_transports and r in transports:
            continue
                
        capacity = capacities[r.id]
        
        # this value for capacity determines the maximum flow possible through the reaction
        r.upper_bound = capacity

        # capacity due to enzyme expression can work in either direction (for reversible reactions)
        if r.reversible:
            r.lower_bound = -capacity
        else:
            r.lower_bound = 0.0
    
    # growth medium may be set in a binary on/off way
        # but we make any activated exchanges unlimited
    for r in exchanges:
        if r.lower_bound < 0.0:
            r.lower_bound = -np.infty
        if r.upper_bound > 0.0:
            r.upper_bound = np.infty
    
    # our problem can now be solved using the standard FBA algorithm
    FBA(model,show=show,norm=norm)
    
    return model
