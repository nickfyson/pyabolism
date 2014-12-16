
import os,sys

import numpy as np

from .LP import grb, GRB, generate_basic_lp

from .FBA import FBA

import re
import networkx as nx
def _get_capacity(gene_association,expressions):
    """find the upper bound on a reaction with the given gene_association string"""
    gene_association = gene_association.replace('(',' ').replace(')',' ')
    clauses = gene_association.lower().split(' or ')
    
    capacity = 0
    for clause in clauses:
        genes = [gene.strip() for gene in clause.split(' and ')]

        capacity += np.min([expressions.get(gene,np.infty) for gene in genes])

    return capacity

def _get_exchange_reactions(model):
    """exchange reactions are those that convert boundary metabolites into external metabolites"""
    exchanges = []
    
    for r in model.reactions.values():

        if [m for (m,s) in r.participants.items() if m.boundaryCondition]:
            exchanges.append( r )
            continue

        if len(r.participants)==1:
            exchanges.append( r )
            continue

    return exchanges

def EFlux(model,expressions,limit_unpaired=False,norm='L2',show=False,alt_eflux=False):
    
    exchanges = _get_exchange_reactions(model)
    
    # build list of all transport reactions (external to cytoplasm compartments)
    transports = []
    for r in model.reactions.values():
        if set(['e','c'])==set([m.compartment for m in r.participants.keys()]):
            transports.append( r )
    
    capacities = {}

    for r in model.reactions.values():
        gene_association = r.notes.get('GENE_ASSOCIATION','')

        capacities[r.id] = _get_capacity(gene_association,expressions)
    
    norm_capacity = np.median( [v for k,v in capacities.items() if v] )

    for r in model.reactions.values():

        if r in exchanges or r in transports:
            continue
                
        capacity = capacities[r.id]
        
        if not capacity:
            if limit_unpaired:
                capacity = norm_capacity
            else:
                capacity = np.infty
        
        if r.reversible:
            r.lower_bound = -capacity
        else:
            r.lower_bound = 0.0

        r.upper_bound = capacity
    
    # make any activated exchanges unlimited
    for r in exchanges:
        if r.lower_bound < 0.0:
            r.lower_bound = -np.infty
        if r.upper_bound > 0.0:
            r.upper_bound = np.infty

    FBA(model,show=show,norm=norm)
    
    return model
