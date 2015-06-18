
import os,sys

import numpy as np

from .FBA import FBA

from ..tools import get_transport_reactions, get_exchange_reactions

import networkx as nx

def _get_capacity(gene_association,expressions):
    """find the upper bound on a reaction with the given gene_association string"""
    
    string = gene_association.replace('_',' ').replace('(',' ( ').replace(')',' ) ')

    tokenized = string.split()

    graph = nx.DiGraph()
    graph.add_node('root')
    parent = 'root'

    for element in tokenized:

        if element == '(':
            old_parent = parent
            parent     = 'node_%d'%len(graph.nodes())

            graph.add_node(parent)
            graph.add_edge(old_parent,parent)

        elif element == ')':
            parent = graph.predecessors(parent)[0]

        elif element.lower() in ['and','or']:
            if 'operation' in graph.node[parent]:
                graph.node[parent]['operation'].append(element.lower())
            else:
                graph.node[parent]['operation'] = [element.lower()]
        else:
            graph.add_edge(parent,element)

    for node in graph.nodes():
        if 'operation' in graph.node[node]:
            unique_ops = set(graph.node[node]['operation'])
            if len(unique_ops) > 1:
                raise Exception('non-unique operators within a bracket - ambiguous statement!')
            graph.node[node]['operation'] = unique_ops.pop()

    for node in reversed(nx.topological_sort(graph)):

        if graph.out_degree(node) == 0:
            graph.node[node]['capacity'] = expressions.get(node,np.infty)
        else:
            if graph.node[node].get('operation','') == 'or':
                graph.node[node]['capacity'] = sum([graph.node[child]['capacity'] for child in graph.successors(node)])

            elif graph.node[node].get('operation','') == 'and':
                graph.node[node]['capacity'] = min([graph.node[child]['capacity'] for child in graph.successors(node)])

            else:
                if len(graph.successors(node)) > 1:
                    raise Exception('missing operation instructions!')
                graph.node[node]['capacity'] = [graph.node[child]['capacity'] for child in graph.successors(node)][0]
    
    
    return graph.node['root']['capacity']


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
