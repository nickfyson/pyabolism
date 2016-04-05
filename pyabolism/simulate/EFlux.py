
import numpy as np

from .FBA import FBA

from ..tools import get_transport_reactions, get_exchange_reactions, GPR_string2tree

import networkx as nx


def _get_capacity(gene_association, expressions):
    """find the upper bound on a reaction with the given gene_association string"""

    gprTree = GPR_string2tree(gene_association)

    for node in reversed(nx.topological_sort(gprTree)):

        if gprTree.out_degree(node) == 0:
            gprTree.node[node]['capacity'] = expressions.get(node, np.infty)
        else:
            if gprTree.node[node].get('operation', '') == 'or':
                gprTree.node[node]['capacity'] = \
                    sum([gprTree.node[child]['capacity'] for child in gprTree.successors(node)])

            elif gprTree.node[node].get('operation', '') == 'and':
                gprTree.node[node]['capacity'] = \
                    min([gprTree.node[child]['capacity'] for child in gprTree.successors(node)])

            else:
                if len(gprTree.successors(node)) > 1:
                    raise Exception('missing operation instructions!')
                gprTree.node[node]['capacity'] = \
                    [gprTree.node[child]['capacity'] for child in gprTree.successors(node)][0]

    return gprTree.node['root']['capacity']


def EFlux(model, expressions, norm='L2', show=False, unlimited_transports=False):
    """implementation of the EF-Flux algorithm (Colijn et al)"""

    exchanges  = get_exchange_reactions(model)
    transports = get_transport_reactions(model)

    # every reaction gets a maxiumum capacity, dictated by its GPR string
    capacities = {}
    for r in model.reactions():
        gene_association = r.notes.get('GENE_ASSOCIATION', '')

        capacities[r.id] = _get_capacity(gene_association, expressions)

    for r in model.reactions():

        if r in exchanges:
            continue

        if unlimited_transports and r in transports:
            continue

        capacity = capacities[r.id]

        # in order to avoid destroying limits determined by curators of the model
        # we only allow our expression data to constrict existing constraints, not
        # relax them
        if (0.0 < capacity < r.upper_bound):
            # this value for capacity determines the maximum flow possible through the reaction
            r.upper_bound = capacity

        # similarly, for lower bound we only allow expression to change the limit if it
        # shifts the bound closer to 0.0
        if r.lower_bound < -capacity < 0.0:
            r.lower_bound = -capacity

    # our problem can now be solved using the standard FBA algorithm
    FBA(model, show=show, norm=norm)

    return
