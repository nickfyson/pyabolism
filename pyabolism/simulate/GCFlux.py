
from copy import copy, deepcopy

import numpy as np
import networkx as nx

from ..model import Reaction, Gene
from ..tools import get_transport_reactions, GPR_string2tree

from .LP import grb, GRB, generate_basic_lp

from .simtools import irreversify, deirreversify


def _get_sufficient_complexes(gene_association):
    """from a complex gene association string,
    retrieve a list of all feasible gene complexes (as a list)"""

    gprTree = GPR_string2tree(gene_association)

    for node in reversed(nx.topological_sort(gprTree)):

        if gprTree.out_degree(node) == 0:
            complexes = [[node, ], ]

        elif gprTree.node[node].get('operation', '') == 'or':

            complexes = []
            for child in gprTree.successors(node):
                for c_element in gprTree.node[child]['complexes']:
                    complexes.append(c_element)

        elif gprTree.node[node].get('operation', '') == 'and':

            for i, child in enumerate(gprTree.successors(node)):

                if i == 0:
                    complexes = gprTree.node[child]['complexes']

                else:
                    new_complexes = []
                    for j in range(len(complexes)):
                        for k in range(len(gprTree.node[child]['complexes'])):
                            new_complexes.append(complexes[j] + gprTree.node[child]['complexes'][k])

                    complexes = new_complexes

        else:
            if len(gprTree.successors(node)) > 1:
                raise Exception('branching internal node with no operation label! This is wrong!!')
            child     = gprTree.successors(node)[0]
            complexes = gprTree.node[child]['complexes']

        gprTree.node[node]['complexes'] = complexes

    return gprTree.node['root']['complexes']


def _convert_model(model):

    model = irreversify(model)

    # to ensure that all orginal flux constraints present in the model
    # are preserved in the GC-Flux representation, we keep a record of
    # which reactions are duplicated, for later use
    duplicates_and_bounds = []

    # exchanges = get_exchange_reactions(model)

    # now we duplicate reactions such that each gene association clause
    # is attached to a different reaction
    for r in model.reaction.values():

        # # exchanges are left unaltered
        # if r in exchanges:
        #     continue

        gene_association = r.notes.get('GENE_ASSOCIATION', '')

        suff_complexes = _get_sufficient_complexes(gene_association)

        duplicates = []
        # we iterate over all the feasible complexes, each becoming an individual reaction
        for suff_complex in suff_complexes:

            # the duplicates each get a unique identifier
            new_r = Reaction(r.id + '%03d' % len(duplicates))

            new_r.name                   = r.name
            new_r.reversible             = False
            new_r.lower_bound            = 0.0
            new_r.upper_bound            = np.infty
            new_r.default_bounds         = r.default_bounds
            new_r.objective_coefficient  = r.objective_coefficient
            new_r.flux_value             = r.flux_value
            new_r.notes                  = copy(r.notes)
            new_r.participants           = copy(r.participants)

            new_r.notes['GENE_ASSOCIATION'] = '( %s )' % (' and '.join(suff_complex))

            # all genes in the 'and' clause are required, and the list is recorded
            new_r.genes = [Gene(g_string.strip()) for g_string in suff_complex]

            new_r.notes['original_rid']             = r.notes['original_rid']
            new_r.notes['original_rid_multiplier']  = r.notes['original_rid_multiplier']

            model.reaction.add(new_r)
            duplicates.append(new_r)

        # the duplicates set and the original bounds are recorded
        duplicates_and_bounds.append((duplicates, r.lower_bound, r.upper_bound))

        # finally, once the reaction has been replaced by duplicates, we remove the original
        model.reaction.remove(r)

    model.duplicates_and_bounds = duplicates_and_bounds

    return model


def _build_GCFlux_lp(model, expressions, add_to_existing=False, unlimited_transports=False):

    generate_basic_lp(model, add_to_existing=add_to_existing)

    from collections import defaultdict
    # to set the bounds we need to know the list of reactions catalysed by each gene
    gene2reaction = defaultdict(list)
    for r in model.reactions():
        for gid in [gene.id for gene in r.genes]:
            gene2reaction[gid].append(r.id)

    # build list of all transport reactions (external to cytoplasm compartments)
    transports = get_transport_reactions(model)

    for gid, rids in gene2reaction.items():
        if unlimited_transports:
            # in some cases we may to restrict constraints to internal reactions,
            # and leave transports unrestricted
            r_vars = [model.reaction[rid].lp_var
                      for rid in rids if model.reaction[rid] not in transports]
        else:
            # typically, we need to get a list of the gurobi variables for *every*
            # reaction that is catalysed by the gene
            r_vars = [model.reaction[rid].lp_var for rid in rids]

        # each gene in the model becomes a constraint on the associated reactions
        # with upper bound on total flux dictated by expression of the gene (where available)
        model.lp.addConstr(grb.LinExpr(np.ones(len(r_vars)), r_vars),
                           GRB.LESS_EQUAL,
                           expressions.get(gid, np.infty),
                           'sum_flux_%s' % gid)

    # as well as the expression-related bounds, we must use the duplicate sets stored
    # earlier to re-apply the bounds from the original model definition
    for duplicates, lower_bound, upper_bound in model.duplicates_and_bounds:

        # print duplicates[0].id[:-3], lower_bound, upper_bound

        r_vars = [model.reaction[r.id].lp_var for r in duplicates]

        model.lp.addConstr(grb.LinExpr(np.ones(len(r_vars)), r_vars),
                           GRB.GREATER_EQUAL,
                           lower_bound,
                           'duplicate_lower_%s' % duplicates[0].id[:-3])

        model.lp.addConstr(grb.LinExpr(np.ones(len(r_vars)), r_vars),
                           GRB.LESS_EQUAL,
                           upper_bound,
                           'duplicate_upper_%s' % duplicates[0].id[:-3])

    model.lp.update()


def GCFlux(model, expressions, limit_unpaired=False,
           norm='L2', show=False, unlimited_transports=False):
    """implementation of the GC-Flux algorithm
    Gene complex-centric simulation of cellular metabolism"""

    # we make a duplicate model in order that the original remain in tact
    model = _convert_model(model)

    _build_GCFlux_lp(model, expressions)

    # finally, with the complete model we can run the optimization
    model.lp.optimize()

    # we extract results from the linear program and add them to the GC-Flux model
    for reaction in model.reactions():
        try:
            reaction.flux_value  = reaction.lp_var.X
        except Exception as e:
            print model.lp.status
            raise e

    # we don't yet have a unique solution to our problem,
    # but rather an abitrary flux vector from the solution space
    # to deal with this, we generally minimise the norm while maintaining total objective value
    if norm:
        # tax-cab or 'L1' norm
        if norm == 'L1':
            objective = grb.LinExpr()
            for reaction in model.reactions():

                var        = reaction.lp_var
                flux_value = reaction.flux_value

                if reaction.objective_coefficient != 0:
                    # for those reactions that are part of the objective,
                        # we constrain them to have (almost exactly) the same value
                    var.lb = flux_value * (1. - 1e-12)
                    var.ub = np.infty
                else:
                    # all other reactions are permitted to vary between zero and their current value
                    # to perform the L1 norm, we require that all fluxes maintain their sign from
                    # the original solution
                    var.lb = 0.0
                    var.ub = np.infty  # flux_value

                # all fluxes are in the objective function, such that we can minimise the magnitudes
                objective += var

        elif norm == 'L2':
            # euclidean or 'L2' norm
            objective = grb.QuadExpr()
            for reaction in model.reactions():
                var        = reaction.lp_var
                flux_value = reaction.flux_value

                if reaction.objective_coefficient != 0:
                    # only for the objective reactions to we add additional constraints
                    var.lb = flux_value * (1. - 1e-12)
                    var.ub = np.infty
                # all fluxes are in the objective function, such that we can minimise the magnitudes
                objective += var * var
        else:
            raise Exception('Unknown norm type...')

        # we add the objective function that was built in the last section
        model.lp.setObjective(objective)

        # we now wish to minimise the norm of the flux vector
        model.lp.modelSense = GRB.MINIMIZE

        model.lp.optimize()

        count = 0
        # for some problems, we need to relax the BarConvTol parameter
        # in order to achieve an optimal solution
        while model.lp.status == 13:
            model.lp.setParam('BarConvTol', model.lp.getParamInfo('BarConvTol')[2] * 10.0)
            model.lp.optimize()
            count += 1
            if count > 6:
                raise Exception('Error: relaxing BarConvTol failed to permit optimal solution')

        if model.lp.status != 2:
            raise Exception('non-optimal solution... ' + str(model.lp.status))

    # with the optimisation complete, we need to extract the results
    # and add them back to the original model
    # the fluxes from this LP are transferred to the relevant reactions in model
    for reaction in model.reactions():
        try:
            reaction.flux_value  = reaction.lp_var.X
        except Exception as e:
            print model.lp.status
            raise e

    model = deirreversify(model)

    model.total_objective = 0.0
    for reaction in [r for r in model.reactions() if r.objective_coefficient != 0]:
        model.total_objective += reaction.flux_value * reaction.objective_coefficient
        if show:
            print '%s flux = %18.10f\n' % (reaction.id, reaction.flux_value)

    return model
