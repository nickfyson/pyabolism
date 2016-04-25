
from copy import copy, deepcopy

from ..model import Reaction


def irreversify(original):

    original.lp = None
    for r in original.reactions():
        r.lp_var = None
    for m in original.metabolites():
        m.lp_constr = None

    model = deepcopy(original)

    model.original = original

    # our approach requires that reversible reactions be treated separately,
    # in forward and backward direction
    for r in model.reactions():

        # a new reaction is created
        new_r = Reaction(r.id + '_f')

        new_r.name                   = r.name
        new_r.reversible             = False  # all reactions in our altered model are irreversible
        new_r.lower_bound            = 0.0
        new_r.upper_bound            = r.upper_bound
        # all reactions in our altered model are irreversible
        new_r.default_bounds         = (0.0, r.default_bounds[1])
        new_r.objective_coefficient  = r.objective_coefficient
        new_r.flux_value             = r.flux_value
        new_r.notes                  = copy(r.notes)
        new_r.genes                  = copy(r.genes)

        new_r.participants           = copy(r.participants)

        # for retrieving the eventual results, we require the following information
        new_r.notes['original_rid']            = r.id
        new_r.notes['original_rid_multiplier'] = 1.0

        model.reaction.add(new_r)

        # if the original reaction is reversible, we add a second (flipped) version to our new model
        if r.reversible:
            new_r = Reaction(r.id + '_b')

            new_r.name                   = r.name
            new_r.reversible             = False  # again, all reactions are irreversible
            new_r.lower_bound            = 0.0
            new_r.upper_bound            = abs(r.lower_bound)
            # the upper bound is the absolute value of the original lower
            new_r.default_bounds         = (0.0, abs(r.default_bounds[0]))
            new_r.objective_coefficient  = r.objective_coefficient
            new_r.flux_value             = r.flux_value
            new_r.notes                  = copy(r.notes)
            new_r.genes                  = copy(r.genes)

            new_r.participants           = copy(r.participants)

            new_r.notes['original_rid']            = r.id
            new_r.notes['original_rid_multiplier'] = -1.0

            # here we negate the stoichiometry of all participants
            for m, s in new_r.participants.items():
                new_r.participants[m] = -1.0 * s

            model.reaction.add(new_r)

        # now we have replaced the original reaciton with duplicate(s) we remove it
        model.reaction.remove(r)

    return model


def deirreversify(model):

    original = model.original

    # we now sum the fluxes of each reaction in the model to the relevant reaction in original
    # we first store the loaded value, and set the flux to zero
    for r in original.reactions():
        r.loaded_flux = r.flux_value
        r.flux_value  = 0.0

    # then sum up the flux found in each associated reaction in the GC-Flux model
    for r in model.reactions():
        if isinstance(r.flux_value, float):
            original.reaction[r.notes['original_rid']].flux_value += \
                r.notes['original_rid_multiplier'] * r.flux_value

    return original
