
from collections import defaultdict

import matplotlib.pyplot as plt

from matplotlib.path import Path
import matplotlib.patches as patches

try:
    import seaborn as sns
    sns.set_context("poster")
except ImportError:
    pass

from pyabolism.tools import get_exchange_reactions
from pyabolism.tools import get_transport_reactions


def get_flux_line(reaction, plot_type='value'):
    # flux_value might contain a range
    if plot_type == 'value':
        line_start = 0.0
        line_end   = reaction.flux_value
    if plot_type == 'range':
        line_start = min(reaction.flux_range[0], reaction.flux_range[1])
        line_end   = max(reaction.flux_range[0], reaction.flux_range[1])

    return line_start, line_end


def _extract_axis_limits(reactions):

    breathing_space = 1.1

    # if plot_type == 'value':
    ymin = min([r.flux_value for r in reactions if r.flux_value <= 0.0])
    ymax = max([r.flux_value for r in reactions if r.flux_value >= 0.0])

    return breathing_space * ymin, breathing_space * ymax
    # return ymin, ymax


def plot_flux_distribution(model, reactions=None, title=None):

    exchanges  = set(get_exchange_reactions(model) + get_transport_reactions(model))
    objectives = [r for r in model.reactions() if r.objective_coefficient != 0.0]

    if not reactions:
        reactions = model.reactions()

    for r in reactions:
        if r in exchanges:
            r.category = 'exchange'
        elif r in objectives:
            r.category = 'objective'
        else:
            r.category = 'main'

    reactions.sort(key=lambda x: (x.category, x.id))

    fig = plt.figure(figsize=(20, 5))

    axis = fig.add_subplot(111)

    _plot_bounds(reactions, ax=None)
    try:
        _plot_ranges(reactions, ax=None)
    except AttributeError:
        _plot_fluxes(reactions, ax=None)

    axis.set_xlim(-0.5, len(reactions))
    axis.set_xticks([])
    axis.set_ylim(_extract_axis_limits(reactions))

    axis.set_xlabel('reactions')
    axis.set_ylabel('flux value')

    if title:
        axis.set_title(title)

    return


def _plot_bounds(reactions, ax=None):

    if not ax:
        ax = plt.gca()

    verts = []
    codes = []
    for i, r in enumerate(reactions):
        verts += [(i - 0.4, r.lower_bound),
                  (i - 0.4, r.upper_bound),
                  (i + 0.4, r.upper_bound),
                  (i + 0.4, r.lower_bound),
                  (i - 0.4, r.lower_bound)]

        codes += [Path.MOVETO,
                  Path.LINETO,
                  Path.LINETO,
                  Path.LINETO,
                  Path.CLOSEPOLY]

    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor=(0.8, 0.8, 0.8), lw=0, alpha=0.7)
    ax.add_patch(patch)


def _plot_fluxes(reactions, ax=None):

    if not ax:
        ax = plt.gca()

    colors = {'exchange': 'g', 'main': 'b', 'objective': 'r'}

    verts = defaultdict(list)
    codes = defaultdict(list)

    for i, r in enumerate(reactions):
        verts[r.category] += [(i - 0.25, 0.0),
                              (i - 0.25, r.flux_value),
                              (i + 0.25, r.flux_value),
                              (i + 0.25, 0.0),
                              (i - 0.25, 0.0)]

        codes[r.category] += [Path.MOVETO,
                              Path.LINETO,
                              Path.LINETO,
                              Path.LINETO,
                              Path.CLOSEPOLY]
    for category in verts.keys():

        path = Path(verts[category], codes[category])
        patch = patches.PathPatch(path, facecolor=colors[category], lw=0)
        ax.add_patch(patch)


def _plot_ranges(reactions, ax=None):

    if not ax:
        ax = plt.gca()

    colors = {'exchange': 'g', 'main': 'b', 'objective': 'r'}

    verts = defaultdict(list)
    codes = defaultdict(list)

    for i, r in enumerate(reactions):

        verts[r.category] += [(i - 0.25, r.flux_range[0]),
                              (i - 0.25, r.flux_range[1]),
                              (i + 0.25, r.flux_range[1]),
                              (i + 0.25, r.flux_range[0]),
                              (i - 0.25, r.flux_range[0])]

        codes[r.category] += [Path.MOVETO,
                              Path.LINETO,
                              Path.LINETO,
                              Path.LINETO,
                              Path.CLOSEPOLY]
    for category in verts.keys():

        path = Path(verts[category], codes[category])
        patch = patches.PathPatch(path, facecolor=colors[category], lw=0)
        ax.add_patch(patch)

    fluxes = [r.flux_value for r in reactions]
    ax.scatter(range(len(fluxes)), fluxes, color='r', zorder=3)
