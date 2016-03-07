

import matplotlib.pyplot as plt
try:
    import seaborn as sns
except ImportError:
    pass

from pyabolism.tools import get_exchange_reactions
from pyabolism.tools import get_transport_reactions


def plot_flux_distribution(model):

    plot_vals = []

    exchanges  = set(get_exchange_reactions(model) + get_transport_reactions(model))

    for r in model.reactions():

        value = r.flux_value
        upper = r.upper_bound

        if r.reversible:
            lower = r.lower_bound
        else:
            lower = 0.0

        if r in exchanges:
            category = 'exchange'
        # elif r in biomass:
        #     category = 'biomass'
        else:
            category = 'main'

        plot_vals.append((lower, value, upper, category, r.id))

    category_order = {'exchange': 1, 'main': 2, 'objective': 3}

    plot_vals.sort(key=lambda x: (category_order[x[-2]], x[-1]))

    colors = {'exchange': 'g', 'main': 'b', 'objective': 'r'}

    fig = plt.figure(figsize=(20, 5))

    sns.set_context("poster")

    axis = fig.add_subplot(111)

    for i, (lower, value, upper, category, rid) in enumerate(plot_vals):
        plt.plot([i, i], [lower, upper], color=(0.8, 0.8, 0.8), zorder=0)
        plt.plot([i, i], [0.0, value], colors[category], zorder=1)

    axis.set_xlim(0, len(plot_vals))

    max_magnitude = max([abs(r.flux_value) for r in model.reactions()])

    axis.set_ylim(-max_magnitude, +max_magnitude)

    axis.set_xticks([])

    axis.set_xlabel('reactions')
    axis.set_ylabel('flux value')

    return


def _plot_flux_points(flux_points, axis, linewidth=1.0):
    """docstring for _plot_flux_points"""

    y_min = y_max = 0
    for (index, lower, value, upper, color) in flux_points:
        # if lower and upper:
        axis.plot([index, index], [lower, upper], color=color,
                  alpha=0.2, linewidth=linewidth, zorder=1)

        if value:
            if value < 0:
                axis.plot(index, value, color + 'v',
                          linewidth=linewidth, zorder=3, markeredgewidth=0)
                axis.vlines(index, value, 0, colors=color, linewidth=linewidth, zorder=3)
            else:
                axis.plot(index, value, color + '^',
                          linewidth=linewidth, zorder=3, markeredgewidth=0)
                axis.vlines(index, 0, value, colors=color, linewidth=linewidth, zorder=3)

            y_min = min(value, y_min)
            y_max = max(value, y_max)

    axis.set_ylim([y_min * 1.1, y_max * 1.1])

