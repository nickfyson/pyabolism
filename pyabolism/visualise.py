

import matplotlib.pyplot as plt
import seaborn as sns

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

    axis = plt.subplot(111)

    for i, (lower, value, upper, category, rid) in enumerate(plot_vals):

            # axis.plot()
        plt.plot([i, i], [lower, upper], color=(0.8, 0.8, 0.8), zorder=0)
        plt.plot([i, i], [0.0, value], colors[category], zorder=1)

    axis.set_xlim(0,len(plot_vals))

    max_magnitude = max([abs(r.flux_value) for r in model.reactions()])

    axis.set_ylim(-max_magnitude, +max_magnitude)

    axis.set_xticks([])

    axis.set_xlabel('reactions')
    axis.set_ylabel('flux value')


    return


def _plot_flux_points(flux_points,axis,linewidth=1.0):
    """docstring for _plot_flux_points"""

    y_min = y_max = 0
    for (index,lower,value,upper,color) in flux_points:
        # if lower and upper:
        axis.plot([index,index],[lower,upper],color=color,alpha=0.2,linewidth=linewidth,zorder=1)

        if value:
            if value < 0:
                # axis.plot([index,index],[value,0],color,linewidth=linewidth,zorder=3)
                axis.plot(index,value,color+'v',linewidth=linewidth,zorder=3,markeredgewidth=0)
                axis.vlines(index,value,0,colors=color,linewidth=linewidth,zorder=3)
            else:
                # axis.plot([index,index],[0,value],color,linewidth=linewidth,zorder=3)
                axis.plot(index,value,color+'^',linewidth=linewidth,zorder=3,markeredgewidth=0)
                axis.vlines(index,0,value,colors=color,linewidth=linewidth,zorder=3)

            y_min = min(value,y_min)
            y_max = max(value,y_max)

    axis.set_ylim([y_min*1.1,y_max*1.1])


def defunct_plot_fluxes(model,axis=None,reaction_labels=False,title='',filename='flux_vis'):
    """docstring for plot_fluxes"""
    exchange = []
    internal = []
    biomass  = []
    obj_inds = []
    y_lim = (0,0)
    for (index,reaction) in enumerate(model.reactions.values()):
        if reaction.reversible == False:
            plot_vals = (0,reaction.flux_value,reaction.upper_bound)
        else:
            plot_vals = (reaction.lower_bound,reaction.flux_value,reaction.upper_bound)
        # y_lim = ( min(y_lim[0],min(plot_vals)) , max(y_lim[1],max(plot_vals)) )
        y_lim = ( min(y_lim[0],plot_vals[1]) , max(y_lim[1],plot_vals[1]) )

        if 'ex_' in reaction.id.lower():
            exchange.append( plot_vals )
        elif 'bm_' in reaction.id.lower() or 'mass' in reaction.name.lower():
            biomass.append( plot_vals )
        else:
            if reaction.objective_coefficient != 0:
                obj_inds.append(len(internal))
            internal.append( plot_vals )

    show = False
    if axis==None:
        show = True
        axis = plot.sub (figsize=(10,5))

    x_exchange = range(0,len(exchange))
    x_internal = range(len(exchange),len(exchange)+len(internal))
    obj_inds   = np.array(obj_inds) +len(exchange)
    x_biomass  = range(len(exchange)+len(internal),len(exchange)+len(internal)+len(biomass))
    axis.hold(True)


    for (x,y) in zip(x_exchange,exchange):
        axis.plot([x,x],[y[0],y[2]],color=(0.9,0.9,0.9),linewidth=1.0,zorder=1)
        axis.plot([x,x],[min(y[1],0),max(y[1],0)],'b',linewidth=1.0,zorder=3)
    for (x,y) in zip(x_biomass,biomass):
        axis.plot([x,x],[y[0],y[2]],color=(0.9,0.9,0.9),linewidth=1.0,zorder=1)
        axis.plot([x,x],[min(y[1],0),max(y[1],0)],'c',linewidth=1.0,zorder=3)
    for (x,y) in zip(x_internal,internal):
        axis.plot([x,x],[y[0],y[2]],color=(0.9,0.9,0.9),linewidth=1.0,zorder=1)
        if x in obj_inds:
            axis.plot([x,x],[min(y[1],0),max(y[1],0)],'r',linewidth=1.0,zorder=3)
            axis.plot(x,y[1],'r*',markersize=10.0,linewidth=1.0,zorder=3)
        else:
            axis.plot([x,x],[min(y[1],0),max(y[1],0)],'g',linewidth=1.0,zorder=3)

    axis.set_ylim((y_lim[0]*1.1,y_lim[1]*1.1))
    axis.set_ylabel('flux')
    if reaction_labels:
        axis.xaxis.set_ticks([np.mean(x_exchange),np.mean(x_internal),np.mean(x_biomass)])
        axis.set_xticklabels(['exchange','internal','biomass'])
    else:
        axis.xaxis.set_ticks([])

    axis.set_title(title)

    if show:
        plot.savefig("%s.pdf"%filename)
        os.system("open %s.pdf"%filename)
        plot.close('all')


