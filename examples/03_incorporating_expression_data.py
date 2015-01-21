
# We demonstrate how to use the EFlux function in the simulate module...

from pyabolism import io

model = io.load_model('ecoli_core.xml')


# we have a list of some genes found in the ecoli_core model
gene_ids = ['b2926', 'b2925', 'b0008', 'b3734', 'b3735', 'b3736', 'b3737', 'b3731', 'b3732', 'b0767', 'b3738', 'b3739', 'b1702', 'b2914', 'b3919', 'b3916', 'b2587', 'b3117', 'b0722', 'b0723', 'b0720', 'b0721', 'b0726', 'b0727', 'b0724', 'b0728', 'b1380', 'b1136', 'b1241', 'b0474', 'b1417', 'b1416', 'b2276', 'b2277', 'b2278', 'b2279', 'b2465', 'b2463', 'b2779', 'b1276', 'b4025', 'b2975', 'b2976', 'b0729', 'b3603', 'b2416', 'b2417', 'b2415', 'b1479', 'b2287', 'b1779']

# we create some synthetic expression data, but building a dictionary with a random (0,100) number associated with each of the genes
from random import random
expressions = {}
for gid in gene_ids:
    expressions[gid] = 100*random()


# FBA and EFlux are both found in the 'simulate' module
from pyabolism.simulate import FBA,EFlux

# we first determine the growth rate from basic FBA
FBA(model,show=True)

# and then compare it to that predicted by the EFlux algorithm
EFlux(model,expressions,show=True)


# additionally, if the predicted flux distribution of interest,
    # it is import to use a norm, for example the L2 or 'euclidean' norm
EFlux(model,expressions,norm='L2',show=True)
