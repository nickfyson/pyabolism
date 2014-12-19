
def get_exchange_reactions(model):
    """exchange reactions are those that convert boundary metabolites into external metabolites"""
    exchanges = []
    
    for r in model.reactions():

        if [m for (m,s) in r.participants.items() if m.boundaryCondition]:
            exchanges.append( r )
            continue

        if len(r.participants)==1:
            exchanges.append( r )
            continue

    return exchanges

def get_transport_reactions(model):
    """return all reactions that move metabolites across compartment boundaries"""
    transports = []
    for r in model.reactions():
        if set(['e','c'])==set([m.compartment for m in r.participants.keys()]):
            transports.append( r ) 
    return transports


