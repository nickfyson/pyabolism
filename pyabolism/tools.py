
def get_exchange_reactions(model):
    """exchange reactions are those that convert boundary metabolites into external metabolites"""
    exchanges = []
    
    for r in model.reactions():

        if [m for (m, s) in r.participants.items() if m.boundaryCondition]:
            exchanges.append(r)
            continue

        if len(r.participants) == 1:
            exchanges.append(r)
            continue

    return exchanges


def get_transport_reactions(model):
    """return all reactions that move metabolites across compartment boundaries"""
    transports = []
    for r in model.reactions():
        if set(['e', 'c']) == set([m.compartment for m in r.participants.keys()]):
            transports.append(r)
    return transports


def GPR_string2tree(gene_association):
    """parse a standard format GPR string into an annotated tree"""
    import networkx as nx

    string = gene_association.replace('_', ' ').replace('(', ' ( ').replace(')', ' ) ')

    tokenized = string.split()

    graph = nx.DiGraph()
    graph.add_node('root')
    parent = 'root'

    for element in tokenized:

        if element == '(':
            old_parent = parent
            parent     = 'node_%d' % len(graph.nodes())

            graph.add_node(parent)
            graph.add_edge(old_parent, parent)

        elif element == ')':
            parent = graph.predecessors(parent)[0]

        elif element.lower() in ['and', 'or']:
            if 'operation' in graph.node[parent]:
                graph.node[parent]['operation'].append(element.lower())
            else:
                graph.node[parent]['operation'] = [element.lower()]
        else:
            graph.add_edge(parent, element)

    for node in graph.nodes():
        if 'operation' in graph.node[node]:
            unique_ops = set(graph.node[node]['operation'])
            if len(unique_ops) > 1:
                raise Exception('non-unique operators within a bracket - ambiguous statement!')
            graph.node[node]['operation'] = unique_ops.pop()

    return graph


def find_config_folder():
    """searches path from cwd to $HOME, searching for .pyabolism directory
        if cwd not within home directory, searches cwd only
        if finds nothing, returns None"""

    from os import getcwd
    from os.path import sep, expanduser, isdir

    folder_name = '.pyabolism'
    
    home_path = expanduser('~')
    cwd       = getcwd()

    if cwd[0:len(home_path)] != home_path:

        test_path = sep.join([cwd, folder_name])
        if isdir(test_path):
            return test_path

    else:
    
        split_home_path = home_path.split(sep)
        split_cwd       = cwd.split(sep)

        for i in range(len(split_cwd) - len(split_home_path) + 1):
            
            test_path = sep.join(split_cwd[:len(split_cwd) - i] + [folder_name])

            if isdir(test_path):
                return test_path

    return None
