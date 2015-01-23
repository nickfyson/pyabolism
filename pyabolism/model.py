import sys
from collections import OrderedDict

class Unit(object):
    """docstring for Unit"""
    def __init__(self, arg, **kwargs):
        self.kind       = arg
        self.scale      = kwargs.get('scale',0.0)
        self.exponent   = kwargs.get('exponent',1.0)
        self.offset     = kwargs.get('offset',0.0)
        self.multiplier = kwargs.get('multiplier',1.0)
    
    def __str__(self):
        return self.id
    
    def __repr__(self):
        return self.id
    


class UnitDefinition(object):
    """docstring for Unit"""
    
    def __init__(self, arg):
        self.id    = arg
        self.units = []
    
    def __str__(self):
        return self.id
    
    def __repr__(self):
        return self.id
    


class Compartment(object):
    """docstring for Compartment"""
    def __init__(self, arg, **kwargs):
        self.id      = arg
        self.name    = kwargs.get('name',None)        
        self.outside = kwargs.get('name',None)
    
    def __str__(self):
        return self.name
    
    def __repr__(self):
        return self.name
    


class Metabolite(object):
    """docstring for Metabolite"""
    def __init__(self, arg, **kwargs):
        self.id                 = arg
        self.name               = kwargs.get('name',self.id)
        self.readable           = ''
        self.compartment        = kwargs.get('compartment',None)
        self.formula            = kwargs.get('formula',None)
        self.charge             = kwargs.get('charge',None)
        self.boundaryCondition  = kwargs.get('boundaryCondition',False)
        self.participations     = OrderedDict()
        self.notes              = {}
    
    def __str__(self):
        return self.id
    
    def __repr__(self):
        return self.id
    


class Reaction(object):
    """docstring for Reaction"""
    def __init__(self, arg, **kwargs):
        self.id                     = arg
        self.name                   = kwargs.get('name',self.id)
        self.reversible             = kwargs.get('reversible',None)
        self.participants           = OrderedDict()
        self.lower_bound            = -1e4
        self.upper_bound            = 1e4
        self.default_bounds         = (self.lower_bound,self.upper_bound)
        self.objective_coefficient  = 0.0
        self.flux_value             = None
        self.notes                  = {}
        self.genes                  = []
    
    def __str__(self):
        reactants = []
        products  = []
        for p in self.participants:
            if self.participants[p] < 0:
                reactants.append( ' %03.1f %s '%(abs(self.participants[p]),p.name) )
            else:
                products.append( ' %03.1f %s '%(abs(self.participants[p]),p.name) )
        if self.reversible:
            arrow = '<===>'
        else:
            arrow = '====>'
        return self.id+' : '+'+'.join(reactants)+'   '+arrow+'   '+'+'.join(products)+'\n'
    
    def __repr__(self):
        return self.id

    def clear_participants(self):
        """docstring for clear_participants"""
        for metabolite in self.participants:
            metabolite.participations.pop(self.id)
        self.participants.clear()
    
    def remove_participant(self,metabolite):
        """docstring for remove_participant"""
        metabolite.participations.pop(self.id)
        self.participants.pop(metabolite)
    
    def add_participant(self,metabolite,stoichiometry):
        """docstring for add_participant"""
        self.participants[metabolite]      = stoichiometry
        metabolite.participations[self.id] = stoichiometry
    
    def set_default_bounds(self):
        """docstring for set_default_bounds"""
        self.default_bounds = (self.lower_bound,self.upper_bound)
    
    def reset_bounds(self):
        """docstring for reset_bounds"""
        self.lower_bound = self.default_bounds[0]
        self.upper_bound = self.default_bounds[1]
    


class Gene(object):
    """docstring for Gene"""
    def __init__(self, arg, **kwargs):
        self.id         = arg
        self.name       = kwargs.get('name','')
        self.expression = kwargs.get('expression',1)
    
    def __str__(self):
        return str(self.expression > 0)
    
    def __repr__(self):
        return str(self.expression > 0)
    

class _CompartmentDict(OrderedDict):
    """docstring for Compartments"""
    def __init__(self, *arg, **kwargs):
        super(_CompartmentDict, self).__init__(*arg,**kwargs)
    
    def add(self,compartment):
        """docstring for add"""
        if compartment.id in self:
            sys.exit('Error! The compartment id %s already exists!'%compartment.id)
        self[compartment.id] = compartment
    
    def remove(self,compartment):
        """docstring for remove"""
        self.pop(compartment.id)
    


class _MetaboliteDict(OrderedDict):
    """docstring for Metabolites"""
    def __init__(self, *arg, **kwargs):
        super(_MetaboliteDict, self).__init__(*arg,**kwargs)
    
    def add(self,metabolite):
        """docstring for _add_metabolite"""
        if metabolite.id in self:
            raise Exception('Error! The metabolite id %s already exists!'%metabolite.id)
        self[metabolite.id] = metabolite
    
    def remove(self,metabolite):
        """docstring for remove"""
        self.pop(metabolite.id)

class _ReactionDict(OrderedDict):
    """docstring for Reactions"""
    def __init__(self, *arg, **kwargs):
        super(_ReactionDict, self).__init__(*arg,**kwargs)
    
    def add(self,reaction):
        """docstring for _add_metabolite"""
        if reaction.id in self:
            raise Exception('Error! That reaction id %s already exists!'%reaction.id)
        self[reaction.id] = reaction
        for metabolite in reaction.participants:
            metabolite.participations[reaction.id] = reaction.participants[metabolite]
    
    def remove(self,reaction):
        """docstring for remove"""
        reaction.clear_participants()
        self.pop(reaction.id)

    def get_by_contains(self,metabolite):
        """docstring for contains"""
        return [self[r_id] for r_id in metabolite.participations]
    
    def get_by_consumes(self,metabolite):
        """docstring for contains"""
        return [self[r_id] for r_id in metabolite.participations if metabolite.participations[r_id] < 0]
    
    def get_by_produces(self,metabolite):
        """docstring for contains"""
        return [self[r_id] for r_id in metabolite.participations if metabolite.participations[r_id] > 0]


class _GeneDict(OrderedDict):
    """docstring for Genes"""
    def __init__(self, *arg, **kwargs):
        super(_GeneDict, self).__init__(*arg,**kwargs)
    
    def add(self,gene):
        """docstring for add"""
        if gene.id in self:
            sys.exit('Error! The gene id %s already exists!'%gene.id)
        self[gene.id] = gene
    
    def remove(self,gene):
        """docstring for remove"""
        self.pop(gene.id)
    


class MetaModel(object):
    """docstring for MetaModel"""
    def __init__(self, **kwargs):
        self.id               = kwargs.get('id',None)
        self.name             = kwargs.get('name',None)
        self.metabolite       = kwargs.get('metabolites',_MetaboliteDict())
        self.reaction         = kwargs.get('reactions',_ReactionDict())
        self.compartment      = kwargs.get('compartments',_CompartmentDict())
        self.gene             = kwargs.get('genes',_GeneDict())
        self.unit_definition  = kwargs.get('unit_definitions',OrderedDict())   

    def metabolites(self):
        """docstring for metabolites"""
        return self.metabolite.values()

    def reactions(self):
        """docstring for metabolites"""
        return self.reaction.values()

    def compartments(self):
        """docstring for metabolites"""
        return self.compartment.values()

    def genes(self):
        """docstring for metabolites"""
        return self.gene.values()

    def unit_definitions(self):
        """docstring for unit_definitions"""
        return self.unit_definition.values()