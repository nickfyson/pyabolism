
import os
import re
from libsbml import SBMLReader, SBMLWriter, SBMLDocument, UnitKind_toString, UnitKind_forName

from .model import MetaModel, Compartment, Metabolite, Reaction, Unit, UnitDefinition


def load_model(filename, filetype='sbml'):

    if filetype.lower() == 'sbml':
        return load_sbml(filename)
    elif filetype.lower() == 'pickle':
        return load_pickle(filename)
    else:
        raise Exception('Unknown filetype!')


def save_model(model, filename, filetype='sbml'):

    if filetype.lower() == 'sbml':
        save_sbml(model, filename)
    elif filetype.lower() == 'pickle':
        save_pickle(model, filename)
    else:
        raise Exception('Unknown filetype!')


def load_sbml(filename):
    """docstring for load_model"""

    notes_pattern = re.compile('\<\w*:?\w*\>([^<>]*)\<\/\w*:?\w*\>')

    sbml_document = SBMLReader().readSBML(filename)
    sbml_model    = sbml_document.getModel()

    if sbml_model is None:
        raise IOError('Failed to load model.')

    # initialise model

    model = MetaModel(id=sbml_model.getId())
    model.name = sbml_model.getName()

    # add all compartments
    for sbml_compartment in sbml_model.getListOfCompartments():
        compartment = Compartment(sbml_compartment.getId(), name=sbml_compartment.getName(),
                                  outside=sbml_compartment.getOutside())
        model.compartment.add(compartment)

    # add all the unit definitions
    for sbml_unitDefinition in sbml_model.getListOfUnitDefinitions():
        unit_definition = UnitDefinition(sbml_unitDefinition.getId())
        for sbml_unit in sbml_unitDefinition.getListOfUnits():
            unit = Unit(UnitKind_toString(sbml_unit.getKind()))
            unit.multiplier = sbml_unit.getMultiplier()
            unit.scale      = sbml_unit.getScale()
            unit.exponent   = sbml_unit.getExponent()
            unit.offset     = sbml_unit.getOffset()
            unit_definition.units.append(unit)
        model.unit_definition[unit_definition.id] = unit_definition

    # add all species
    for sbml_species in sbml_model.getListOfSpecies():

        metabolite = Metabolite(sbml_species.getId(),
                                name=sbml_species.getName(),
                                compartment=sbml_species.getCompartment(),
                                charge=sbml_species.getCharge(),
                                boundaryCondition=sbml_species.getBoundaryCondition())

        metabolite.raw_notes = sbml_species.getNotesString()
        model.metabolite.add(metabolite)

    # add all reactions
    for sbml_reaction in sbml_model.getListOfReactions():

        reaction = Reaction(sbml_reaction.getId(), name=sbml_reaction.getName(),
                            reversible=sbml_reaction.getReversible())

        for string in notes_pattern.findall(sbml_reaction.getNotesString()):
            reaction.notes[string.split(':')[0].strip()] = string.split(':')[-1].strip()

        reaction.lower_bound = \
            sbml_reaction.getKineticLaw().getParameter('LOWER_BOUND').getValue() \
            if sbml_reaction.getKineticLaw().getParameter('LOWER_BOUND') else None
        reaction.upper_bound = \
            sbml_reaction.getKineticLaw().getParameter('UPPER_BOUND').getValue() \
            if sbml_reaction.getKineticLaw().getParameter('UPPER_BOUND') else None

        reaction.default_bounds = (reaction.lower_bound, reaction.upper_bound)

        reaction.objective_coefficient = \
            sbml_reaction.getKineticLaw().getParameter('OBJECTIVE_COEFFICIENT').getValue() \
            if sbml_reaction.getKineticLaw().getParameter('OBJECTIVE_COEFFICIENT') else 0.0

        for sbml_reactant in sbml_reaction.getListOfReactants():

            try:
                metabolite      = model.metabolite[sbml_reactant.getSpecies()]
            except KeyError:
                m_id = sbml_reactant.getSpecies()
                metabolite = Metabolite(m_id,
                                        compartment=m_id.split('_')[-1],
                                        boundaryCondition=False)
                model.metabolite.add(metabolite)
            stoichiometry   = -1.0 * sbml_reactant.getStoichiometry()

            reaction.add_participant(metabolite, stoichiometry)

        for sbml_reactant in sbml_reaction.getListOfProducts():

            try:
                metabolite      = model.metabolite[sbml_reactant.getSpecies()]
            except KeyError:
                m_id = sbml_reactant.getSpecies()
                metabolite = Metabolite(m_id,
                                        compartment=m_id.split('_')[-1],
                                        boundaryCondition=False)
                model.metabolite.add(metabolite)
            stoichiometry   = 1.0 * sbml_reactant.getStoichiometry()

            reaction.add_participant(metabolite, stoichiometry)

        model.reaction.add(reaction)

    return model


def save_sbml(model, filename):
    """docstring for save_model"""

    sbml_document = SBMLDocument(2, 1)
    sbml_model    = sbml_document.createModel(model.id)

    sbml_model.getNamespaces().add("http://www.w3.org/1999/xhtml", "html")

    if model.name:
        sbml_model.setName(model.name)

    for compartment in model.compartments():
        sbml_compartment = sbml_model.createCompartment()
        sbml_compartment.setId(compartment.id)
        sbml_compartment.setName(compartment.name)

    for unit_definition in model.unit_definitions():
        sbml_unitDefinition = sbml_model.createUnitDefinition()
        sbml_unitDefinition.setId(unit_definition.id)
        for unit in unit_definition.units:
            sbml_unit = sbml_unitDefinition.createUnit()
            sbml_unit.setKind(UnitKind_forName(unit.kind))
            sbml_unit.setMultiplier(unit.multiplier)
            sbml_unit.setOffset(unit.offset)
            sbml_unit.setExponent(unit.exponent)
            sbml_unit.setScale(unit.scale)

    for metabolite in model.metabolites():
        species = sbml_model.createSpecies()
        species.setId(metabolite.id)
        species.setName(metabolite.name)
        if metabolite.charge:
            species.setCharge(metabolite.charge)
        species.setCompartment(metabolite.compartment)
        species.setBoundaryCondition(metabolite.boundaryCondition)
        if metabolite.notes:
            for (key, value) in metabolite.notes.items():
                species.appendNotes('\n<html:p>%s: %s</html:p>' % (key, value))
            species.appendNotes('\n')

    for reaction in model.reactions():
        sbml_reaction = sbml_model.createReaction()
        sbml_reaction.setId(reaction.id)
        sbml_reaction.setName(reaction.name)
        sbml_reaction.setReversible(reaction.reversible)

        for (key, value) in reaction.notes.items():
            sbml_reaction.appendNotes('\n<html:p>%s: %s</html:p>' % (key, value))
        sbml_reaction.appendNotes('\n')

        kineticLaw = sbml_reaction.createKineticLaw()
        kineticLaw.setFormula('FLUX_VALUE')

        sbml_lower = kineticLaw.createParameter()
        sbml_lower.setId('LOWER_BOUND')
        sbml_lower.setValue(reaction.lower_bound)

        sbml_upper = kineticLaw.createParameter()
        sbml_upper.setId('UPPER_BOUND')
        sbml_upper.setValue(reaction.upper_bound)

        sbml_upper = kineticLaw.createParameter()
        sbml_upper.setId('OBJECTIVE_COEFFICIENT')
        sbml_upper.setValue(reaction.objective_coefficient)

        if reaction.flux_value:
            sbml_upper = kineticLaw.createParameter()
            sbml_upper.setId('FLUX_VALUE')
            sbml_upper.setValue(reaction.flux_value)

        for (metabolite, coefficient) in reaction.participants.items():
            if coefficient < 0:
                speciesReference = sbml_reaction.createReactant()
                coefficient      = abs(coefficient)
            else:
                speciesReference = sbml_reaction.createProduct()
            speciesReference.setSpecies(metabolite.id)
            speciesReference.setStoichiometry(float(coefficient))

    writer = SBMLWriter()
    writer.writeSBML(sbml_document, filename)


def uncache_model(bug_name):
    """Searches .pyabolism file for matching bug name
        if available will load pickle file
        otherwise, will load from SBML and store pickle for future use

        NB Pyabolism bugs are intended to speed up loading models, not for storage
            of edits made. These should be written back to external SBML files for
            safekeeping.

            Bugs will not in general survive upgrades to Pyabolism code base!
            Use with care!
        """
    raise Exception("Apologies, in it's current state this functionality is best avoided!")

    from .tools import find_config_folder
    config_folder = find_config_folder()

    if not config_folder:
        raise Exception("Unable to load bug, can't find a config folder!")

    import pickle
    from os.path import sep

    try:
        return pickle.load(open(sep.join([config_folder, 'bugs', '%s.pickle' % bug_name]), 'r'))
    except IOError:
        raise IOError('Sorry, unable to find a bug of that name...')


def cache_model(model, bug_name, overwrite=False):
    """cache a model as pickle inside the pyabolism config folder

        NB Pyabolism bugs are intended to speed up loading models, not for storage
            of edits made. These should be written back to external SBML files.

            Bugs will not in general survive upgrades to Pyabolism code base!
            Use with care!!"""

    raise Exception("Apologies, in it's current state this functionality is best avoided!")

    from .tools import find_config_folder
    config_folder = find_config_folder()

    if not config_folder:
        raise Exception("Unable to save bug, can't find a config folder!")

    from os.path import isdir, sep, isfile

    if not isdir(sep.join([config_folder, 'bugs'])):
        os.mkdir(sep.join([config_folder, 'bugs']))

    pickle_name = sep.join([config_folder, 'bugs', '%s.pickle' % bug_name])

    if isfile(pickle_name) and not overwrite:
        raise Exception('Bug already exists! Pass overwrite=True to replace existing pickle.')

    save_pickle(model, open(pickle_name, 'w'))

    return


import pickle


def save_pickle(model, filename):

    model.lp = None
    for r in model.reactions():
        r.lp_var = None
    for m in model.metabolites():
        m.lp_constr = None

    pickle.dump(model, open(filename, 'w'))


def load_pickle(filename):

    return pickle.load(open(filename, 'r'))
