
# Version 0.2dev introduced a change to the way data is stored in the model, and the new approach is outlined here.


from pyabolism import io

model = io.load_model('ecoli_core.xml')


# To get all reactions as a list, we call the function 'reactions'...

for reaction in model.reactions()[:4]:
    print reaction


# To access reactions by ID we can use the ordered dictionary 'reaction'...

print model.reaction['R_AKGDH']
print model.reaction['R_ACONT']



# Metabolites can be accessed in the same manner...

for metabolite in model.metabolites()[:4]:
    print metabolite
