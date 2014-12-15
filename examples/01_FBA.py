


from pyabolism import io,simulate

model = io.load_model('ecoli_core.xml')

simulate.FBA(model,show=True)


