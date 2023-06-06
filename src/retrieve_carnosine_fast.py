import pandas as pd
import os

def ret_ccs(infile):
    with open(infile, 'r') as f:
        lines = f.readlines()
    ccsline = lines[0]
    return float(ccsline.split(':')[1])

# record results in lists
ids = []
charges = []
models = []
ccs_vals = []

for i in range(20):
    for charge in [-1, 1]:
        for model in os.listdir('results/carnosine_' + str(i) + '/' + str(charge) + '/ensemble_fast'):
            inpath = 'results/carnosine_' + str(i) + '/' + str(charge) + '/ensemble_fast/' + model + '/'
            ids.append('carnosine_' + str(i))
            charges.append(charge)
            models.append(model)
            ccs_vals.append(ret_ccs(inpath + 'ccs.txt'))

df = pd.DataFrame({'id': ids, 'charge': charges, 'model': models, 'ccs': ccs_vals})
df.to_csv('carnosine_fast.csv')
print(df)

ids = []
charges = []
models = []
ccs_vals = []
which = []

for theid in [x for x in os.listdir('results') if x[:3] == 'pos']:
    for which_workflow in ['fast', 'standard']:
        for model in os.listdir('results/' + theid + '/1/ensemble/'):
            inpath = 'results/' + theid + '/1/ensemble/' + model + '/ccs.txt'
            if which_workflow == 'fast':
                inpath = inpath.replace('ensemble', 'ensemble_fast')
            if not os.path.isfile(inpath):
                continue
            ids.append(theid)
            charges.append(1)
            models.append(model)
            ccs_vals.append(ret_ccs(inpath))
            which.append(which_workflow)

for theid in [x for x in os.listdir('results') if x[:3] == 'neg']:
    for which_workflow in ['fast', 'standard']:
        for model in os.listdir('results/' + theid + '/-1/ensemble/'):
            inpath = 'results/' + theid + '/-1/ensemble/' + model + '/ccs.txt'
            if which_workflow == 'fast':
                inpath = inpath.replace('ensemble', 'ensemble_fast')
            if not os.path.isfile(inpath):
                continue
            ids.append(theid)
            charges.append(-1)
            models.append(model)
            ccs_vals.append(ret_ccs(inpath))
            which.append(which_workflow)

df = pd.DataFrame({'id': ids, 'charge': charges, 'workflow': which, 'model': models, 'ccs': ccs_vals})

exp_df = pd.DataFrame({'id': ['pos' + str(x) for x in range(13)] + ['neg' + str(x) for x in range(10)],
                       'ccs_exp': [165.0156, 141.7, 150.3, 124.9, 139.8, 184.2, 163.7, 127.6, 143.064906, 151.6, 125.4, 134.7, 189.412707, 127.5, 158.6, 120.7, 150.9, 140.9, 128.5, 122.5, 170.4, 196.3, 178.5]})

df = df.merge(exp_df, how = 'left', on = 'id')

df.to_csv('fast_workflow_evalutation.csv')
print(df)
