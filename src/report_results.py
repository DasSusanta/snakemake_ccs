# 2021-09-02

# Kiyoto Aramis Tanemura

# The CCS values are stored in ccs.txt, but are located in individual directories. Consolidate the results into one file.

import os
import pandas as pd

filepaths = []

# record all the paths to ccs.txt
for mol_id in [x for x in os.listdir('results/') if os.path.isdir('results/' + x)]:
    if not os.path.isdir('results/' + mol_id):
        continue
    if mol_id[0] != 'M':
        continue
    for charge in os.listdir('results/' + mol_id):
        if not os.path.isdir('results/' + mol_id + '/' + charge + '/ensemble_fast/'):
            continue
        for themodel in os.listdir('results/' + mol_id + '/' + charge + '/ensemble_fast/'):
            if os.path.isfile('results/' + mol_id + '/' + charge + '/ensemble_fast/' + themodel + '/ccs.txt'):
                filepaths.append('results/' + mol_id + '/' + charge + '/ensemble_fast/' + themodel + '/ccs.txt')

id_list = []
charge_list = []
model_list = []
ccs_list = []
uncert_list = []

# record all variables
for thefile in filepaths:
    x, mol_id, charge, x, themodel, x = thefile.split('/')
    with open(thefile, 'r') as f:
        cont = f.readlines()
    ccs, abs_err, rel_err = [float(x.split()[-1]) for x in cont]
    id_list.append(mol_id)
    charge_list.append(charge)
    model_list.append(themodel)
    ccs_list.append(ccs)
    uncert_list.append(abs_err)

df = pd.DataFrame({'charge': charge_list, 'model': model_list, 'ccs': ccs_list, 'uncertainty': uncert_list}, index = id_list)
df.sort_index(inplace=True)
df.to_csv('results/results.csv')
print(df)
