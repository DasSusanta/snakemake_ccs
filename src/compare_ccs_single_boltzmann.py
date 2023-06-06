# 2021-11-30

# Kiyoto Aramis Tanemura

# Reviewer asked to justify use of Boltzmann weighting of CCS values. Let's compare Boltzmann weighted CCS against lowest energy CCS.

import os
import pandas as pd

# define paths
inpath = 'results/'
outpath = 'results/comparison.csv'

# initialize lists to store results
identifier = [] # metabolite ID, charge, model
single = []
boltzmann = []

dirs = [x for x in os.listdir(inpath) if os.path.isdir(inpath + x)]
for theID in dirs:
    for charge in os.listdir(inpath + theID):
        for themodel in os.listdir(inpath + theID + '/' + charge + '/ensemble/'):
            resultpath = inpath + theID + '/' + charge + '/ensemble/' + themodel + '/'
            if not os.path.isfile(resultpath + 'ccs.txt'):
                continue
            identifier.append([theID, charge, themodel])
            df = pd.read_csv(resultpath + 'results.csv')
            single.append(df.at[0, 'ccs'])
            with open(resultpath + 'ccs.txt', 'r') as f:
                cont = f.readlines()
            boltzmann_ccs_value = cont[0]
            boltzmann_ccs_value = boltzmann_ccs_value.split(':')[1]
            boltzmann.append(float(boltzmann_ccs_value))

ids = [x[0] for x in identifier]
charge = [x[1] for x in identifier]
model = [x[2] for x in identifier]

df = pd.DataFrame({'id': ids, 'charge': charge, 'model': model, 'single_ccs': single, 'boltzmann_ccs': boltzmann}).sort_values(['id', 'charge', 'model']).reset_index(drop = True)

df.to_csv(outpath)
print(df)
