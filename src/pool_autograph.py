#!/usr/bin/env python

# 2021-08-26

# Kiyoto Aramis Tanemura

# Upon completion of ten iterations of AutoGraph conformational clustering, combine centroids in one output.
# Arg1: path to directory containing various AutoGraph outputs

import os, sys
import pandas as pd
from shutil import copyfile
import json
with open('paths.json', 'r') as f:
    paths = json.load(f)
sys.path.append(paths['AutoGraph_path'])
from AutoGraph import AutoGraph

if len(sys.argv) < 2:
    print('Upon completion of ten iterations of AutoGraph conformational clustering, combine centroids in one output.')
    print('Arg1: path to directory containing various AutoGraph outputs')
    quit()

# format paths
inpath = '/'.join(sys.argv[1].split('/')) + '/'

df = pd.DataFrame()

for i in range(10):
    df = df.append(pd.read_csv(inpath + str(i) + '/cluster_summary.csv', index_col = 0))

centroids = set(df.loc[df['center'] == 1].index)
if not os.path.isdir(inpath + 'centroids/'):
    os.mkdir(inpath + 'centroids/')
for conf_file in centroids:
    copyfile(inpath.replace('clustered', 'optimized') + conf_file, inpath + 'centroids/' + conf_file)
    
df.loc[df['center'] == 1].to_csv(inpath + 'centroids.csv')
