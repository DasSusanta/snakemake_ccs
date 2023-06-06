#!/usr/bin/env python
# 2021-08-10

# Kiyoto Aramis Tanemura

# gather centroid energy data
# Arg1: path to cluster_summary.csv

import sys
import pandas as pd

if len(sys.argv) == 1:
    print('gather centroid energy data')
    print('Arg1: path to cluster_summary.csv')
    quit()

clus_path = sys.argv[1].split('/')
identifier = clus_path[1]
adduct = clus_path[2]
model = clus_path[4]

df = pd.read_csv('results/' + identifier + '/' + adduct + '/optimized_conformers/' + model + '/energy.csv', index_col = 0)
centroids = pd.read_csv(sys.argv[1], index_col = 0)
df.loc[centroids.index].to_csv('results/' + identifier + '/' + adduct + '/ensemble_fast/' + model + '/energy.csv')
