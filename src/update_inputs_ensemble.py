# 2021-12-18

# Kiyoto Aramis Tanemura

# update inputs for unfinished QUICK optimization jobs

import os

completion_keyword = 'Normal Termination'
geometry_keyword = 'GEOMETRY INPUT'

charge = 1

# I be only using the positive adduct for now.
for i in range(13):
    inpath = 'results/pos' + str(i) + '/' + str(charge) + '/ensemble/'
    for model in os.listdir(inpath):
        outfiles = [x for x in os.listdir(inpath + model) if x[-3:] == 'out']
        for thefile in outfiles:
            with open(inpath + model + '/' + thefile, 'r') as f:
                cont = f.read()
            if completion_keyword in cont:
                continue
            print(inpath + model + '/' + thefile)
            os.remove(inpath + model + '/' + thefile)
            if geometry_keyword not in cont:
                continue
            cont = cont.split(geometry_keyword)
            cont = cont[-1]
            lines = cont.split('\n')[1:]
            end_index = lines.index('')
            coords = lines[1:end_index]
            # now we have indices in the xyz format (element x y z)
            # update the input file
            coords = ['DFT B3LYP BASIS=6-31++G(d,p) CUTOFF=1.0d-10 DENSERMS=1.0d-6 OPTIMIZE DIPOLE CHARGE=' + str(charge) + '\n', '\n'] + [x + '\n' for x in coords]
            with open(inpath + model + '/' + thefile.replace('out', 'in'), 'w') as g:
                g.writelines(coords)
