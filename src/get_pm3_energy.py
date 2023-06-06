import numpy as np
import pandas as pd
from ase.io import read
from ase.calculators.mopac import MOPAC

elist = []

for thefile in ['2', '6', '1', '0']:
    geom = read('results/HMDB0000001/protonated/site_screening/' + thefile + '.xyz')
    geom[1].charge = 1
    geom.set_calculator(MOPAC(label = 'mopac', charge = 1, method='PM3', NSPA=60, EPS=78.4))
    elist.append(geom.get_potential_energy())

elist = np.array(elist)
elist = elist / 27.211386245988 * 627.5094740631
elist -= np.min(elist)
print(pd.DataFrame({'energy': elist.tolist()}, index = ['2', '6', '1', '0']))
