import numpy as np
from ase.io import read, write 
vaspobj = read('vasprun.xml', index=":")
confnam = "liquid"
for ind, ib in enumerate(vaspobj):
    xx,yy,zz,yz,xz,xy=-vaspobj[ind].calc.results['stress']*vaspobj[ind].get_volume()
    vaspobj[ind].info['virial']= np.array([(xx, xy, xz), (xy, yy, yz), (xz, yz, zz)])
    del vaspobj[ind].calc.results['stress']
    vaspobj[ind].pbc=True
    vaspobj[ind].info['config_type']=confnam
    write(f"{confnam}.xyz",vaspobj[ind],append=True)
