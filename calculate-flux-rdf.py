import os, sys, time
from utils import *


def write_flux_function(path, bins, values):
    with open(path, 'w') as f:
         f.write('R   F\n')
         for r, value in zip(bins, values):
             f.write('%s    %s\n' % (r, value))


path_data = '.'
try:
   os.mkdir(path_data)
except:
   pass



# READ TRAJ
path = '..'
traj = Traj(path)
traj.read_traj_dlpoly()

# CALCULATE FLUX
length = traj.steps - 1
#length = 3
start = time.time()
atoms = ['C', 'O', 'Li', 'K']
pairs = []
for atom in atoms:
    for atom2 in atoms:
        pairs.append( tuple(sorted( (atom, atom2) )) )
pairs = set(pairs)        
traj.determine_flux(freq0 = 0.0, binwidth = 0.01, length = length, pairs = pairs)
end = time.time()
print 'Calculate flux lasts %s second' % (end - start)

# WRITE FLUX AND RDF
for label in traj.rdfs:
    write_flux_function('%s/rdf-label-%s-length-%s.dat' % (path_data, '-'.join(label), length) , traj.bins, traj.rdfs[label])
    write_flux_function('%s/jrad-label-%s-length-%s.dat'   % (path_data, '-'.join(label),  length) , traj.bins, traj.j_rads[label])
    write_flux_function('%s/jtheta-label-%s-length-%s.dat' % (path_data, '-'.join(label), length) , traj.bins, traj.j_thetas[label])
