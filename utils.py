import re
import numpy as np
import os
import gzip
import scipy.linalg
import time, sys
#import psutil
#from guppy import hpy
#test

def calculate_fft(time, serie):
    serie = np.append(serie, np.zeros(serie.shape))
    fft = np.fft.rfft(serie)
    length = len(fft)
    xaxis = np.arange(length) * (np.pi / time[-1])
    return xaxis, fft

class Molecule(object):

    def __init__(self, label, atoms):
        self.label = label
        self.atoms = atoms


class Traj(object):
    
    def __init__(self, path):
        self.path = path
        self.natoms = 0
        self.mass = 28 
        self.is_read = False
        self.did_cm = False
        self.times = {}
        self.bins = {}
        self.mass = {}

    #@profile
    def read_traj_dlpoly(self, freqstep = 1, printstep = 0):
        if self.is_read:
            print
            'Traj already read once'
        else:
            print "Start to read trajectory"
            with open('%s/HISTORY' % self.path) as f:
                self.labels = {}
                self.forces = []; self.positions = []; self.velocities = []
                list_array = [self.positions, self.velocities, self.forces]
                step = -1; my_step = -1
                for iline, line in enumerate(f):
                    if "timestep" in line:
                        my_step += 1
                        iatom = -1
                        if (my_step % freqstep) == 0:
                            step += 1
                            for array in list_array:
                                array.append([])
                        #if step % 100 == 0:
                        #    print step, psutil.Process(os.getpid()).memory_info()[0]  /1E9
                        #    print sys.getsizeof(self) /1E9
                        #    for array in list_array:
                        #        print sys.getsizeof(array)  /1E9
                        #    h = hpy()
                        #    print h.heap()
                    if iline == 1:
                        self.natoms = int(line.split()[2])
                        print 'Number of atoms: %s' % self.natoms
                    elif iline == 2:
                        self.timestep = float(line.split()[5])
                        print 'Timestep : %s ps' % (self.timestep)
                        self.times['first'] = float(line.split()[1])
                        self.times['saved'] = []
                    elif iline == 3:
                        self.box_length = float(line.split()[0])
                        print 'Box length: %s A' % (self.box_length)
                    if (my_step % freqstep) == 0:
                        if ( np.abs(iline-2) % (4*self.natoms + 4) ) == 0:
                            self.times['saved'].append(float(line.split()[1]) - self.times['first'])
                        if ( np.abs(iline-2) % (4*self.natoms + 4) ) > 3:
                            my_line = (iline - 6 - 4*step) % 4
                            if  my_line == 0:
                                iatom += 1
                                self.labels[iatom] = line.split()[0]
                                for array in [self.forces, self.velocities, self.positions]:
                                    array[step].append([])
                            else:
                                for mod, array in enumerate(list_array):
                                    if (my_line -1) == mod:
                                        array[step][iatom] += map(float, line.split())
                self.steps = len(self.positions)
                self.printed_steps = step
                print 'Steps: %s' % self.steps
                for step in range(self.steps):
                    if len(self.positions[step]) != self.natoms:
                        print "error for step", step, len(self.positions[step])
                self.positions = np.array(self.positions)
                self.velocities = np.array(self.velocities)
                self.forces= np.array(self.forces)
                self.is_read = True
                self.times['saved'] = np.array(self.times['saved'])*self.timestep
                #self.times['printed'] = printstep * np.arange(self.printed_steps)
                #self.times['MD'] = self.timestep * np.arange(self.total_steps)
                self.shift_com()
                self.unwrap()
                self.build_molecules()
                self.calculate_cm()
                print "Trajectory is read"

    def build_molecules(self):
        self.nmols = self.natoms
        self.molecules = {}
        for imol in range(self.nmols):
            self.molecules[imol] = Molecule(label = self.labels[imol], atoms = [imol,])
        self.associate_atoms_mol()

    def associate_atoms_mol(self):
        ii = -1
        self.labels_mol = {}
        self.find_molecules = {}
        for imol in self.molecules:
            self.labels_mol[imol] = self.molecules[imol].label
            for iatom in self.molecules[imol].atoms:
                ii += 1
                self.find_molecules[ii] = self.molecules[imol]



    def unwrap(self):
        #print 0, self.positions[0, 410, :]
        self.positions_wrap = np.copy(self.positions)
        for step in range(1,self.steps):
            #print step, self.positions[step, 410, :]
            vect = self.positions[step, :, :] - self.positions[step -1, : ,: ]
            vect = np.rint(vect/self.box_length)
            self.positions[step, :, :] -= np.multiply(vect,  np.array([[self.box_length,]*3,]*self.natoms))
            #print step, self.positions[step, 410, :], "after"
            if np.isnan(self.positions[step, :, :]).any():
                print "Problem unwrap"
                stop

    def shift_com(self):
        com_vel = np.mean(self.velocities, axis = 0)
        self.vel_shifted = self.velocities - com_vel

    def calculate_cm(self):
        self.positions_cm = np.zeros((self.steps, self.nmols, 3 ))
        self.velocities_cm = np.zeros((self.steps, self.nmols, 3))
        self.forces_cm = np.zeros((self.steps, self.nmols, 3))
        for imol in range(self.nmols):
            mass_mol = 0.0
            for iatom in self.molecules[imol].atoms:
                mass_mol += self.mass.get(self.labels[iatom], 1)
                self.positions_cm[:,imol,:] += self.mass.get(self.labels[iatom], 1) * self.positions[:,iatom,:]
                self.velocities_cm[:, imol, :] += self.mass.get(self.labels[iatom], 1) * self.velocities[:, iatom, :]
                self.forces_cm[:, imol, :] +=  self.forces[:, iatom, :]
            self.positions_cm[:,imol,:] = self.positions_cm[:,imol,:] / mass_mol
            self.velocities_cm[:,imol,:] = self.velocities_cm[:, imol, :]  / mass_mol
        self.did_cm = True


    def calculate_spectral_density(self, length = -1, ):
        if length < 0:
            length = self.steps
        Z_ws = {}
        for atom in set(self.labels.values()):
            Z_ws[atom] = np.zeros(length + 1)
        for i in range(3):
            for iatom in range(self.natoms):
                xaxis, fft = calculate_fft(self.times['saved'][:length], self.vel_shifted[:length, iatom, i])
                Z_ws[self.labels[iatom]] += np.abs(fft)**2
        self.Z_ws = {}
        for atom in set(self.labels.values()):
            self.Z_ws[atom] = (xaxis, Z_ws[atom]  * length / ( np.sum(Z_ws[atom] ) * 3 * self.labels.values().count(atom)))

    def info_matrix(self, matrix):
            rank =  np.linalg.matrix_rank(matrix)
            print 'Rank velvel:', np.linalg.matrix_rank(matrix)
            D, Pvel = scipy.linalg.eigh(matrix)
            is_pos =  np.all(D > 0)
            print 'Is velvel positive definite?', np.all(D > 0)
            print min(D), max(D)
            #velvel = np.matmul(P, np.matmul(np.diag(np.abs(D)), P.T))      

            #plt.subplot(121)
            #plt.hist(np.abs(D),  bins=np.logspace(np.log(min(np.abs(D))),np.log(max(np.abs(D))), 100))
            #plt.xlim([min(np.abs(D)), max(np.abs(D))])
            #plt.xscale('log')

            #plt.subplot(122)
            #plt.plot(np.abs(D), [1/np.sum(np.power(P[:,i], 4)) for i in range(len(D))]  )
            #plt.xlim([min(np.abs(D)), max(np.abs(D))])
            #plt.xscale('log')

            #plt.tight_layout()
            #plt.show()
            
    def calculate_freq_eigen(self, length, start):
            self.velvel = correlation_matrix(self.velocities, length = length, start = start)    
            vv = np.trace(self.velvel) / self.velvel.shape[0]
            print '<v2> = %s' % vv
            temp = self.mass*1E-3 * vv/8.314
            print 'Numerical temperature: %.1f K' % temp
            self.info_matrix(self.velvel)
           
            self.accacc =  correlation_matrix(self.acceleration, length = length, start = start)
            aa = np.trace(self.accacc)/self.accacc.shape[0]
            print '< F2>/M2 = %s' %  aa
            self.info_matrix(self.accacc)            
            print 'Typical frequency: %s ps-1' % (np.sqrt(aa/vv)*1E-12)

            #self.genvp  = scipy.linalg.eigh(self.accacc, b = self.velvel)
            self.genvp  = scipy.linalg.eig(self.accacc, b = self.velvel)

    def count_pair(self, pair):
        result = 1.0
        for atom in pair:
             result = result * self.labels.values().count(atom)
        return result

    def count_pair_mol(self, pair):
        result = 1.0
        for mol in pair:
             result = result * self.labels_mol.values().count(mol)
        return result

    def determine_msd(self, length = -1, atoms = ()):
        if atoms == ():
            atoms = set(self.labels.values())
        if length < 0:
            length = int(self.steps / 2)
        length = min(length, len(self.times['saved']))
        self.msds = {}
        self.times['msd'] = self.times['saved'][:length]
        for atom in atoms:
            self.msds[atom] = np.zeros(length)
        for tau in range(length):
            vect = np.roll(self.positions, -tau, axis=0) - self.positions
            dist = np.sum(np.power(vect, 2), 2)
            for istep in range(self.steps - tau):
                for iatom in range(self.natoms):
                    self.msds[self.labels[iatom]][tau] += dist[istep, iatom]
            for atom in atoms:
                self.msds[atom][tau] = self.msds[atom][tau] / (self.natoms * (self.steps - tau))

    def determine_fluctuations_volume(self, binwidth, atoms = [], nrand = 1):
        start_tot = time.time()
        freqstep = max(1, self.steps / nrand)
        repeat = max(1, nrand/self.steps)
        if atoms == []:
            atoms = set(self.labels.values())
        probas = {}
        self.void = {}
        nbins = (int(self.box_length / (2*binwidth)))
        self.bins['fluct_vol'] = binwidth*np.array(range(1, nbins+1))
        index = 0
        for atom in atoms:
            probas[atom] = np.zeros((nbins, self.natoms))
            self.void[atom] = []
        for rep in range(repeat):
            new_pos = self.positions[::freqstep,:,:]
            rand = np.random.rand(new_pos.shape[0], 3)
            rand = rand*self.box_length
            new_pos = new_pos - rand[:,None,:]
            new_pos += - self.box_length * np.rint(new_pos / self.box_length)
            dist = np.sqrt(np.sum(np.power(new_pos, 2),2))
            idist = np.ceil(dist/binwidth)
            for step in range(idist.shape[0]):
                index += 1
                if (index % 1000) == 0:
                    print "For irand %s finish in:" % index, time.time() - start_tot
                for atom in atoms:
                    natom_pres = 0
                    self.void[atom].append(min(dist[step]))
                    for ibin in range(probas[atom].shape[0]):
                        mask = [x == atom for x in self.labels.values()]
                        bool_list = [d==ibin and lab for (d,lab) in zip(idist[step],mask)]
                        natom_pres += np.sum(bool_list)
                        probas[atom][ibin, natom_pres] += 1
        self.fluct_vol = probas
        print "Total time is :", time.time() - start_tot
        #print probas['C']



    def determine_rdf_cm_forces(self, binwidth, kbT, pairs = [], ):
        if not self.did_cm:
            self.calculate_cm()
        if len(pairs) == 0:
            for mol in set(self.labels_mol.values()):
                for mol2 in set(self.labels_mol.values()):
                    pairs.append( tuple( sorted( (mol, mol2) )))
            pairs = set(pairs)
        start_tot = time.time()
        self.bins['rdf_cm_forces'] = binwidth * np.array(range(int(np.sqrt(2) *self.box_length / (2 * binwidth))))
        self.bins['rdf_cm_forces'] = binwidth * np.array(range(int(self.box_length / binwidth)))
        self.rdf_cm_forces = {}
        count = {}
        for pair in pairs:
            count[pair] = 0
            self.rdf_cm_forces [pair] = np.zeros(len(self.bins['rdf_cm_forces']))
        for imol in range(self.nmols):
            if (imol % 1000) == 0:
                start = time.time()
            for imol2 in [i for i in range(imol + 1, self.nmols)]:
                pair = tuple(sorted([self.molecules[imol].label, self.molecules[imol].label]))
                vect = (self.positions_cm[:, imol2, :] - self.positions_cm[:, imol, :])
                vect = vect - self.box_length * np.rint(vect / self.box_length)
                dist = np.sqrt(np.sum(np.power(vect, 2), 1))
                diff_forces = 0.5*(self.forces_cm[:,imol2, :] - self.forces_cm[:,imol,:])
                dot = (diff_forces * vect).sum(1)
                toadd = dot / dist**3
                for step in range(len(dist)):
                    n =int(dist[step] / binwidth)
                    lim = int((np.sqrt(2) * self.box_length / 2) / binwidth)
                    lim = int(self.box_length/ (2*binwidth))
                    lim = int(self.box_length/binwidth)
                    if  n < lim:
                        count[pair] += 1
                        for k in range(n+1):
                            self.rdf_cm_forces[pair][k] += toadd[step]
                        #for k in range(n,len(self.bins['rdf_forces'])):
                        #    self.rdf_forces[pair][k] += toadd[step]
                    else:
                        print n, lim
            if (imol % 1000) == 0:
                print "For mol %s finish in:" % imol, time.time() - start
        for pair in pairs:
            print pair
            print count[pair]
            print self.count_pair_mol(pair) * self.steps / ( ( 1 + int(pair[0] == pair[1]) ))
            eps = ( 1 + int(pair[0] == pair[1]) )* self.box_length**3 / self.count_pair_mol(pair)
            self.rdf_cm_forces[pair] =  - eps* self.rdf_cm_forces[pair] / (4*np.pi*kbT*self.steps)
        print "Total time is :", time.time() - start_tot

    def determine_rdf_forces(self, binwidth, kbT, pairs = [], ):
        if not self.is_read:
            self.read_traj_dlpoly()
        start_tot = time.time()
        if len(pairs) == 0:
            for atom in set(self.labels.values()):
                for atom2 in set(self.labels.values()):
                    pairs.append( tuple( sorted( (atom, atom2))))
            pairs = set(pairs)
        self.bins['rdf_forces'] = binwidth * np.array(range(int(np.sqrt(2) *self.box_length / (2 * binwidth))))
        self.bins['rdf_forces'] = binwidth * np.array(range(int(self.box_length / binwidth)))
        self.rdf_forces = {}
        count = {}
        for pair in pairs:
            count[pair] = 0
            self.rdf_forces [pair] = np.zeros(len(self.bins['rdf_forces']))
        for iatom in range(self.natoms):
            if (iatom % 1000) == 0:
                start = time.time()
            for iatom2 in [i for i in range(iatom + 1, self.natoms)]:
                pair = tuple(sorted([self.labels[iatom], self.labels[iatom2]]))
                if pair in pairs and iatom2 not in self.find_molecules[iatom].atoms:
                    vect = (self.positions[:, iatom2, :] - self.positions[:, iatom, :])
                    vect = vect - self.box_length * np.rint(vect / self.box_length)
                    dist = np.sqrt(np.sum(np.power(vect, 2), 1))
                    diff_forces = 0.5*(self.forces[:,iatom2, :] - self.forces[:,iatom,:])
                    dot = (diff_forces * vect).sum(1)
                    toadd = dot / dist**3
                    for step in range(len(dist)):
                        n =int(dist[step] / binwidth)
                        lim = int((np.sqrt(2) * self.box_length / 2) / binwidth)
                        lim = int(self.box_length/ (2*binwidth))
                        lim = int(self.box_length/binwidth)
                        if  n < lim:
                            count[pair] += 1
                            for k in range(n+1):
                                self.rdf_forces[pair][k] += toadd[step]
                            #for k in range(n,len(self.bins['rdf_forces'])):
                            #    self.rdf_forces[pair][k] += toadd[step]
                        else:
                            print n, lim
            if (iatom % 1000) == 0:
                print "For atom %s finish in:" % iatom, time.time() - start
        for pair in pairs:
            print pair
            print count[pair]
            print self.count_pair(pair) * self.steps / ( ( 1 + int(pair[0] == pair[1]) ))
            eps = ( 1 + int(pair[0] == pair[1]) )* self.box_length**3 / self.count_pair(pair)
            self.rdf_forces[pair] =  - eps* self.rdf_forces[pair] / (4*np.pi*kbT*self.steps)
        print "Total time is :", time.time() - start_tot


    def determine_rdf(self, binwidth,  pairs = [], wrap = False):
        if not self.is_read:
            self.read_traj_dlpoly()
        if len(pairs) == 0:
            for atom in set(self.labels.values()):
                for atom2 in set(self.labels.values()):
                    pairs.append( tuple( sorted( (atom, atom2))))
            pairs = set(pairs)
        start_tot = time.time()
        self.bins['rdf'] = binwidth * np.array(range(int(np.sqrt(2) *self.box_length / (2 * binwidth))))
        rdfs = {}
        for pair in pairs:
            rdfs[pair] = np.zeros(len(self.bins['rdf']))
        for iatom in range(self.natoms):
            if (iatom % 1000) == 0:
                start = time.time()
            for iatom2 in [i for i in range(iatom + 1, self.natoms)]:
                pair = tuple(sorted([self.labels[iatom], self.labels[iatom2]]))
                if pair in pairs:
                    vect = (self.positions[:,iatom2, :] - self.positions[:, iatom, :])
                    if wrap:
                        vect = (self.positions_wrap[:,iatom2, :] - self.positions_wrap[:, iatom, :])
                    vect = vect - self.box_length * np.rint(vect / self.box_length)
                    dist = np.sqrt(np.sum(np.power(vect, 2), 1))
                    for step in range(len(dist)):
                        if int(dist[step]/binwidth) <  int((np.sqrt(2) *self.box_length / 2)/binwidth):
                           rdfs[pair][int(dist[step]/binwidth)] += 1
            if (iatom % 1000) == 0:
                print "For atom %s finish in:" % iatom, time.time() - start
        vol = np.zeros(len(self.bins['rdf'])+1)
        for i, rr in enumerate(np.append(self.bins['rdf'],self.bins['rdf'][-1]+binwidth)):
            vol[i] = ((4.0/3.0) * np.pi ) * rr**3
            if rr > self.box_length / 2:
                x = self.box_length / (2 * rr)
                vol[i] = vol[i] * ( - 2 + 4.5*x  - 1.5 * x**3)
        for pair in pairs:
            dens = self.count_pair(pair) / self.box_length ** 3
            histo_id = np.zeros(len(self.bins['rdf']))
            for i in range(1, len(self.bins['rdf'])+1):
                histo_id[i-1] = dens*(vol[i] - vol[i-1])
            for i in range(1, len(self.bins['rdf'])):
                rdfs[pair][i] = rdfs[pair][i] / histo_id[i]
            rdfs[pair] =  ( 1 + int(pair[0] == pair[1]) )* rdfs[pair] / self.steps
        self.rdfs = rdfs
        print "Total time is :", time.time() - start_tot

    def calculate_fpt(self, binwidth, freq_tau = 1, length = -1, target_size = 1):
        if not self.is_read:
            self.read_traj_dlpoly()
        start_tot = time.time()
        if length < 0:
            length = self.box_length
        self.bins['fpt'] = binwidth * np.array(range(int(length /  binwidth )))
        local_step = range(0, self.steps, freq_tau)
        self.times['fpt'] = np.array([self.times['saved'][k] for k in local_step ])
        count = np.zeros(len(self.times['fpt']))
        prop = np.zeros( (len(self.bins['fpt']), len(self.times['fpt']) ) )
        for iatom in range(self.natoms):
            for itau, tau in enumerate(local_step):
                for subitau, subtau in enumerate(local_step[:itau]):
                    vect = self.positions[subtau:tau,iatom,:] - self.positions[tau, iatom, :]
                    vect += - self.box_length * np.rint(vect / self.box_length)
                    chi = 1
                    if any( np.sum(np.power(vect, 2),1) < target_size ):
                        chi = 0
                    dist = np.sqrt(np.sum(np.power(vect[0],2)))
                    if int(dist/binwidth) < len(self.bins['fpt']):
                        prop[int(dist/binwidth), itau - subitau ] += chi
                    count[itau- subitau] += 1
        for itau in range(len(self.times['fpt'])):
            if count[itau] > 0:
                prop[:,itau] = prop[:, itau] /  count[itau]
            if np.sum(prop[:,itau]) > 0:
                prop[:,itau] = prop[:, itau] / (binwidth*np.sum(prop[:,itau]))
        self.fpt = prop
        print "Total time is :", time.time() - start_tot


    def calculate_propagator(self, binwidth, freq_tau = 1, length = -1, wrap = False ):
        if not self.is_read:
            self.read_traj_dlpoly()
        start_tot = time.time()
        if length < 0 or wrap:
            length = self.box_length
        txtwrap = ''
        if wrap:
            txtwrap = '_wrap'
        self.bins['prop' + txtwrap] = binwidth * np.array(range(int(length /  binwidth )))
        local_step = range(0, self.steps, freq_tau)
        self.times['prop' + txtwrap] = np.array([self.times['saved'][k] for k in local_step ])
        prop = np.zeros( (len(self.bins['prop' + txtwrap]), len(self.times['prop' + txtwrap])  ) )
        max_ = 0
        for itau, tau in enumerate(local_step):
            count = 0
            #print tau, self.steps, len(range(self.steps - tau))
            l_ = []
            for istep in range(self.steps - tau):
                l_.append((istep+tau, istep))
                vect = self.positions[istep + tau, :, : ] - self.positions[istep, :, :]
                if wrap:
                    vect += - self.box_length * np.rint(vect / self.box_length)
                dist = np.sqrt(np.sum(np.power(vect, 2), 1))
                max_ = max(max(dist), max_)
                for iatom in range(self.natoms):
                    count += 1
                    if int(dist[iatom]/binwidth) < len(self.bins['prop' + txtwrap]):
                        prop[int(dist[iatom]/binwidth), itau] += 1
            prop[:,itau] = prop[:, itau] / ((self.steps - tau))
            prop[:,itau] = prop[:, itau] / (binwidth*np.sum(prop[:,itau]))
        if wrap:
            self.prop_wrap = prop
        else:
            self.prop = prop
        print "Max distance is", max_
        print "Total time is :", time.time() - start_tot

    def determine_flux(self, freq0, length,  binwidth, pairs = []):
        # sum over length
        # frame lab/part
        # check it works
        # add frequency
        # speed-up ?
        if not self.is_read:
            self.read_traj_dlpoly()
        start_tot = time.time()
        self.bins['flux'] = binwidth * np.array(range(int(np.sqrt(2) *self.box_length / (2 * binwidth))))
        j_rads = {}; j_thetas = {}; rdfs = {}
        for pair in pairs:
            j_rads[pair] = np.zeros(len(self.bins['flux']))
            j_thetas[pair] = np.zeros(len(self.bins['flux']))
            rdfs[pair] = np.zeros(len(self.bins['flux']))
        nsteps_mean = self.steps - length
        for iatom in range(self.natoms):
            deltatild = np.roll(self.positions[:,iatom,:], -length, axis=0) - self.positions[:,iatom,:]
            deltatild += - self.box_length * np.rint(deltatild / self.box_length)
            for iatom2 in [i for i in range(iatom + 1, self.natoms)]:
                pair = tuple(sorted([self.labels[iatom], self.labels[iatom2]]))
                if pair in pairs:
                    vel2 = self.velocities[:,iatom2, :]
                    vect = (self.positions[:,iatom2, :] - self.positions[:, iatom, :])
                    vect = vect - self.box_length * np.rint(vect / self.box_length)
                    dist = np.sqrt(np.sum(np.power(vect, 2), 1))
                    longitudinal = np.multiply(vect, deltatild).sum(-1) \
                                       * np.multiply(vect, vel2).sum(-1)
                    transversal = np.multiply(deltatild, vel2).sum(-1)
                    for step in range(nsteps_mean):
                        if int(dist[step] / binwidth) < int((np.sqrt(2) * self.box_length / 2) / binwidth):
                            rdfs[pair][int(dist[step] / binwidth)] += 1
                            j_rads[pair][int(dist[step] / binwidth) ] += longitudinal[step] / dist[step]**4
                            j_thetas[pair][int(dist[step] / binwidth) ] += longitudinal[step] / dist[step]**4 - transversal[step]/dist[step]**2
        vol = np.zeros(len(self.bins['flux']))
        for i, rr in enumerate(self.bins['flux']):
            vol[i] = ((4.0/3.0) * np.pi ) * rr**3
            if rr > self.box_length / 2:
                x = self.box_length / (2 * rr)
                vol[i] = vol[i] * ( - 2 + 4.5*x  - 1.5 * x**3)
        for pair in pairs:
            for i in range(1, self.bins['flux']):
              #  j_rads[pair][i] = j_rads[pair][i]
              #  j_thetas[pair][i] = j_thetas[pair][i]
                rdfs[pair][i] = rdfs[pair][i] / (vol[i] - vol[i-1])
            rdfs[pair] =  ( 1 + int(pair[0] == pair[1]) )* rdfs[pair] * self.box_length ** 3 / (self.count_pair(pair) * nsteps_mean)
            j_rads[pair] = j_rads[pair] / nsteps_mean
            j_thetas[pair] = j_thetas[pair] / nsteps_mean
        self.rdfs = rdfs
        self.j_rads = j_rads
        self.j_thetas = j_thetas
        print "Total time is :", time.time() - start_tot

class CO2(Traj):
    def __init__(self, path):
        super(CO2, self). __init__(path=path)
        self.mass = { 'C' : 12.011, 'O' : 15.99}

    def build_molecules(self):
        self.nmols = self.natoms / 3
        self.molecules = {}
        for imol in range(self.nmols):
            self.molecules[imol] = Molecule(label='CO2', atoms = [imol*3, imol*3+1, imol*3+2])
        self.associate_atoms_mol()


    def determine_map(self, binwidth ):
        start_tot = time.time()
        nbins = int(self.box_length/(2*binwidth))
        self.bins['map'] = np.meshgrid(binwidth*np.array(range(nbins)),binwidth*np.array(range(nbins)))
        self.maps = {}
        for mol in ['CO2']:
            self.maps[mol] = {}
            for atom in set(self.labels.values()):
                self.maps[mol][atom] = np.zeros((nbins, nbins))
        for imol in range(self.nmols):
            if (imol % 500) == 0:
                start = time.time()
            c, o1, o2 = self.molecules[imol].atoms
            u = self.positions[:,o2,:] - self.positions[:,o1,:]
            u += - self.box_length*np.rint(u/self.box_length)
            row_sum =  np.sqrt((u**2).sum(1, keepdims=True))
            u = u / row_sum
            vect1 = self.positions[:,o1,:]
            vect1 += - self.box_length*np.rint(vect1/self.box_length)
            vect2 = self.positions[:,o2,:]
            vect2 += - self.box_length*np.rint(vect1/self.box_length)
            mean = (vect1 + vect2)/2
            latoms = range(self.natoms)
            for iatom in self.molecules[imol].atoms:
                latoms.remove(iatom)
            for iatom in latoms:
                pos = self.positions[:,iatom,:] - mean
                pos += - self.box_length*np.rint(pos/self.box_length)
                dot = np.multiply(u,pos).sum(1)
                row_dot = np.multiply(u,pos).sum(1, keepdims = True)
                #rhos = np.rint(np.sqrt(np.abs((pos**2).sum(1) -dot**2   )) / binwidth).astype(int)
                rhos =  np.rint( np.sqrt( ((pos - row_dot * u)**2).sum(1) ) / binwidth).astype(int)
                zetas = np.rint( dot /binwidth).astype(int)
                for rho, zeta in zip(rhos, zetas):
                    if (rho < nbins) and (0 < zeta + int(np.rint(nbins / 2)) < nbins):
                        self.maps[mol][self.labels[iatom]][rho, zeta+int(np.rint(nbins/2)) ] += 1
            if (imol % 500) == 0:
                print "For atom %s finish in:" % imol, time.time() - start
        vol = (2*np.pi*binwidth**3*np.arange(1, nbins+1, 1))[:, np.newaxis]
        for mol in ['CO2']:
            for atom in set(self.labels.values()):
                self.maps[mol][atom] = self.maps[mol][atom] / vol
        print "Total time is :", time.time() - start_tot


def correlation_matrix(array, length = -1, start = 0):
    new = array.reshape( array.shape[0], array.shape[1]*array.shape[2] )
    covar = np.zeros( (array.shape[1]*array.shape[2], array.shape[1]*array.shape[2]))
    mean = np.zeros(array.shape[1]*array.shape[2])
    count = 0    
    if length == -1:
        length = new.shape[0]
    for step in range(start, start + length):
        line = new[step]
        for i in range(len(line)):
            covar[i,i] += line[i]*line[i]
            mean[i] += line[i]
            for j in range(i+1,len(line)):
                covar[i,j] += line[i]*line[j]
                covar[j,i] += line[j]*line[i]
        count += 1
    covar = covar/count
    mean = mean / count
    for i in range(covar.shape[0]):
        for j in range(covar.shape[1]):
            covar[i,j] += - mean[i]*mean[j]
    return covar 
    
def write_flux_function(path, bins, values):
    with open(path, 'w') as f:
         f.write('R   F\n')
         for r, value in zip(bins, values):
             f.write('%s    %s\n' % (r, value))

def write_map(path, map_):
    with open(path, 'w') as f:
        for line in map_:
            f.write('%s\n' % ' '.join(map(str, line)))

def write_list(path, list_ ):
    with open(path, 'w') as f:
        for line in list_:
            f.write('%s\n' % str(line))