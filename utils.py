import re
import numpy as np
import os
import gzip
import scipy.linalg
import time, sys
#import psutil
#from guppy import hpy

def calculate_fft(time, serie):
    serie = np.append(serie, np.zeros(serie.shape))
    fft = np.fft.rfft(serie)
    length = len(fft)
    xaxis = np.arange(length) * (np.pi / time[-1])
    return xaxis, fft

class Traj(object):
    
    def __init__(self, path):
        self.path = path
        self.natoms = 0
        self.mass = 28 
        self.is_read = False
        self.times = {}
        self.bins = {}

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
                #self.times['printed'] = printstep * np.arange(self.printed_steps)
                #self.times['MD'] = self.timestep * np.arange(self.total_steps)
                self.shift_com()
                self.unwrap()
                print "Trajectory is read"

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

    def determine_rdf_forces(self, binwidth, pairs = [],):
        if not self.is_read:
            self.read_traj_dlpoly()
        start_tot = time.time()
        self.bins['rdf_forces'] = binwidth * np.array(range(int(np.sqrt(2) *self.box_length / (2 * binwidth))))
        self.rdf_forces = {}
        for pair in pairs:
            self.rdf_forces [pair] = np.zeros(len(self.bins['rdf_forces']))
        for iatom in range(self.natoms):
            if (iatom % 1000) == 0:
                start = time.time()
            for iatom2 in [i for i in range(iatom + 1, self.natoms)]:
                pair = tuple(sorted([self.labels[iatom], self.labels[iatom2]]))
                if pair in pairs:
                    vect = (self.positions[:, iatom2, :] - self.positions[:, iatom, :])
                    vect = vect - self.box_length * np.rint(vect / self.box_length)
                    dist = np.sqrt(np.sum(np.power(vect, 2), 1))
                    diff_forces = (self.forces[:,iatom, :] - self.forces[:,iatom2,:])
                    dot = (diff_forces * vect).sum(1)
                    toadd = dot / dist**3
                    for step in range(len(dist)):
                        n =int(dist[step] / binwidth)
                        if  n < int((np.sqrt(2) * self.box_length / 2) / binwidth):
                            for k in range(n+1):
                                self.rdf_forces[pair][k] += toadd[step]
            if (iatom % 1000) == 0:
                print "For atom %s finish in:" % iatom, time.time() - start
        vol = np.zeros(len(self.bins['rdf_forces'])+1)
        for i, rr in enumerate(np.append(self.bins['rdf_forces'],self.bins['rdf_forces'][-1]+binwidth)):
            vol[i] = ((4.0/3.0) * np.pi ) * rr**3
            if rr > self.box_length / 2:
                x = self.box_length / (2 * rr)
                vol[i] = vol[i] * ( - 2 + 4.5*x  - 1.5 * x**3)
        for pair in pairs:
            dens = self.count_pair(pair) / self.box_length ** 3
            histo_id = np.zeros(len(self.bins['rdf_forces']))
            for i in range(1, len(self.bins['rdf_forces'])+1):
                histo_id[i-1] = dens*(vol[i] - vol[i-1])
            for i in range(1, len(self.bins['rdf_forces'])):
                self.rdf_forces[pair][i] = self.rdf_forces[pair][i] / histo_id[i]
            self.rdf_forces[pair] =  ( 1 + int(pair[0] == pair[1]) )* self.rdf_forces[pair] / self.steps
        print "Total time is :", time.time() - start_tot


    def determine_rdf(self, binwidth,  pairs = [], wrap = False):
        if not self.is_read:
            self.read_traj_dlpoly()
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
        self.times['fpt'] = range(0, self.steps, freq_tau)
        count = np.zeros(len(self.times['fpt']))
        prop = np.zeros( (len(self.bins['fpt']), len(self.times['fpt']) ) )
        for iatom in range(self.natoms):
            for itau, tau in enumerate(self.times['fpt']):
                for subitau, subtau in enumerate(self.times['fpt'][:itau]):
                    vect = self.positions[subtau:tau,iatom,:] - self.positions[tau, iatom, :]
                    vect += - self.box_length * np.rint(vect / self.box_length)
                    chi = 1
                    if any( np.sum(np.power(vect, 2),1) < target_size ):
                        chi = 0
                    dist = np.sum(np.power(vect[0],2))
                    if int(dist/binwidth) < len(self.bins['fpt']):
                        prop[int(dist/binwidth), itau - subitau ] += chi
                    count[itau- subitau] += 1
        for itau in range(len(self.times['fpt'])):
            prop[:,itau] = prop[:, itau] /  count[itau]
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
        self.times['prop' + txtwrap] = range(0, self.steps, freq_tau)
        prop = np.zeros( (len(self.bins['prop' + txtwrap]), len(self.times['prop' + txtwrap])  ) )
        max_ = 0
        for itau, tau in enumerate(self.times['prop' + txtwrap]):
            count = 0
            #print tau, self.steps, len(range(self.steps - tau))
            l_ = []
            for istep in range(self.steps - tau):
                l_.append((istep+tau, istep))
                vect = self.positions[istep + tau, :, : ] - self.positions[istep, :, :]
                if wrap:
                    vect += - self.box_length * np.rint(vect / self.box_length)
                dist = np.sqrt(np.sum(np.power(vect, 2), 1))
                if np.isnan(dist).any():
                    print "problem in dist"
                    stop
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

    def build_molecules(self):
        self.nmols = self.natoms / 3
        self.mol_indexes = range(self.nmols)
        for imol in range(self.nmols):
            self.mol_indexes[imol] = [imol*3, imol*3+1, imol*3+2]

    def determine_map(self, binwidth ):
        nbins = int(self.box_length/(2*binwidth))
        self.bins['map'] = np.meshgrid(binwidth*np.array(range(nbins)),binwidth*np.array(range(nbins)))
        self.maps = {}
        for mol in ['CO2']:
            self.maps[mol] = {}
            for atom in set(self.labels.values()):
                self.maps[mol][atom] = np.zeros((len(self.bins['map']),len(self.bins['map'])))
        for imol in range(self.nmols):
            c, o1, o2 = self.mol_indexes[imol]
            u = self.positions[:,o2,:] - self.positions[:,o1,:]
            u += - self.box_length*np.rint(u/self.box_length)
            u = u / np.linalg.norm(u)
            vect1 = self.positions[:,o1,:]
            vect1 += - self.box_length*np.rint(vect1/self.box_length)
            vect2 = self.positions[:,o2,:]
            vect2 += - self.box_length*np.rint(vect1/self.box_length)
            mean = (vect1 + vect2)/2
            for iatom in range(self.natoms):
                pos = self.positions[:,iatom,:] - mean
                pos += - self.box_length*np.rint(pos/self.box_length)
                rhos = np.zeros(self.steps)
                for step in range(self.steps):
                    rhos[step] = np.rint( np.linalg.norm( np.cross(u[step], pos[step]) )/binwidth )
                rhos = rhos.astype(int)
                zetas = np.rint(np.multiply(u,pos).sum(1) /binwidth).astype(int)
                for rho, zeta in zip(rhos, zetas):
                    if (rho < len(self.bins['map'])) and (0 < zeta + int(np.rint(len(self.bins['map']) / 2)) < len(self.bins['map'])):
                        self.maps[mol][self.labels[iatom]][rho, zeta+int(np.rint(len(self.bins['map'])/2)) ] += 1



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
