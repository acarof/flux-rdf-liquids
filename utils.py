import re
import numpy as np
import os
import gzip
import scipy.linalg
import time, sys
#import psutil
#from guppy import hpy


class Traj(object):
    
    def __init__(self, path):
        self.path = path
        self.natoms = 0
        self.mass = 28 
        self.is_read = False

    #@profile
    def read_traj_dlpoly(self, freqstep = 1):
        if self.is_read:
            print
            'Traj already read once'
        else:
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
                        print 'Timestep : %s fs' % (self.timestep * 1E3)
                    elif iline == 3:
                        self.box_length = float(line.split()[0])
                        print 'Box length: %s A' % (self.box_length)
                    if (my_step % freqstep) == 0:
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
                print 'Steps: %s' % self.steps
                for step in range(self.steps):
                    if len(self.positions[step]) != self.natoms:
                        print "error for step", step, len(self.positions[step])
                self.positions = np.array(self.positions)
                self.velocities = np.array(self.velocities)
                self.forces= np.array(self.forces)
                self.is_read = True


    def read_traj_dlpoly_old(self):
        if self.is_read:
            print 'Traj already read once'
        else:
            with open('%s/HISTORY' % self.path) as f:
                if self.natoms == 0:
                    self.natoms = int(lines[1].split()[2])
                print 'Number of atoms: %s' % self.natoms
                self.steps = int( (len(lines) - 2) / float(self.natoms*4 + 4) )
                print 'Steps: %s' % self.steps
                self.box_length = float(lines[3].split()[0])
                print 'Box length: %s A' % (self.box_length)
                self.timestep = float(lines[2].split()[5])
                print 'Timestep : %s fs' % (self.timestep*1E3)

                self.labels = {}
                self.forces = []; self.positions = []; self.velocities = []
                for array in [self.forces, self.velocities, self.positions]:
                    print "new array"
                    for istep in range(self.steps):
                        if istep % 1000 == 0:
                            print "istep =", istep
                            print sys.getsizeof(array)
                            print psutil.Process(os.getpid()).memory_info()[0] * 8 / 1E9
                        array.append([])
                        for i in range(self.natoms):
                            array[istep].append( [] )
                for istep in range(self.steps):
                    for iatom in range(self.natoms):
                        iline = 6 + (4 + self.natoms*4)*istep + 4*iatom
                        if self.labels.get(iatom) is None:
                            self.labels[iatom] = lines[iline].split()[0]
                        line = map(float, lines[iline+1].split())
                        self.positions[istep][iatom] += line
                        line = map(float, lines[iline+1].split())
                        self.velocities[istep][iatom] += line
                        line = map(float, lines[iline+1].split())
                        self.forces[istep][iatom] += line
                self.positions = np.array(self.positions)
                self.velocities = np.array(self.velocities)
                self.forces= np.array(self.forces)
                self.is_read = True



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

    def determine_rdf(self, binwidth,  pairs = []):
        if not self.is_read:
            self.read_traj_dlpoly()
        start_tot = time.time()
        nbins = int(np.sqrt(2) *self.box_length / (2 * binwidth))
        print "Use nbins for RDF", nbins
        bins = binwidth * np.array(range(nbins))
        rdfs = {}
        for pair in pairs:
            rdfs[pair] = np.zeros(nbins)
        for iatom in range(self.natoms):
            if (iatom % 1000) == 0:
                start = time.time()
            for iatom2 in [i for i in range(iatom + 1, self.natoms)]:
                pair = tuple(sorted([self.labels[iatom], self.labels[iatom2]]))
                if pair in pairs:
                    vect = (self.positions[:,iatom2, :] - self.positions[:, iatom, :])
                    vect = vect - self.box_length * np.rint(vect / self.box_length)
                    dist = np.sqrt(np.sum(np.power(vect, 2), 1))
                    for step in range(len(dist)):
                        if int(dist[step]/binwidth) <  int((np.sqrt(2) *self.box_length / 2)/binwidth):
                           rdfs[pair][int(dist[step]/binwidth)] += 1
            if (iatom % 1000) == 0:
                print "For atom %s finish in:" % iatom, time.time() - start
        vol = np.zeros(nbins)
        for i, rr in enumerate(bins):
            vol[i] = ((4.0/3.0) * np.pi ) * rr**3
            if rr > self.box_length / 2:
                x = self.box_length / (2 * rr)
                vol[i] = vol[i] * ( - 2 + 4.5*x  - 1.5 * x**3)
        for pair in pairs:
            for i in range(1, nbins):
                rdfs[pair][i] = rdfs[pair][i] / (vol[i] - vol[i-1])
            rdfs[pair] =  ( 1 + int(pair[0] == pair[1]) )* rdfs[pair] * self.box_length ** 3 / (self.count_pair(pair) * self.steps)
        self.bins = bins
        self.rdfs = rdfs
        print "Total time is :", time.time() - start_tot

    def calculate_propagator(self, binwidth, freq_tau = 1, length = -1, wrap = False ):
        if not self.is_read:
            self.read_traj_dlpoly()
        start_tot = time.time()
        if length < 0 or wrap:
            length = self.box_length
        nbins = int(length /  binwidth )
        print "Use nbins for RDF", nbins
        bins = binwidth * np.array(range(nbins))
        taus = range(0, self.steps, freq_tau)
        prop = np.zeros( (nbins, len(taus)  ) )
        max_ = 0
        count = 0
        for itau, tau in enumerate(taus):
            for istep in range(self.steps - tau):
                vect = self.positions[istep + tau, :, : ] - self.positions[istep, :, :]
                if wrap:
                    vect += - self.box_length * np.rint(vect / self.box_length)
                dist = np.sqrt(np.sum(np.power(vect, 2), 1))

                max_ = max(max(dist), max_)
                for iatom in range(self.natoms):
                    count += 1
                    if int(dist[iatom]/binwidth) < nbins:
                        prop[int(dist[iatom]/binwidth), itau] += 1
        if wrap:
            self.bins_prop_wrap = bins
            self.taus_prop_wrap = taus
            self.prop_wrap = prop
        else:
            self.bins_prop = bins
            self.taus_prop = taus
            self.prop = prop
        print "Max distance is", max_
        print count
        print "Total time is :", time.time() - start_tot

    def determine_flux(self, freq0, length,  binwidth, pairs = []):
        if not self.is_read:
            self.read_traj_dlpoly()
        start_tot = time.time()
        nbins = int(np.sqrt(2) *self.box_length / (2 * binwidth))
        print "Use nbins for RDF", nbins
        bins = binwidth * np.array(range(nbins))
        j_rads = {}; j_thetas = {}; rdfs = {}
        for pair in pairs:
            j_rads[pair] = np.zeros(nbins) 
            j_thetas[pair] = np.zeros(nbins) 
            rdfs[pair] = np.zeros(nbins)
        for iatom in range(self.natoms):
            deltatild = np.roll(self.positions[:,iatom,:], length, axis=0) - self.positions[:,iatom,:]
            deltatild += - self.box_length * np.rint(deltatild / self.box_length)
            print deltatild
            for iatom2 in [i for i in range(iatom + 1, self.natoms)]:
                pair = tuple(sorted([self.labels[iatom], self.labels[iatom2]]))
                if pair in pairs:
                    vel2 = self.velocities[:,iatom2, :]
                    vect = (self.positions[:,iatom2, :] - self.positions[:, iatom, :])
                    vect = vect - self.box_length * np.rint(vect / self.box_length)
                    dist = np.sqrt(np.sum(np.power(vect, 2), 1))
                    longitudinal = np.multiply(vect, deltatild).sum(-1) \
                                       * np.multiply(vect, vel2).sum(-1)
                    for step in range(len(dist)):
                        if int(dist[step] / binwidth) < int((np.sqrt(2) * self.box_length / 2) / binwidth):
                            rdfs[pair][int(dist[step] / binwidth)] += 1
                            j_rads[pair][int(dist[step] / binwidth) ] += longitudinal[step]
                    stop
 #                   #try:
 #
 #                       j_thetas[pair][int(dist/binwidth)] += np.dot(deltatild, vect) * np.dot(vel2, vect) / dist**2
 #                       j_thetas[pair][int(dist/binwidth)] += - np.dot(deltatild, vel2)
 #                       rdfs[pair][int(dist/binwidth)] += 1.0
 #                   except:
 #                        pass
 #               deltatild += - self.velocities[step][iatom]
        for pair in pairs:
            for i in range(nbins):
                if bins[i] != 0:
                    j_rads[pair][i] = j_rads[pair][i] /rdfs[pair][i]
                    j_thetas[pair][i] = j_thetas[pair][i] / rdfs[pair][i]
                    rdfs[pair][i] = rdfs[pair][i]  / bins[i]**2
                    #j_rad[i] = j_rad[i] / bins[i]**2
            rdfs[pair] = rdfs[pair]*self.box_length**3 / (4*np.pi*binwidth*count[pair])
        self.bins = bins
        self.rdfs = rdfs
        self.j_rads = j_rads
        self.j_thetas = j_thetas

    
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
