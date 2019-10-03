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

                step = 0
                self.labels = {}
                self.forces = []; self.positions = []; self.velocities = []
                print "I use the correct version 2"
                for array in [self.forces, self.velocities, self.positions]:
                #for array in [self.forces, self.velocities]:
                    print "new array"
                    for istep in range(self.steps):
                        if istep % 1000 == 0:
                            print "istep =", istep
                            print sys.getsizeof(array)
                            print psutil.Process(os.getpid()).memory_info()[0] * 8 / 1E9
                        array.append([])
                        for i in range(self.natoms):
                            array[istep].append( [] )
                print "I create the array"
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
    
    def determine_flux(self, freq0, binwidth, length, pairs = []):
        if not self.is_read:
            self.read_traj_dlpoly()
        nbins = int(self.box_length/(2*binwidth))
        print "Use nbins", nbins
        bins = binwidth * np.array(range(nbins))
        j_rads = {}; j_thetas = {}; rdfs = {}; count = {}
        for pair in pairs:
            j_rads[pair] = np.zeros(nbins) 
            j_thetas[pair] = np.zeros(nbins) 
            rdfs[pair] = np.zeros(nbins)
            count[pair] = 0
        for iatom in range(self.natoms):
            deltatild = np.zeros(3) 
            deltatild  = (self.positions[length][iatom] - self.positions[0][iatom])
            deltatild += - self.box_length * np.floor( 2*deltatild/self.box_length) / 2
            for step in range(length):
                for iatom2 in  [i for i in range(self.natoms) if i !=iatom]:
                    pair = tuple(sorted([self.labels[iatom], self.labels[iatom2]]))
                    if pair in pairs:
                        vel2 = self.velocities[step][iatom2]
                        vect = (self.positions[step][iatom2] - self.positions[step][iatom])
                        vect = np.array([x - self.box_length * np.rint(x / self.box_length) for x in vect])
                        dist = np.sqrt(np.sum(np.power(vect, 2)))
                        count[pair] += 1
                        try:
                            j_rads[pair][int(dist/binwidth)] += np.dot(deltatild, vect) * np.dot(vel2, vect)   / dist**2  
                            j_thetas[pair][int(dist/binwidth)] += np.dot(deltatild, vect) * np.dot(vel2, vect) / dist**2
                            j_thetas[pair][int(dist/binwidth)] += - np.dot(deltatild, vel2)
                            rdfs[pair][int(dist/binwidth)] += 1.0

                        except:
                             pass
                        
                    #else:
                    #    print 'Pair is not in pairs:', pair
                deltatild += - self.velocities[step][iatom]
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