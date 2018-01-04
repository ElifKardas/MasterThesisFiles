# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 08:45:16 2017

@author: Jurg Spaak
"""

from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# fitting the monocultures
def fitting_mono_treat(N0,Nt):
    def model(t, lamb, alpha):
        def iterator(N,i):
            if i>1:
                return iterator(lamb*N/(1+alpha*N),i-1)
            else:
                return lamb*N/(1+alpha*N)
        return np.array([iterator(N0,s) for s in t])
    time = np.tile(np.arange(1,len(Nt)//3+1),3)
    lamb_pre = np.mean(Nt[:3])/N0
    alpha_pre = -(1-lamb_pre)/np.amax(Nt)
    popts = curve_fit(model, time, Nt, [lamb_pre,alpha_pre], 
                         bounds = [[1,0], [100,1]])[0]
    '''plt.plot(time,Nt,'go')
    plt.plot(0,N0,'o')

    
    time2 = [1,2,3]
    plt.plot([0]+list(time2),[N0]+list(model(time2, popts[0], popts[1])),)
    plt.show()
    print(popts)'''
    return popts
    
def fitting_mono_rep(N0,Nt):
    def model(t, lamb, alpha):
        def iterator(N,i):
            if i>1:
                return iterator(lamb*N/(1+alpha*N),i-1)
            else:
                return lamb*N/(1+alpha*N)
        return np.array([iterator(N0,s) for s in t])
    time = np.array([1,2,3])
    lamb_pre = np.mean(Nt[0])/N0
    alpha_pre = -(1-lamb_pre)/np.amax(Nt)
    if alpha_pre == 100:
        print(alpha_pre)
        alpha_pre = 50
    popts = curve_fit(model, time, Nt, [lamb_pre,alpha_pre], 
                          bounds = [[1,0], [100,1]])[0]
    """plt.plot(time,Nt,'go')
    plt.plot(0,N0,'o')
    
    
    time2 = [1,2,3]
    plt.plot([0]+list(time2),[N0]+list(model(time2, popts[0], popts[1])),)
    plt.show()
    print(popts)"""
    return popts
    
densities_mono = np.genfromtxt("densities_mono.csv", delimiter = ",", skip_header = 1)

popts_treat_mono = np.empty((12,4))
popts_rep_mono = np.empty((36,5))

for i,conc in list(enumerate([0,1,10,100])):
    dens_1 = densities_mono[densities_mono[:,0] == conc]
    for j,spec in list(enumerate([1,2,3])):
        dens = dens_1[dens_1[:,1]==spec]
        N0 = np.average(dens[:,2])
        Nt = dens[:,3:].reshape(-1)
        popts_treat_mono[j+i*3] = conc, spec, *fitting_mono_treat(N0,Nt)
        for k in range(3):
            N0 = dens[k,2]
            Nt = dens[k,3:]
            popts_rep_mono[k+3*(j+3*i)] = conc, spec, k+1, *fitting_mono_rep(N0,Nt)
        
np.savetxt("alpha, lambda,treatments.csv", popts_treat_mono,
           header ="Conc;Species;Lambda; alpha",delimiter = ";")
np.savetxt("alpha, lambda,replica.csv", popts_rep_mono,
           header ="Conc;Species;Replica;Lambda; alpha", delimiter = ";")

###############################################################################
# fitting competition

def fitting_comp(N0,Nt, lambs, alphamono):
    def model(t, alpha12, alpha21):
        alpha = np.array([[alphamono[0], alpha12],[alpha21,alphamono[1]]])
        def iterator(N,i):
            if i>1:
                return iterator(lambs*N/(1+np.sum(alpha*N, axis = 1)),i-1)
            else:
                return lambs*N/(1+np.sum(alpha*N, axis = 1))
        return np.array([iterator(N0,s) for s in t]).reshape(-1)
    time = np.array([1,2,3])
    popts = curve_fit(model, time, Nt.reshape(-1), [0,0], 
                          bounds = [0, alphamono])[0]
    '''plt.plot(time, Nt, 'go')
    plt.plot([0,0], N0, 'ro')
    Nt_fit = model(time, *popts).reshape(Nt.shape)
    plt.plot([0,1,2,3],[N0[0]]+list(Nt_fit[:,0]))
    plt.plot([0,1,2,3],[N0[1]]+list(Nt_fit[:,1]))
    plt.semilogy()
    plt.show()'''
    return popts
    
densities_comp = np.genfromtxt("densities_comp.csv", delimiter = ",", skip_header = 1)
popts_treat_comp = np.empty((12,4))        
popts_rep_comp = np.empty((36,5))
asoc_spec = np.array([[1,2], [1,3], [2,3]])
for i,conc in list(enumerate([0,1,10,100])):
    dens_1 = densities_comp[densities_comp[:,1] == conc]
    for j,asoc in list(enumerate([1,2,3])):
        dens = dens_1[dens_1[:,2]==asoc]
        for k in range(3): # different replica
            #print(k+3*(j+3*i), i,j,k)
            N0 = dens[k,2+asoc_spec[j]]
            times = (3*np.arange(1,4)+k)[:,np.newaxis]
            lamb = popts_treat_mono[i*3+asoc_spec[j]-1, 2]
            alpha_mono = popts_treat_mono[i*3+asoc_spec[j]-1, 3]
            Nt = dens[times, 2+asoc_spec[asoc-1]]
            #print(Nt)
            popts_rep_comp[k+3*(j+3*i)] = conc, asoc, k+1, *fitting_comp(N0,Nt, lamb, alpha_mono)
        
        
np.savetxt("alpha, competition.csv", popts_rep_comp,
           header ="Conc;Association;Replica;alpha12;alpha21", delimiter = ";")        
        
        