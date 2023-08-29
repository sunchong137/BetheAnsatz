# Copyright 2016-2023 Chong Sun (sunchong137@gmail.com)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import numpy as np
from matplotlib import pyplot as plt
import sys
import os
sys.path.append("./gs_knizia/")

class BetheAnsatzFT(object):
    # half-filling, Sz = 0
    def __init__(self, U, T, ngrid,gaussian=False):
        self.U = U/4.
        self.T = T
        self.ngrid = ngrid
        self.Q = np.pi
        self.B = 10
        self.kgrid = np.zeros(ngrid)
        self.lgrid = np.zeros(ngrid)
        self.kwts  = np.ones(ngrid) #weights of the k point for integration
        self.lwts  = np.ones(ngrid)
        #self.dk = 2.*self.Q / ngrid
        #self.dl = 2.*self.B / ngrid
        #some functions used in integration functions
        self.kappa0 = np.zeros(ngrid)
        self.kappa = np.zeros(ngrid)
        self.eps1m = np.zeros(ngrid)
        self.sig0 = np.zeros(ngrid)
        self.rho0 = np.zeros(ngrid)
        self.e0 = 0.
        self.grandpot = 0.
        self.gaussian_grids = gaussian
    
    def generate_grids(self):
        # symmetrized partition
        if self.gaussian_grids:
            self.kgrid, self.kwts = _get_scaled_legendre_roots(-self.Q,self.Q,self.ngrid)
            self.lgrid, self.lwts = _get_scaled_legendre_roots(-self.B,self.B,self.ngrid)
        else:
            self.kgrid = np.linspace(-self.Q, self.Q, self.ngrid, endpoint=False)
            self.lgrid = np.linspace(-self.B, self.B, self.ngrid, endpoint=False)
            dk = 2.*self.Q / self.ngrid
            dl = 2.*self.B / self.ngrid
            self.kgrid += dk/2.
            self.lgrid += dl/2.
            self.kwts  *= dk
            self.lwts  *= dl
        
    def get_sig0(self):
        sig0 = np.zeros(self.ngrid)
        i = 0
        for l in self.lgrid:
            sig0[i] = np.sum(self.sfunc(l, np.sin(self.kgrid))*self.kwts)/(2.*np.pi)
            i += 1
        self.sig0 = sig0
        return sig0

    def get_rho0(self):
        sig0 = self.sig0
        if(np.linalg.norm(sig0)<1e-9):
            sig0 = self.get_sig0()
        rho0 = np.zeros(self.ngrid)
        i = 0
        for k in self.kgrid:
            rho0[i] = (self.U*np.cos(k)/np.pi) * np.dot(sig0*self.lwts, 1./((self.lgrid-np.sin(k))**2.+self.U**2.))
            
            i += 1
        rho0 += 1./(2.*np.pi)
        self.rho0 = rho0
        return rho0

    def get_e0(self):
        rho0 = self.rho0
        if(np.linalg.norm(rho0) < 1e-9):
            rho0 = self.get_rho0()
        e0 = -2.*np.dot(np.cos(self.kgrid), rho0*self.kwts)
        self.e0 = e0
        return e0
        
    def get_grandpot(self):
        part1 = np.dot(self.rho0*self.kwts, self.Gfunc(self.kappa))
        part2 = np.dot(self.sig0*self.lwts, self.eps1m-self.T*np.log(1./4.))
        Omega = self.e0 - self.U*2. - part1 - part2
        self.grandpot = Omega
        return Omega
    
    def sfunc(self, l, k):
        self.s = 1./(4.*self.U*np.cosh(np.pi*(l-k)/(2*self.U)))
        return self.s
    def kappa0func(self):
        kappa0 = np.zeros(self.ngrid)
        i = 0
        for k in self.kgrid:
            kappa0[i] = -2.*np.cos(k) - 4.*np.dot(self.lwts*self.sfunc(self.lgrid,np.sin(k)),np.sqrt(1.-(self.lgrid-1.j*self.U)**2).real)
            i+=1
        self.kappa0 = kappa0
        return kappa0
    def Gfunc(self, y):
        return self.T * np.log(1.+np.exp(y/self.T))
    def sconvfunc(self, fl):
        convf = np.zeros(self.ngrid)
        i = 0
        for l in self.lgrid:
            convf[i] = np.dot(self.sfunc(l, self.lgrid), fl*self.lwts)
            i += 1
        return convf
    # functions to be used in self-consistency
    def kappafunc(self, eps1p, eps1m):
        kappa = np.zeros(self.ngrid)
        self.kappa0 = self.kappa0func()
        i = 0
        for k in self.kgrid:
            kappa[i] = np.dot(self.sfunc(self.lgrid, np.sin(k))*self.lwts, eps1p-eps1m)
            i += 1
        self.kappa = kappa + self.kappa0
        return self.kappa
    
    def epsnfunc(self, n, En):
        bn = 1./((n+1.)**2)
        return self.T * np.log(bn + (1-bn)*np.exp(En/self.T))

    def E1func(self, eps2, p=1.):
        # p : parity
        E1 = np.zeros(self.ngrid)
        i = 0
        for l in self.lgrid:
            E1[i] = -np.dot(np.cos(self.kgrid)*self.sfunc(l,np.sin(self.kgrid)), self.Gfunc(p*self.kappa)*self.kwts)
            i += 1

        E1 += self.sconvfunc(eps2)
        return E1

    def Enfunc(self, epsm, epsp):
        return self.sconvfunc(epsm + epsp)
   
    def solve_epsilon_scf(self, nstep=150, tol = 1e-10):
    
        print "Calculating epsilon and kappa..."
        nmax  = 20
        ngrid = self.ngrid
        Ep   = np.zeros((nmax+1, ngrid))
        Em   = np.zeros((nmax+1, ngrid))
        epsp = np.zeros((nmax+1, ngrid))
        epsm = np.zeros((nmax+1, ngrid))
        # initial guess: eps1p = eps1m, high orders are zero
        kappa = self.kappa.copy()
        for itr in range(nstep):
            self.kappa = self.kappafunc(epsp[0], epsm[0])
            Ep[0] = self.E1func(epsp[1], 1.)
            Em[0] = self.E1func(epsm[1], -1.)
            epsp[0] = self.epsnfunc(1, Ep[0])
            epsm[0] = self.epsnfunc(1, Em[0])
            for i in range(1, nmax):
                Ep[i] = self.sconvfunc(epsp[i-1]+epsp[i+1])
                Em[i] = self.sconvfunc(epsm[i-1]+epsm[i+1])
                epsp[i] = self.epsnfunc(i+1, Ep[i])
                epsm[i] = self.epsnfunc(i+1, Em[i])
            diff = np.linalg.norm(self.kappa - kappa)
            #print diff
            if diff < tol:
                #print "The convergence is achieved after %d loops!"%(itr+1), "Final diff:", diff
                break
            else:
                kappa = self.kappa.copy()
        if(itr == nstep-1):
            print "The convergence is NOT acheived after %d loops!"%nstep, "Final diff:", diff
        self.eps1m = np.copy(epsm[0])

    def solve_grandpot(self):
        self.generate_grids()
        self.solve_epsilon_scf()
        self.e0 = self.get_e0()
        grandpot = self.get_grandpot()
        return grandpot

def _get_scaled_legendre_roots(wl, wh, nw):
    """
    Scale nw Legendre roots, which lie in the
    interval [-1, 1], so that they lie in [wl, wh]

    Returns:
        freqs : 1D ndarray
        wts : 1D ndarray
    """
    freqs, wts = np.polynomial.legendre.leggauss(nw)
    freqs += 1
    freqs *= (wh - wl) / 2.
    freqs += wl
    wts *= (wh - wl) / 2.
    return freqs, wts

def solve_energy_curve(U, Tgrid, outdir='./data', dT=0.006, ngrid=100,savefile=False, gaussian=False, save_beta=True, **kwargs):
    
    entropy = []
    grandpot = []
    Tmin = Tgrid[0]
    Tmax = Tgrid[-1]

    FiniteTgrid = Tgrid
    
    if Tmin <1e-5:
        import bethe_ansatz as gsba
        
        e,_ = gsba.CalcBetheEnergy_UandN(U, 1.0)
        grandpot.append([0, e-U/2.])
        entropy.append([0, 0])
        FiniteTgrid = Tgrid[1:]
        save_beta = False


    for T in FiniteTgrid:
        obj = BetheAnsatzFT(U, T, ngrid,gaussian=gaussian)
        g = obj.solve_grandpot()
        objp = BetheAnsatzFT(U, T+dT, ngrid,gaussian=gaussian)
        objm = BetheAnsatzFT(U, T-dT, ngrid,gaussian=gaussian)
        s = -(objp.solve_grandpot() - objm.solve_grandpot())/(2*dT)
        if save_beta:
            grandpot.append([1./T, g])
            entropy.append([1./T, s])
        else:
            grandpot.append([T, g])
            entropy.append([T, s])
            
        e = g + U/2. + T*s
        print "T: %0.4f       GrandPot: %0.6f"%(T, g)
        print "T: %0.4f       Entropy:  %0.6f"%(T, s)
        print "T: %0.4f       Energy:  %0.6f"%(T, e)
    entropy = np.asarray(entropy)
    grandpot = np.asarray(grandpot)
    energy = grandpot.copy()
    energy[:,1] += U/2. + Tgrid*entropy[:,1]
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    edir = outdir+"/energy/"
    gdir = outdir+"/grandpot/"
    sdir = outdir+"/entropy/"

    for d in [edir, gdir, sdir]:
        if not os.path.exists(d):
            os.mkdir(d)
        
    if(savefile):
        np.savetxt(edir+"/energy_BA_U%d.txt"%U, energy)
        np.savetxt(gdir+"/grandpot_BA_U%d.txt"%U, grandpot)
        np.savetxt(sdir+"/entropy_BA_U%d.txt"%U, entropy)
    else:
        print energy

    return energy

if __name__ == "__main__":
    U = float(sys.argv[1])
    ngrid = 200
    write2file = True
    gaussian = True
    outdir = "data/"
    beta = np.array([10.])#np.linspace(0.1,10,100,endpoint=True)
    beta = np.arange(10,20,1)
    Tgrid = 1./beta
    #Tgrid = np.linspace(0.00,2.0,41,endpoint=True)
    #Tgrid = np.linspace(1.05,2.0,20, endpoint=True)
    solve_energy_curve(U,Tgrid,ngrid=ngrid,savefile=write2file,gaussian=gaussian)


