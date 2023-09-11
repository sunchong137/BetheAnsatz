# This program is released under a "I hope no one has to deal with this
# ever again"-license. You may do with it as you wish, including
# distributing it in source- or binary form, as long as you do not blame
# me if anything goes wrong.
# 
# If you find the program useful, I would appreciate a citation to
# 
#    Knizia, Chan - Density Matrix Embedding: A Simple Alternative
#       to Dynamical Mean-Field Theory
#    http://dx.doi.org/10.1103/PhysRevLett.109.186404
# 
# in published work using the program. That is the paper for which it
# was written. In any case you should cite
# 
#    Shiba - Magnetic Susceptibility at Zero Temperature for the
#       One-Dimensional Hubbard Model
#    http://dx.doi.org/10.1103/PhysRevB.6.930
# 
# which describes the algorithm which was employed in this program.
# -- Gerald Knizia, 2014-06-03

from __future__ import print_function
import numpy as np
import scipy.optimize as optimize
from subprocess import getoutput
from math import pi

CALL_DIR = './'

def CalcBetheEnergy_UandQ(U, Q_):
   """invoke the C++ program for given parameters U, Q. Return energy for
   those parameters. Q is given in units of radians and must be in range
   of (0...2) (corresponding to just-above zero to just below-full filling)
   """
   if Q_ > 1:
      # "code" for N > 1. Transform to the equivalent system with N < 1.
      # (iterative algorithm works only for N <= 1)
      E, N = CalcBetheEnergy_UandQ(U, 2. - Q_)
      return 2*((1.-N/2)-.5)*U + E, 2. - N

   # clamp to valid range. Q=0 is empty, Q=1.0 is half-filling.
   Q = np.clip(Q_, 0.0, 1.0)
   Text = getoutput(CALL_DIR+"bethe_ansatz %.15f %.15f" % (U, Q))
   Lines = Text.splitlines()[-1]
   ls = Lines.split()
   E = float(ls[1])
   N = float(ls[3])
   return E, N

def FindQForN(U, TargetN):
   """Find the Q(U,n) parameter which will result in a filling fraction
   of <n> = TargetN for the given U. Note: This is done by numerical root
   search."""
   assert(0. <= TargetN <= 2.)
   if TargetN > 1:
      return 2. - FindQForN(U, 2.-TargetN)
   if TargetN == 1.0:
      return 1.  # no need to search for that.
   Tol = 1e-10
   def Err(Q):
      E, N = CalcBetheEnergy_UandQ(U, Q)
      return (N - TargetN)**2
   try:
      Q = optimize.brent(Err, brack = (0.0, 0.95, 1.0), tol=Tol, maxiter=1000)
   except ValueError:
      Q = optimize.fminbound(Err, 0.0, 1.0, xtol=Tol, maxfun=1000)
   E, N = CalcBetheEnergy_UandQ(U, Q)
   assert(abs(N - TargetN) <= 1e-5)
   return Q

def CalcBetheEnergy_UandN(U, N):
   """as CalcBetheEnergy_UandQ(U,Q), but use U and <n> as input. This requires first
   solving for Q such that <n> = n(U,Q) and is consequently much slower
   than CalcBetheEnergy_UandQ."""
   return CalcBetheEnergy_UandQ(U, FindQForN(U,N))

def GetLeftMuAtHalfFilling(U):
   """Compute derivative d/d[n] E(U,n) for limit n -> 1.0 - 0
   (i.e., on left side of n=1)"""
   Ql = FindQForN(U, 1.0 - 1e-7)
   Em, Nm = CalcBetheEnergy_UandQ(U, Ql)
   E0, N0 = CalcBetheEnergy_UandQ(U, 1.0)
   return (E0-Em)/(N0-Nm)

def FindGaps():
   """Compute only the one-particle gaps at half filling."""
   TargetUs = [4.]
   for U in TargetUs:
      MuLeft = GetLeftMuAtHalfFilling(U)
      MuHalf = U/2
      Gap = 2 * (MuHalf - MuLeft)
      print("U = {:.4f}  MuLeft = {:.6f}  MuHalf = {:.6f}  Gap = {:.6f}".format(U, MuLeft, MuHalf, Gap))

def MakeSpinCorrelationCurves():
   """Make table of spin correlation functions."""
   ResultLines = []

   TargetUs = []
   TargetUs = [4.,2.,1.,0.5,0.25,0.125]
   #Ns = [.25,1./3.,.5,2./3.,1.0]
   #Ns = [1./3.,2./3., .25,.5,.75, 1.]
   Ns = [1./3.,2./3.,1.]

   def SzSz_for_NaNb(N, NaNb):
      """find <Sz(i) Sz(i)> correlation function for given <nA(i) nB(i)> correlation
      function and filling <n>."""
      SzSz = 3*(N - 2*NaNb)/4
      # ^- PRB 37 7541 (1988) eq 33 (in vector form in that paper)
      #    cites Shiba PRB 6 930, eq 2.12

      # FIXME: hm...  Shibas' correlation functions for <\vec S^2> go from 3/8 to 3/4;
      # mine for <Sz.Sz> from .5 to 1.0. The numbers appear to agree perfectly,
      # but I am not sure why I should divide by (3/4)!
      # That seems to be some normalization factor for S^2. SzSz = N - 2*NaNb makes
      # sense to me (w/o the 3/4.). I just accept it here...
      SzSz /= (3./4.)
      return SzSz

   if (tuple(Ns) == (1./3., 2./3., 1.)):
      ResultLines.append("       0.000   0.333333   0.027778 {:10.6f} {:10.6f}  * note: manual. <na.nb> = 1/(9.*4)".format(1/(9.*4), SzSz_for_NaNb(1/3., 1/(9.*4))))
      ResultLines.append("       0.000   0.666667   0.111111 {:10.6f} {:10.6f}  * note: manual. <na.nb> = 1/9.".format(1/9., SzSz_for_NaNb(2/3.,1/9.)))
      ResultLines.append("       0.000   1.000000   0.250000 {:10.6f} {:10.6f}  * note: manual. <na.nb> = 1/4.".format(1/4., SzSz_for_NaNb(1., 1/4.)))
      # ^- exact solutions for some boundary cases without Us.

   for N in Ns:
      for U in TargetUs:
         if U == 0.:
            continue

         # <nalpha nbeta> = d/dU E(U,n)
         h = 0.001
         Ep,Np = CalcBetheEnergy_UandN(U+h, N)
         Em,Nm = CalcBetheEnergy_UandN(U-h, N)
         E0,N0 = CalcBetheEnergy_UandN(U,N)
         NaNb = (Ep-Em)/(2*h)
         SzSz = SzSz_for_NaNb(N, NaNb)
         print("U = {:10.3f}  <n> = {:10.4f}  <na.nb> = {:10.4f}  <Sz.Sz> = {:10.4f}   {:.2e} {:.6f} {:.6f}".format(U, N, NaNb, SzSz, Ep-Em, Np, Nm))
         ResultLines.append("  {:10.3f} {:10.6f} {:10.6f} {:10.6f} {:10.6f}".format(U, N, E0, NaNb, SzSz))
   print()
   print("  {:>10} {:>10} {:>10} {:>10} {:>10}".format("U  ", "<n>  ", "E  ", "<na.nb>", "<Sz.Sz>"))
   print("   " + 5 * "-----------")
   ResultLines.sort()
   print("\n".join(ResultLines))

def MakeEnergyCurvesAtHalfFilling():
   print("     U         Filling       E/Site         Q/pi")
   print("---------------------------------------------------")
   for U in np.arange(0.25,128.01,0.25):
      for Q in [1.]:
         E, N = CalcBetheEnergy_UandQ(U, Q)
         print("{:8.3f} {:13.6f} {:13.6f} {:13.6f}".format(U, .5*N, E, Q))
      #print

def MakeCurvesForQ():
   """Compute E(U,n) curves, with N given implicitly as a function of Q. (i.e.,
   the curves will be a representation of the full E(U,n) range, but they will have
   different N coordinates for different U).
   Note that the Energy is supposed to have a spike and mu a gap """
   Cap = "   {:^7s} {:^13s} {:^13s} {:^13s} {:^13s}".format("U", "Filling", "E/Site", "Chem.Pot", "Q/pi")
   print(Cap + "\n" + len(Cap)*"-")
   #for U in np.arange(0.25,32.01,0.25):
   for U in [4.]:
      for Q in np.linspace(0.0, 2.0, 161):
      #for Q in np.arange(0.98,1.02+1e-8,0.001):
         dQ = 1e-3
         if Q <= dQ or (Q >= 1.0 - dQ and Q <= 1.0 + dQ):
            continue # can't make numerical derivatives there.
         Ep, Np = CalcBetheEnergy_UandQ(U, Q+dQ)
         Ec, Nc = CalcBetheEnergy_UandQ(U, Q)
         Em, Nm = CalcBetheEnergy_UandQ(U, Q-dQ)
         Mu = (Ep-Em)/(Np-Nm)
         print("{:8.3f} {:13.6f} {:13.6f} {:13.6f} {:13.6f}".format(U, .5*Nc, Ec, Mu, Q))
      print()

def MakeRefDataAtHalfFilling():
   """Compute various information relevant at and around half filling."""
   print("     U        Filling       E/Site        HL-Gap        Chem.Pot      <na.nb>")
   print("------------------------------------------------------------------------------")
   L2 = []
   for U in np.arange(1.,16.+1e-8,1.):
      Q = 1.
      dQ = 1e-3
      Eq, Nq = CalcBetheEnergy_UandQ(U, Q+2*dQ)
      Ep, Np = CalcBetheEnergy_UandQ(U, Q+dQ)
      Ec, Nc = CalcBetheEnergy_UandQ(U, Q)
      Em, Nm = CalcBetheEnergy_UandQ(U, Q-dQ)
      En, Nn = CalcBetheEnergy_UandQ(U, Q-2*dQ)

      E = Ec
      MuN = (Ep-Eq)/(Np-Nq)
      MuP = (En-Em)/(Nn-Nm)
      Gap = MuN - MuP
      Mu = (Ep-Em)/(Np-Nm)

      h = 0.001
      Fp,Mp = CalcBetheEnergy_UandQ(U+h, Q)
      Fm,Mm = CalcBetheEnergy_UandQ(U-h, Q)
      F0,M0 = CalcBetheEnergy_UandQ(U,Q)
      D = (Fp-Fm)/(2*h)

      print("{:8.3f}  {:12.6f}                              {:12.6f}".format(U, Nn/2., MuP))
      print("{:8.3f}  {:12.6f}  {:12.6f}  {:12.6f}  {:12.6f}  {:12.6f}".format(U, Nc/2., Ec, Gap, Mu, D))


# choose the computation to run.
Fn = [FindGaps, MakeSpinCorrelationCurves, MakeEnergyCurvesAtHalfFilling, MakeCurvesForQ, MakeRefDataAtHalfFilling][-1]
Fn()
