from bethe_ansatz import CalcBetheEnergy_UandQ
import numpy as np
import sys



def print_mat(mat):
    s = mat.shape
    for i in range(s[0]):
        for j in range(s[1]):
            sys.stdout.write('%2.6f     '%E[i,j])
        sys.stdout.write('\n')
    sys.stdout.flush()

E = np.zeros([13+1, 7+1])
U = np.linspace(0.0, 12.0, 13)
nf = np.linspace(0.4,1.6,7)
E[0, 1:] = nf
E[1:, 0] = U.T
for i in range(13):
    for j in range(7):
        e, n =  CalcBetheEnergy_UandQ(U[i], nf[j])
        E[i+1,j+1] = e

print_mat(E)
