import bethe_ansatz as ba
import numpy as np

E = []
U = 8
for N in np.linspace(0.8, 1.2, 11):
    e, n = ba.CalcBetheEnergy_UandN(U, N)
    E.append('%1.2f         %2.10f'%(N, e))

for i in range(len(E)):
    print E[i]
