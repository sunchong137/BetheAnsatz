import bethe_ansatz as ba
import numpy as np

def test1():
    N = 1.0
    E = []
    for U in np.linspace(0.05, 1.00, 20):
        e, n = ba.CalcBetheEnergy_UandN(U, N)
        E.append('%1.2f         %2.10f'%(U, e))
    
    for i in range(len(E)):
        print E[i]
def test2():
    N = 1.0
    U = 4.0
    e, n = ba.CalcBetheEnergy_UandN(U, N)
    print e

test2()
