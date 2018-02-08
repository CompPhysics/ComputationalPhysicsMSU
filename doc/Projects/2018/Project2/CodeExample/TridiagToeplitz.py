# Diagonalizing a Tridiagonal Toeplitz matrix
from math import cos, pi, log10
import numpy as np
#Function for initialization of parameters
def initialize():
    RMin = 0.0
    RMax = 1.0
    Dim = 20
    return RMin, RMax, Dim

#Get the boundary, orbital momentum and number of integration points
RMin, RMax, Dim = initialize()

#Initialize constants
Step    = RMax/(Dim+1)
DiagConst = 2.0 / (Step*Step)
NondiagConst =  -1.0 / (Step*Step)

#Calculate array of potential values
v = np.zeros(Dim)
r = np.linspace(RMin,RMax,Dim)

#Setting up a tridiagonal matrix and finding eigenvectors and eigenvalues
Hamiltonian = np.zeros((Dim,Dim))
Hamiltonian[0,0] = DiagConst
Hamiltonian[0,1] = NondiagConst
for i in range(1,Dim-1):
    Hamiltonian[i,i-1]  = NondiagConst
    Hamiltonian[i,i]    = DiagConst
    Hamiltonian[i,i+1]  = NondiagConst
Hamiltonian[Dim-1,Dim-2] = NondiagConst
Hamiltonian[Dim-1,Dim-1] = DiagConst
# diagonalize and obtain eigenvalues, not necessarily sorted
EigValues, EigVectors = np.linalg.eig(Hamiltonian)
# sort eigenvectors and eigenvalues
permute = EigValues.argsort()
EigValues = EigValues[permute]
EigVectors = EigVectors[:,permute]
# now plot the results for the three lowest lying eigenstates
for i in range(Dim):
    lambda_i = DiagConst+2*NondiagConst*cos((i+1)*pi*Step)
    print(EigValues[i], lambda_i)
