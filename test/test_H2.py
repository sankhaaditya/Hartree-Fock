import numpy as np
import sys
sys.path.insert(0, '../core')
import CGTOs
import basis
import rhf
import pp

H1s_exp = 1.24
H1s = CGTOs.getH1s()

title = 'H2, spacing: 0.1 a.u.'
label = 'STO-3G\nexp(H1s): '+str(H1s_exp)

basis_set = []
basis_set.append(basis.basis(H1s))
basis_set.append(basis.basis(H1s))

basis_per_nuclei = [1, 1]

N_charge = [1, 1]

e_count = 2

x = np.arange(0.1, 6.0, 0.1)

rhf = rhf.rhf(basis_set, basis_per_nuclei, N_charge, e_count, x)

E = rhf.solve()

pp = pp.pp(title, label)

print('\nBond Length: '+str(pp.min_energy(x, E)[0]))
print('\nTotal Energy: '+str(pp.min_energy(x, E)[1]))

pp.plot(x, E)
