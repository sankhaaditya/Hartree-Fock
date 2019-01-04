import numpy as np
import integration as integ
import constants as ct

class rhf():

    def __init__(self, basis_set, basis_per_nuclei, N_charge, e_count, x):

        self.basis_set = basis_set
        self.basis_per_nuclei = basis_per_nuclei
        self.N_charge = N_charge
        self.e_count = e_count
        self.x = x

    def solve(self):
        
        basis_per_nuclei = self.basis_per_nuclei
        N_charge = self.N_charge
        e_count = self.e_count
        x = self.x

        E = []

        for rad in x:
            
            N_coord = [[0.0, 0.0, 0.0], [rad, 0.0, 0.0]]    # ONLY DIATOMIC!
            self.N_coord = N_coord

            if len(N_coord) != len(basis_per_nuclei):
                print('Length mismatch: N_coord, basis_per_nuclei')
                break

            count = 0
            for i in range(0, len(basis_per_nuclei)):
                for j in range(0, basis_per_nuclei[i]):
                    self.basis_set[count].setCoords(N_coord[i])
                    count += 1

            energy = self.calc_energy()

            nuclear = 1.0
            for i in range(0, len(N_charge)):
                nuclear *= N_charge[i]

            energy_total = energy + nuclear / rad

            E.append(energy_total)

        return E

    def calc_energy(self):

        e_count = self.e_count

        C, E, H = self.scf()

        energy = 0.0
        for level in range(0, int(e_count / 2)):
            for i in range(0, len(C)):
                for j in range(0, len(C)):
                    energy += C[i][level] * C[j][level] * H[i][j]
            energy += E[level]

        return energy

    def scf(self):

        basis_set = self.basis_set
        lbs = len(basis_set)
        e_count = self.e_count
        self.Val_ee = np.zeros([lbs, lbs, lbs, lbs])

        S, H = self.init_matrix()
        F = H
        
        ## Courtesy of Pu Du <https://github.com/ipudu/SCFpy>
        s, U = np.linalg.eig(S)
        s_mhalf = np.diag(s**(-0.5))
        S_mhalf = np.dot(U, np.dot(s_mhalf, U.T))
        X = S_mhalf
        ##

        count = 0
        error = 1.0
        while error > ct.tolerance:

            ## Courtesy of Pu Du <https://github.com/ipudu/SCFpy>
            F_prime = np.dot(X.T, np.dot(F, X))
            E, C_prime = np.linalg.eigh(F_prime)
            C = np.dot(X, C_prime)
            ##
            
            if count != 0:
                error = np.linalg.norm(E-E_prev) / np.linalg.norm(E_prev)

            E_prev = E
            count += 1

            if count > ct.max_iterations:
                break;

            V_ee = self.corr_matrix(C)
            F = H + V_ee

        return C, E, H

    def init_matrix(self):

        basis_set = self.basis_set
        lbs = len(basis_set)
        N_coord = self.N_coord
        N_charge = self.N_charge

        S = np.zeros([lbs, lbs])
        H = np.zeros([lbs, lbs])

        for i in range(0, lbs):
            for j in range(0, lbs):
                for k in range(0, ct.basis_size):
                    for l in range(0, ct.basis_size):
                        a = basis_set[i].CGTO.infoList[k][0]
                        b = basis_set[j].CGTO.infoList[l][0]
                        pa = basis_set[i].CGTO.powers
                        pb = basis_set[j].CGTO.powers
                        ca = basis_set[i].N_coords
                        cb = basis_set[j].N_coords

                        coeff = basis_set[i].CGTO.infoList[k][1] * \
                                  basis_set[j].CGTO.infoList[l][1] * \
                                  np.power(a * b, 0.75) * np.power(2 / np.pi, 1.5)
                        
                        overlap = integ.overlap(a, b, ca, cb, pa, pb)
                        S[i][j] += overlap * coeff

                        kinetic = integ.kinetic(a, b, ca, cb, pa, pb)
                        H[i][j] += kinetic * coeff

                        for m in range(0, len(N_coord)):
                            nuclear = integ.nuclear(a, b, ca, cb, pa, pb, N_coord[m])
                            H[i][j] -= nuclear * coeff * N_charge[m]

        return S, H

    def corr_matrix(self, C):

        basis_set = self.basis_set
        lbs = len(basis_set)
        V_ee = np.zeros([lbs, lbs])
        Val_ee = self.Val_ee
        e_count = self.e_count
            
        for level in range(0, int(e_count / 2)):
            for p in range(0, lbs):
                for q in range(0, lbs):
                    for r in range(0, lbs):
                        for s in range(0, lbs):

                            # <pr|qs>
                            if Val_ee[p, q, r, s] == 0.0:
                                Val_ee[p, q, r, s] = self.twoe(p, q, r, s)
                                Val_ee[q, p, r, s] = Val_ee[p, q, r, s]
                                Val_ee[p, q, s, r] = Val_ee[p, q, r, s]
                                Val_ee[q, p, s, r] = Val_ee[p, q, r, s]
                                Val_ee[r, s, p, q] = Val_ee[p, q, r, s]
                                Val_ee[s, r, p, q] = Val_ee[p, q, r, s]
                                Val_ee[r, s, q, p] = Val_ee[p, q, r, s]
                                Val_ee[s, r, q, p] = Val_ee[p, q, r, s]

                            V_ee[p][q] += 2 * C[r][level] * C[s][level] * Val_ee[p, q, r, s]
                            V_ee[p][q] -= C[r][level] * C[s][level] * Val_ee[p, r, q, s]

        return V_ee

    def twoe(self, p, q, r, s):

        basis_set = self.basis_set
        val = 0.0

        for i in range(0, ct.basis_size):
            for j in range(0, ct.basis_size):
                for k in range(0, ct.basis_size):
                    for l in range(0, ct.basis_size):

                        a = basis_set[p].CGTO.infoList[i][0]
                        b = basis_set[q].CGTO.infoList[j][0]
                        c = basis_set[r].CGTO.infoList[k][0]
                        d = basis_set[s].CGTO.infoList[l][0]
                        pa = basis_set[p].CGTO.powers
                        pb = basis_set[q].CGTO.powers
                        pc = basis_set[r].CGTO.powers
                        pd = basis_set[s].CGTO.powers
                        ca = basis_set[p].N_coords
                        cb = basis_set[q].N_coords
                        cc = basis_set[r].N_coords
                        cd = basis_set[s].N_coords

                        coeff = basis_set[p].CGTO.infoList[i][1] * \
                            basis_set[q].CGTO.infoList[j][1] * \
                            basis_set[r].CGTO.infoList[k][1] * \
                            basis_set[s].CGTO.infoList[l][1] * \
                            np.power(a * b * c * d, 0.75) * \
                            np.power(2 / np.pi, 3)
                        
                        val += integ.electron(a, pa, ca, b, pb, cb, c, pc, cc, d, pd, cd) \
                               * coeff

        return val
