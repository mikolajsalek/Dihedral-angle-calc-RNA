import Bio
import numpy as np
from Bio.PDB import *
#import matplotlib.pyplot as plt
import sys

pdb_list = PDBList()
print(pdb_list)

arg1 = sys.argv[1] #nazwa pliku

#pobieranie pliku
pdb_id = arg1

pdb_filename = pdb_list.retrieve_pdb_file(pdb_id, pdir="data/PDB_files", file_format="pdb")
structure = Bio.PDB.PDBParser().get_structure(pdb_id, pdb_filename)

#zapis szukanych kątów
beta = ['P', "O5'", "C5'", "C4'"]
gamma = ["O5'", "C5'", "C4'", "C3'"]
delta = ["C5'", "C4'", "C3'", "O3'"]
chi_pur = ["O4'", "C1'", 'N9', 'C4']
chi_pir = ["O4'", "C1'", 'N1', 'C2']
alpha = ['P', "O5'", "C5'"]
alpha_prev = ["O3'"]
epsilon = ["C4'", "C3'", "O3'"]
epsilon_next = ['P']
zeta = ["C3'", "O3'"]
zeta_next = ["P", "O5'"]

#stworzenie listy reszt i atomow
res_list = Selection.unfold_entities(structure, "R")
atm_list = []
residue_list = []

#usuwanie hetatm
for r in res_list:
    tags = r.get_full_id()
    if tags[3][0] == " ":
        #print(r)
        residue_list.append(r)

#tworzenie list gdzie zostaną zapisane atomy do liczenia kątów torsyjnych
beta_selected = []
gamma_selected = []
delta_selected = []
chi_pur_selected = []
chi_pir_selected = []
alpha_selected = []
epsilon_selected = []
zeta_selected = []

matrix_row = [] #do tworzenia matrixa
matrix = []

#stworzenie list atomów dla danej reszty
for i in range(0, len(residue_list)):
    for atm in residue_list[i]:
       atm_list.append(atm)
    #tworzenie list z koordynatami atomow do danych kątów
    if i != 0:
        r = residue_list[i - 1]
        for x in r:
            if x.get_id() in alpha_prev:
                alpha_selected.append(x.get_vector())

    for atom in atm_list:
        if atom.get_id() in beta:
            beta_selected.append(atom.get_vector())
        if atom.get_id() in gamma:
            gamma_selected.append(atom.get_vector())
        if atom.get_id() in delta:
            delta_selected.append(atom.get_vector())
        if atom.get_id() in chi_pur:
            chi_pur_selected.append(atom.get_vector())
        if atom.get_id() in chi_pir:
            chi_pir_selected.append(atom.get_vector())
        if atom.get_id() in alpha:
            alpha_selected.append(atom.get_vector())
        if atom.get_id() in epsilon:
            epsilon_selected.append(atom.get_vector())
        if atom.get_id() in zeta:
            zeta_selected.append(atom.get_vector())

    if i != len(residue_list) - 1:
        r = residue_list[i + 1]
        for x in r:
            if x.get_id() in epsilon_next:
                epsilon_selected.append(x.get_vector())

    if i != len(residue_list) - 1:
        r = residue_list[i + 1]
        for x in r:
            if x.get_id() in zeta_next:
                zeta_selected.append(x.get_vector())

    #print(zeta_selected)
    #print("-------------")
    matrix_row.append(residue_list[i].get_resname())

    if len(alpha_selected) == 4:
        angle_alpha = calc_dihedral(alpha_selected[0], alpha_selected[1], alpha_selected[2], alpha_selected[3])
        angle_alpha_degree = angle_alpha * (180.0 / np.pi)
        matrix_row.append(angle_alpha_degree)
    else:
        matrix_row.append("brak")

    if len(beta_selected) == 4:
        angle_beta = calc_dihedral(beta_selected[0], beta_selected[1], beta_selected[2], beta_selected[3])
        angle_beta_degree = angle_beta * (180.0/np.pi)
        matrix_row.append(angle_beta_degree)
    else:
        matrix_row.append("brak")

    if len(gamma_selected) == 4:
        angle_gamma = calc_dihedral(gamma_selected[0], gamma_selected[1], gamma_selected[2], gamma_selected[3])
        angle_gamma_degree = angle_gamma * (180.0 / np.pi)
        matrix_row.append(angle_gamma_degree)
    else:
        matrix_row.append("brak")

    if len(delta_selected) == 4:
        angle_delta = calc_dihedral(delta_selected[0], delta_selected[1], delta_selected[2], delta_selected[3])
        angle_delta_degree = angle_delta * (180.0 / np.pi)
        matrix_row.append(angle_delta_degree)
    else:
        matrix_row.append("brak")

    if len(epsilon_selected) == 4:
        angle_epsilon = calc_dihedral(epsilon_selected[0], epsilon_selected[1], epsilon_selected[2], epsilon_selected[3])
        angle_epsilon_degree = angle_epsilon * (180.0 / np.pi)
        matrix_row.append(angle_epsilon_degree)
    else:
        matrix_row.append("brak")

    if len(zeta_selected) == 4:
        angle_zeta = calc_dihedral(zeta_selected[0], zeta_selected[1], zeta_selected[2], zeta_selected[3])
        angle_zeta_degree = angle_zeta * (180.0 / np.pi)
        matrix_row.append(angle_zeta_degree)
    else:
        matrix_row.append("brak")

    if residue_list[i].get_resname() == "G" or residue_list[i].get_resname() == "A":
        if len(chi_pur_selected) == 4:
            angle_chi_pur = calc_dihedral(chi_pur_selected[0], chi_pur_selected[1], chi_pur_selected[2],
                                          chi_pur_selected[3])
            angle_chi_pur_degree = angle_chi_pur * (180.0 / np.pi)
            matrix_row.append(angle_chi_pur_degree)
        else:
            matrix_row.append("brak")

    if residue_list[i].get_resname() == "C" or residue_list[i].get_resname() == "U":
        if len(chi_pir_selected) == 4:
            angle_chi_pir = calc_dihedral(chi_pir_selected[0], chi_pir_selected[1], chi_pir_selected[2],
                                          chi_pir_selected[3])
            angle_chi_pir_degree = angle_chi_pir * (180.0 / np.pi)
            matrix_row.append(angle_chi_pir_degree)
        else:
            matrix_row.append("brak")

    matrix.append(matrix_row.copy())
    #print(matrix_row)
    #czyszczenie tablic
    alpha_selected.clear()
    matrix_row.clear()
    beta_selected.clear()
    gamma_selected.clear()
    delta_selected.clear()
    chi_pur_selected.clear()
    chi_pir_selected.clear()
    atm_list.clear()
    zeta_selected.clear()
    epsilon_selected.clear()

file = open("matrix.txt", "w")
with open("matrix.txt", "w") as file:
    # Write header
    file.write("Residue\tAlpha\tBeta\tGamma\tDelta\tEpsilon\tZeta\tChi_Pur\tChi_Pir\n")

    for row in matrix:
        file.write("\t".join(map(str, row)) + "\n")


file.close()
