# -*- coding:utf-8 -*-

"""
@author: Jin-Cheng Liu
@file: mk
@time: 2018/11/26 22:58
"""

from scipy.integrate import odeint
import numpy as np
from scipy import constants as C
import math,os
masslist= [1.00794,4.002602,6.941,9.012182,10.811,12.011,14.00674,\
     15.9994,18.9984032,20.1797,22.989768,24.3050,26.981539,28.0855,\
     30.97362,32.066,35.4527,39.948,39.0983,40.078,44.955910,47.88,\
      50.9415,51.9961,54.93085,55.847,58.93320,58.69,63.546,65.39,\
      69.723,72.61,74.92159,78.96,79.904,83.80,85.4678,87.62,88.90585,91.224,\
      92.90638,95.94,98,101.07,102.90550,106.42,107.8682,112.411,\
      114.82,118.710,121.75,127.60,126.90447,131.29,132.90543,137.327,138.9055,140.115,\
      140.90765,144.24,145,150.36,151.965,157.25,158.92534,162.50,\
      164.93032,167.26,168.93421,173.04,174.967,178.49,180.9479,\
      183.85,186.207,190.2,192.22,195.08,196.96654,200.59,204.3833,\
      207.2,208.98037,209,210.0,222,223,226.025,227.028,232.0381,\
      231.03588,238.0289,237.048,244.,243.,247.,247.,251.,252.,\
      257.,258.,259.,260.0]

elementlist = ['H','He','Li','Be','B','C','N','O','F','Ne',\
     'Na','Mg','Al','Si','P','S','Cl','Ar',\
     'K','Ca','Sc','Ti','V','Cr','Mn','Fe',\
     'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',\
     'Rr','Sr','Y','Zr','Nb','Mo','Tc','Ru',\
     'Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe',\
     'Cs','Ba','La','Ce','Pr','Nd','Pm',\
     'Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W',\
     'Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',\
     'FR','RA','AC','TH','PA','U','NP',\
     'PU','AM','CM','BK','CF','ES','FM','MD','NO','LR']

massDict = dict(zip(elementlist, masslist))

# Constants
kB = C.physical_constants["Boltzmann constant in eV/K"][0]
h = C.physical_constants["Planck constant in eV s"][0]
mAtom = C.physical_constants["atomic mass constant"][0]

# Initial parameters
Ptot = 1        # Pressure in atm
T = 1363         # Temperature in K
kBT = kB * T

def get_react_k(deltaG, T):
    '''
    calculate rate constant for surface catalytic reaction.
    :param deltaG: deltaG = G(transition state) - G(initial state) or
     deltaG = G(transition state) - G(final state)
    :param T: reaction temperature
    :return: reaction rate constant at T
    '''
    kBT = kB * T
    rate = (kBT / h) * (math.e**(-deltaG/kBT))
    return rate


def get_ads_k(deltaG, T, m, P=101325,  S=1., A=9E-20):
    '''
    calculate rate constant for surface catalytic reaction.
    :param deltaG: G(ads) - G(gas), unit: eV
    :param T: reaction temperature, unit: K
    :param P: reaction pressure, unit:Atm
    :param m: molecular mass, unit； ,
    :param S: sticking coefficient, here 1 as an example, need test
    :param A: active site area, unit:m**2, here 4 Angstrom**2 as an example, need test
    :return: adsorption rate constant at T
    '''
    P = P/Ptot
    m = m*mAtom
    kBT = kB * T
    kads = S*P*A/math.sqrt(2*math.pi*m*kBT)  # no P in k
    return kads


def get_des_k(deltaG, T, m, P=101325, S=1., A=9E-20):
    '''
    :param deltaG: G(gas) - G(ads), unit: eV
    :param T:
    :param P:
    :param m:
    :param S:
    :param A:
    :return: desorption rate constant at T
    '''
    kBT = kB * T
    kads = get_ads_k(-deltaG, T, m, P, S, A)
    kdes = (math.e**(-deltaG/kBT))*kads
    return kdes


def get_reactions():
    '''
    read reactions from file
    :return: reactions list
    '''
    Inp = open("REACTION", "r")
    reactions = []
    for line in Inp.readlines():
        line = line.strip().split("<->")
        while '' in line:
            line.remove('')
        line = list(map(lambda x: x.strip().split("+"), line))
        for i in range(len(line)):
            line[i] = list(map(lambda x: x.strip(), line[i]))
            for j in line[i]:
                if str.isdigit(j[0]):
                    for _ in range(int(j[0])):
                        line[i].append(j[1:])
                    line[i].remove(j)
        reactions.append(line)
    return reactions


def get_species():
    '''
    read species information file
    :return: species list
    '''
    SpeciesInp = open("ENERGY", "r")
    species = {}
    concentration = {}
    for line in SpeciesInp.readlines():
        line = line.strip().split(" ")
        while '' in line:
            line.remove('')
        species[line[0]] = float(line[1])
        concentration[line[0]] = float(line[2])
    return species, concentration


def get_mass(r):
    '''
    From reaction list get gas molecule mass. Note there must be only one gas molecule.
    :param r: reactant list ex: ['H2_g', '*', '*']
    :return: gas mass for rate calculation
    '''
    for i in r:  # calc mass
        if i != "*":
            mass = 0
            for j in range(len(i[:-2])):
                if i[:-2][j] in massDict:
                    mass += massDict[i[:-2][j]]
                else:
                    num = int(i[:-2][j])
                    mass += massDict[i[:-2][j - 1]] * (num - 1)
            print("mass", i, mass)
    return mass


def calc_k(reactions,species,T):
    '''
    calculate Delta G and reaction rate
    :param reactions:
    :param species:
    :param T:
    :return: a list with [DG, reversed DG, rate k, reversed rate k]
    '''
    k_G_list = []
    for reaction in reactions:
        if len(reaction) == 2 and '*' in reaction[0]:   # Adsorption
            g1, g2 = 0, 0
            for i in reaction[0]:
                g1 += species[i]
            for j in reaction[1]:
                g2 += species[j]
            DeltaG1 = g2 - g1
            DeltaG2 = g1 - g2
            mass = get_mass(reaction[0])
            rate1 = get_ads_k(DeltaG1, T, mass)
            rate2 = get_des_k(DeltaG2, T, mass)
        elif len(reaction) == 2 and '*' in reaction[1]: # Desorption
            g1, g2 = 0, 0
            for i in reaction[0]:
                g1 += species[i]
            for j in reaction[1]:
                g2 += species[j]
            DeltaG1 = g2 - g1
            DeltaG2 = g1 - g2
            mass = get_mass(reaction[1])
            rate1 = get_des_k(DeltaG1, T, mass)
            rate2 = get_ads_k(DeltaG2, T, mass)
        elif len(reaction) == 3:                        # Reaction
            g1, gts, g2 = 0, 0, 0
            for i in reaction[0]:
                g1 += species[i]
            for j in reaction[1]:
                gts += species[j]
            for k in reaction[2]:
                g2 += species[k]
            DeltaG1 = gts - g1
            DeltaG2 = gts - g2
            rate1 = get_react_k(DeltaG1, T)
            rate2 = get_react_k(DeltaG2, T)
        else:
            print("Error?????")
        k_G_list.append([DeltaG1,DeltaG2,rate1,rate2])
    return k_G_list


def convertReactions(reactions):
    for i in range(len(reactions)):
        if len(reactions[i]) == 2:
            reactions[i].insert(1, ["->"])
    return reactions


def output_k_g(reactions,k_G_list):
    print("\t###################################################################################################\n\
    Reaction rate and Delta G calculations for elementary reaction steps in microkinetic modeling\n\
    ###################################################################################################")
    for i in range(len(k_G_list)):
        print(reactions[i])
        print("G+:{:10.5f}\tG-:{:10.5f}\nk+:{:10.3e}\tk-:{:10.3e}".format(
            k_G_list[i][0], k_G_list[i][1], k_G_list[i][2], k_G_list[i][3]))


def construct_equations(equations={}):
    for s in species:
        equations[s] = []
        for i in range(len(reactions)):
            for j in reactions[i][0]:
                if s == j:
                    equations[s].append(["-", reactions[i][0], k_G_list[i][2]])
                    equations[s].append(["+", reactions[i][2], k_G_list[i][3]])
            for j in reactions[i][2]:
                if s == j:
                    equations[s].append(["+", reactions[i][0], k_G_list[i][2]])
                    equations[s].append(["-", reactions[i][2], k_G_list[i][3]])
    return equations


def output_equations(equations):
    print("\t###################################################################################################\n\
    Differential equations for elementary reaction steps in microkinetic modeling\n\
    ###################################################################################################")
    for i in equations:
        print("d", i, "/ dt =")
        for j in equations[i]:
            print(j[0], "*", j[2], "*", j[1])


def construct_AdsbatesqList(adsorbates = [], adsorbatesq = []):
    for s in species:
        if s[-1] == "*" and equations[s] != []:
            adsorbates.append(s)
            adsorbatesq.append(adsorbates[-1].replace("*", "q"))
    return adsorbates, adsorbatesq

# get the reactions dict and species dict with concentration from input files
reactions = get_reactions()
species, concentration = get_species()

# calculate k and G
k_G_list = calc_k(reactions,species,T)
# let the reactions list to 3*n matrix
convertReactions(reactions)  # Note: this code is ugly
# print k and G
output_k_g(reactions, k_G_list)

# construct differential equations
equations = construct_equations({})
output_equations(equations)

# convert * to q
adsorbates, adsorbatesq = construct_AdsbatesqList()


# generate test.py for ODE solver
odeSolver = []
wread = ", ".join(adsorbatesq)
dsdt = []

odeSolver.append("# -*- coding:utf-8 -*-\n\nfrom scipy.integrate import odeint\
\nimport numpy as np\nimport matplotlib.pyplot as plt\nfrom scipy import constants as C\nimport math\n\n\
def solver(w, t):\n\t")
odeSolver.append(wread)
odeSolver.append("= w\n\treturn [")

for s in adsorbates:
    if equations[s] == []:  # ts
        continue
    if s == adsorbates[0]:
        continue
    elif s[-2:] == "_g":
        continue
    else:
        dsdt.append([])
        for i in equations[s]:
            dsdt[-1].append(i[0])
            dsdt[-1].append(str(i[2]))
            for j in i[1]:
                if j[-1] == "*":
                    jnew = j.replace("*", "q")
                    dsdt[-1].append("*")
                    dsdt[-1].append(jnew)
                elif j[-2:] == "_g":
                    dsdt[-1].append("*")
                    dsdt[-1].append(str(concentration[j]))
                else:
                    print("Error: species input")
        dsdt[-1] = "".join(dsdt[-1])

odeSolver.append("\n\t\t\n")
for i in dsdt:
    odeSolver.append("-(" + i + ")\n")
odeSolver.append(",")
for i in dsdt:
    odeSolver.append("\t\t\n")
    odeSolver.append(i+"\n,")

odeSolver.append("\t\t]\n")
odeSolver.append("t = np.logspace(-20., 15., 1000001)\n\
track = odeint(solver, [")
concentrationList = []
for s in adsorbates:
    concentrationList.append(str(concentration[s]))
odeSolver.append(",".join(concentrationList))
odeSolver.append("], t)\n\
print(track)\n")

odeSolver.append("conlist = []\n\
for i in track[-1]:\n\
\tconlist.append(i)\n")
odeSolver.append(wread)
odeSolver.append("= conlist\n\n")

for i in range(len(k_G_list)):
    odeSolver.append("\n"+"print(")
    odeSolver.append("+")
    odeSolver.append(str(k_G_list[i][2]))
    for j in reactions[i][0]:
        if j[-1] == "*":
            j = j.replace("*","q")
            odeSolver.append("*" + j)
        elif j[-2:] == "_g":
            odeSolver.append("*"+str(concentration[j]))

    odeSolver.append("-")
    odeSolver.append(str(k_G_list[i][3]))
    for j in reactions[i][2]:
        if j[-1] == "*":
            j = j.replace("*", "q")
            odeSolver.append("*" + j)
        elif j[-2:] == "_g":
            odeSolver.append("*"+str(concentration[j]))
    odeSolver.append(")"+"\n")
odeSolver.append("print(conlist)")
# plt draw picture
for i in range(len(adsorbates)):
    colorlib = ['b','g','r','c','m','y','k','w', 'b', 'b', 'b', 'b', 'b']
    odeSolver.append(
        "\nplt.plot(t, track[:, "+str(i)+"], '"+colorlib[i]+"', label=\""+ adsorbates[i] +"\")")

odeSolver.append(
    "\nplt.xscale('symlog')\
    \nplt.legend(loc='best')\
    \nplt.xlabel('t')\
    \nplt.show()")
with open('test.py','w') as writer:
    writer.writelines(odeSolver)

os.system("python3 ./test.py")
