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

# initial parameters
kB = C.physical_constants["Boltzmann constant in eV/K"][0]
h = C.physical_constants["Planck constant in eV s"][0]
mAtom = C.physical_constants["atomic mass constant"][0]
Ptot = 1


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
    Inp = open("reactions1", "r")
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
    SpeciesInp = open("species1", "r")
    species = {}
    concentration = {}
    for line in SpeciesInp.readlines():
        line = line.strip().split(" ")
        while '' in line:
            line.remove('')
        species[line[0]] = float(line[1])
        concentration[line[0]] = float(line[2])
    return species, concentration


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
        if len(reaction) == 2 and '*' in reaction[0]:
            g1, g2 = 0, 0
            for i in reaction[0]:
                g1 += species[i]
            for j in reaction[1]:
                g2 += species[j]
            DeltaG1 = g2 - g1
            DeltaG2 = g1 - g2
            rate1 = get_ads_k(DeltaG1, T, 28.)
            rate2 = get_des_k(DeltaG2, T, 28.)
        elif len(reaction) == 2 and '*' in reaction[1]:
            g1, g2 = 0, 0
            for i in reaction[0]:
                g1 += species[i]
            for j in reaction[1]:
                g2 += species[j]
            DeltaG1 = g2 - g1
            DeltaG2 = g1 - g2
            rate1 = get_des_k(DeltaG1, T, 28.)
            rate2 = get_ads_k(DeltaG2, T, 28.)
        elif len(reaction) == 3:
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


def output_k_g(reactions,k_G_list):
    for i in range(len(k_G_list)):
        print(reactions[i])
        print("G+:{:.5f}\tG-:{:.5f}\nk+:{:.3e}\tk-:{:.3e}".format(
            k_G_list[i][0], k_G_list[i][1], k_G_list[i][2], k_G_list[i][3]))

def output_equations(equations):
    for i in equations:
        print(i)
        for j in equations[i]:
            print(j)


T = 300
reactions = get_reactions()
species, concentration = get_species()
k_G_list = calc_k(reactions,species,T)
kBT = kB * T


# let the reactions list to 3*n matrix
for i in range(len(reactions)):
    if len(reactions[i]) == 2:
        reactions[i].insert(1,["->"])
output_k_g(reactions, k_G_list)

# construct concentration change
equations = {}
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
output_equations(equations)

adsorbates = []
adsorbatesq = []
for s in species:
    if s[-1] == "*" and equations[s] != []:
        adsorbates.append(s)
        adsorbatesq.append(adsorbates[-1].replace("*", "q"))
print(adsorbatesq)



# generate test.py for ODE solver
open('test.py', 'w').write("# -*- coding:utf-8 -*-\n\nfrom scipy.integrate import odeint\
\nimport numpy as np\nimport matplotlib.pyplot as plt\nfrom scipy import constants as C\nimport math\n\n\
def solver(w, t):\n\t")
wread = ", ".join(adsorbatesq)
open('test.py', 'a').write(wread)
open('test.py', 'a').write("= w\n\treturn [")

dsdt = []
for s in adsorbates:
    if equations[s] == []:  # ts
        continue
    if s == adsorbates[0]:
        continue
    elif s[-2:] == "_g":
        continue
    else:
        #print(s)
        dsdt.append([])
        #open('test.py', 'a').write("\n\t\t")
        for i in equations[s]:
            #print(i[0], i[2], end=" ")
            #open('test.py', 'a').write(i[0])
            dsdt[-1].append(i[0])
            #open('test.py', 'a').write(str(i[2]))
            dsdt[-1].append(str(i[2]))
            for j in i[1]:
                if j[-1] == "*":
                    jnew = j.replace("*", "q")
                    #if jnew == adsorbatesq[0]:
                    #    jnew = "(1 - " + adsorbatesq[1] + " - " + adsorbatesq[2] + ")"
                    #print("*", jnew, end=" ")
                    #open('test.py', 'a').write("*")
                    dsdt[-1].append("*")
                    #open('test.py', 'a').write(jnew)
                    dsdt[-1].append(jnew)
                elif j[-2:] == "_g":
                    #print("*", concentration[j], end=" ")
                    #open('test.py', 'a').write("*")
                    dsdt[-1].append("*")
                    #open('test.py', 'a').write(str(concentration[j]))
                    dsdt[-1].append(str(concentration[j]))
                else:
                    print("Error: species input")
        dsdt[-1] = "".join(dsdt[-1])
    #print("")

open('test.py', 'a').write("\n\t\t\n")
for i in dsdt:
    open('test.py', 'a').write("-(" + i + ")\n")
open('test.py', 'a').write(",")
for i in dsdt:
    open('test.py', 'a').write("\t\t\n")
    open('test.py', 'a').write(i+"\n,")

open('test.py', 'a').write("\t\t]\n")
open('test.py', 'a').write("t = np.logspace(-20., 15., 1000001)\n\
track = odeint(solver, [")
concentrationList = []
for s in adsorbates:
    concentrationList.append(str(concentration[s]))
open('test.py', 'a').write(",".join(concentrationList))
open('test.py', 'a').write("], t)\n\
print(track)\n")

open('test.py', 'a').write("conlist = []\n\
for i in track[-1]:\n\
\tconlist.append(i)\n")
open('test.py', 'a').write(wread)
open('test.py', 'a').write("= conlist\n\n")

for i in range(len(k_G_list)):
    open('test.py', 'a').write("\n"+"print(")
    open('test.py', 'a').write("+")
    open('test.py', 'a').write(str(k_G_list[i][2]))
    for j in reactions[i][0]:
        if j[-1] == "*":
            j = j.replace("*","q")
            open('test.py', 'a').write("*" + j)
        elif j[-2:] == "_g":
            open('test.py', 'a').write("*"+str(concentration[j]))

    open('test.py', 'a').write("-")
    open('test.py', 'a').write(str(k_G_list[i][3]))
    for j in reactions[i][2]:
        if j[-1] == "*":
            j = j.replace("*", "q")
            open('test.py', 'a').write("*" + j)
        elif j[-2:] == "_g":
            open('test.py', 'a').write("*"+str(concentration[j]))
    open('test.py', 'a').write(")"+"\n")
open('test.py', 'a').write("print(conlist)")
# plt draw picture
for i in range(len(adsorbates)):
    colorlib = ['b','g','r','c','m','y','k','w']
    open('test.py', 'a').write(
        "\nplt.plot(t, track[:, "+str(i)+"], '"+colorlib[i]+"', label=\""+ adsorbates[i] +"\")")

open('test.py', 'a').write(
    "\nplt.xscale('symlog')\
    \nplt.legend(loc='best')\
    \nplt.xlabel('t')\
    \nplt.show()")

os.system("python3 ./test.py")
