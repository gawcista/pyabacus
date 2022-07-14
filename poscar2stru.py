import os
from numpy import ones,zeros,array
from ase import io
from ase.units import Bohr,Hartree
from ase.formula import Formula
from ase.data import atomic_numbers, atomic_masses
from ase.constraints import FixScaled,FixAtoms
import argparse

flag_find_upf = False
flag_find_orb = False
parse = argparse.ArgumentParser()
parse.add_argument('-l', '--lcao',default=1,type=int)
parse.add_argument('-f', '--find',default=1,type=int)
parse.add_argument('-du', '--dir_upf',default="./",type=str)
parse.add_argument('-do', '--dir_orb',default="./",type=str)
parse.add_argument('-o', '--output',default=1,type=int)
parse.add_argument('-of', '--outputfile',default="STRU",type=str)
parse.add_argument('-c', '--check',default="POSCAR",type=str)
file_in = parse.parse_args().check
flag_lcao = bool(parse.parse_args().lcao)
flag_find = bool(parse.parse_args().find)
if flag_find:
    flag_find_upf = True
    flag_find_orb = True
UPF_DIR = parse.parse_args().dir_upf
ORB_DIR = parse.parse_args().dir_orb
flag_output = bool(parse.parse_args().output)
file_output = parse.parse_args().outputfile
def read_POSCAR(file='POSCAR'):
    return io.read(file)

def find_file(filename,dir="."):
    result = []
    files = os.listdir(dir)
    for f in files:
        if filename in f:
            result.append(f)
    return result

def find_upf(element,dir="."):
    keytitle = element+'_'
    result = find_file(keytitle,dir)
    return result[0]

def find_orb(element,dir="."):
    keytitle = element+'_'
    result = find_file(keytitle,dir)
    return result[0]

def get_formula(symbol_list):
    symbols,nions = [],[]
    count = 1
    symbols.append(symbol_list[0])
    for i in range(1,len(symbol_list)):
        if symbol_list[i] == symbol_list[i-1]:
            count +=1
        else:
            nions.append(count)
            symbols.append(symbol_list[i])
            count = 1
    nions.append(count)
    return symbols,nions

def print_STRU(structure,filename='STRU'):
    symbols,nions = get_formula(structure.get_chemical_symbols())
    lattice_vec = structure.get_cell()[:]/Bohr
    nions_tot = len(structure)
    constraints = ones([nions_tot,3]).astype(int)
    for constraint in structure.constraints:
        if type(constraint) == FixScaled:
            constraints[constraint.todict()['kwargs']['a']] = array(constraint.todict()['kwargs']['mask']).astype(int)
        elif type(constraint) == FixAtoms:
            constraints[constraint.index[0]] = zeros(3).astype(int)
    with open(filename,'w') as stru_file:
        if flag_output:
            fout = stru_file
        else:
            fout = None
        print("ATOMIC_SPECIES",file=fout)
        for element in symbols:
            if flag_find_upf:
                print("%s %.3f %s"%(element,atomic_masses[atomic_numbers[element]],find_upf(element,UPF_DIR)),file=fout)
            else:
                print("%s %.3f %s_ONCV_PBE-1.0.upf"%(element,atomic_masses[atomic_numbers[element]],element),file=fout)
        print(file=fout)

        if flag_lcao:
            print("NUMERICAL_ORBITAL",file=fout)
            for element in symbols:
                if flag_find_orb:
                    print(find_orb(element,ORB_DIR),file=fout)
                else:
                    print("%s %.3f %s_ONCV_PBE-1.0.upf"%(element,atomic_masses[atomic_numbers[element]],element),file=fout)
        print(file=fout)
        
        print("LATTICE_CONSTANT",file=fout)
        print("1.0",file=fout)
        print(file=fout)
        print("LATTICE_VECTORS",file=fout)
        for vec in lattice_vec:
            print("%f\t%f\t%f"%(vec[0],vec[1],vec[2]),file=fout)
        print(file=fout)
        print("ATOMIC_POSITIONS",file=fout)
        print("Direct",file=fout)
        print(file=fout)
        i = 0
        for k in range(len(symbols)):
            print(symbols[k],file=fout)
            print('0.0',file=fout)
            print(nions[k],file=fout)
            for j in range(nions[k]):
                pos = structure[i].scaled_position
                print("%.10f\t%.10f\t%.10f\t %d %d %d"%(pos[0],pos[1],pos[2],constraints[i,0],constraints[i,1],constraints[i,2]),file=fout)
                i += 1
            print(file=fout)

atom=read_POSCAR(file_in)
print_STRU(atom,file_output) 
