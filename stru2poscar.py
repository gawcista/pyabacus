from argparse import ArgumentParser
from operator import mod
from queue import Empty
from turtle import pos
from ase.units import Bohr
from ase import Atoms
from subprocess import getoutput
from numpy import array

parse = ArgumentParser()
parse.add_argument('-c', '--check',default="STRU",type=str)
parse.add_argument('-o', '--output',default="",type=str)
inputfile = parse.parse_args().check
outputfile = parse.parse_args().output

def readline(stringlist):
    flag = True
    line = ""
    while len(stringlist) > 0 and flag:
        stringline = stringlist.pop(0).split("#")[0]
        if len(stringline.split())>0:
            line=stringline
            flag = False
    return line,stringlist



def read_STRU_grep(filename="STRU"):
    POSCARlines = ["Transformed from %s"%filename]
    POSCARlines.append(float(getoutput('grep "LATTICE_CONSTANT" %s -A 1|tail -1'%(filename))))
    for i in range(3):
        POSCARlines.append('  '.join((array(getoutput('grep "LATTICE_VECTORS" %s -A %d |tail -1'%(filename,i+1)).split()[0:3]).astype(float)*Bohr).astype(str)))
    elements=[]
    numbers=[]
    pos = []
    positionline_number = int(getoutput('grep -n "ATOMIC_POSITIONS" %s'%(filename)).split(":")[0])
    positionline = getoutput('tail %s -n +%d'%(filename,positionline_number+1)).split("\n")
    mode,positionline=readline(positionline)
    while len(positionline)>0:
        element,positionline=readline(positionline)
        if element != "":
            elements.append(element)
            readline(positionline)#read magnetic moment
            natom,positionline=readline(positionline)
            numbers.append(natom)
            for i in range(int(natom)):
                position,positionline=readline(positionline)
                pos.append(position.replace("m",""))
    POSCARlines.append("  ".join(elements))
    POSCARlines.append("  ".join(numbers))
    POSCARlines.append("Selective Dynamics")
    POSCARlines.append(mode)
    for line in pos:
        POSCARlines.append(line)
    return POSCARlines
            
def printlines_list(lines,file=""):
    if file == "":
        for line in lines:
            print(line)
    else:
        with open(file,'w') as fout:
            for line in lines:
                print(line,file=fout)

printlines_list(read_STRU_grep(inputfile),outputfile)
    
