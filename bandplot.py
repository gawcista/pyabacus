from argparse import ArgumentParser
from turtle import color
from numpy import array
import json

parse = ArgumentParser()
parse.add_argument('-f', '--fermi',default=0.0,type=float)
parse.add_argument('-b', '--band',default="BANDS_1.dat",type=str)
parse.add_argument('-k', '--kpt',default="KPT",type=str)
parse.add_argument('-c', '--conf',default="conf.json",type=str)
parse.add_argument('-e','--erange',default=[-1.0,1.0], type=float,nargs='+')
parse.add_argument('-o', '--output',default="band.png",type=str)

E_fermi = parse.parse_args().fermi
E_range = parse.parse_args().erange
BANDS_file = parse.parse_args().band
KPT_file = parse.parse_args().kpt
CONF_file = parse.parse_args().conf
OUTPUT_file = parse.parse_args().output

def readline(stringlist,return_annotation=False):
    flag = True
    line = ""
    while len(stringlist) > 0 and flag:
        strtmp = stringlist.pop(0).split("#")
        stringline = strtmp[0]
        annotation = strtmp[-1]
        if len(stringline.split())>0:
            line=stringline
            flag = False
    if return_annotation:
        return line,annotation,stringlist
    else:
        return line,stringlist

def read_CONF(file):
    with open(file,'r') as fin:
        data = json.load(fin)
    Efermi = data["efermi"]
    energy_range = data["energy_range"]
    kptfile = data["kptfile"]
    bandfile = data["bandfile"]
    return bandfile,kptfile,Efermi,energy_range
    

def read_KPT(file):
    with open(file,'r') as fin:
        lines = fin.readlines()
    # Read the keyword tag
    keyword,lines = readline(lines)
    nk,lines =  readline(lines)
    mode,lines =  readline(lines)
    special_k_labels = []
    special_k_index = []
    n_kpt = 0
    if len(mode)>0 and mode.split()[0] == "Line":
        for i in range(int(nk)):
            line,label,lines = readline(lines,return_annotation=True)
            # Add labels of special k points
            if len(label.split())>0:
                special_k_labels.append(label.split()[0])
            else:
                special_k_labels.append("")
            # 
            special_k_index.append(n_kpt)
            n_kpt += int(line.split()[-1])
    return special_k_index,special_k_labels

def read_BANDS(file,E_fermi=0.0):
    kpath = []
    with open(file,'r') as fin:
        lines = fin.readlines()
    band = []
    for line in lines:
        kpath.append(float(line.split()[1]))
        band.append(array(line.split()[2:]).astype(float))
    band=array(band)
    return kpath,band.T-E_fermi

def plot_matplotlib(kpath,bands,special_k_index=[],special_k_labels=[],erange=[-1.0,1.0],outputfile="band.png"):
    import matplotlib.pyplot as plt
    border_width = 3
    plt.rcParams['font.family']=['Times New Roman']
    plt.rcParams['font.size'] = 18
    plt.rcParams['xtick.major.pad'] = 8
    plt.rcParams['ytick.major.pad'] = 8
    plt.rcParams['figure.autolayout'] = True
    fig, ax = plt.subplots(1, 1)
    fig.set_dpi(600)
    ax.spines['bottom'].set_linewidth(border_width)
    ax.spines['left'].set_linewidth(border_width)
    ax.spines['top'].set_linewidth(border_width)
    ax.spines['right'].set_linewidth(border_width)
    ax.tick_params(direction='in', width=border_width)
    ax.set_xlim(kpath[0],kpath[-1])
    ax.set_ylim(erange[0],erange[1])
    ax.set_ylabel('Energy (eV)',fontsize=24)
    xticks_position=[]
    xticks_label=[]
    for i in range(len(special_k_index)):
        xticks_position.append(kpath[special_k_index[i]])
        ax.vlines(kpath[special_k_index[i]],erange[0],erange[1],linewidth=1.0,edgecolor="gray",linestyles="--")
        if special_k_labels[i]=="G":
            xticks_label.append("Γ")
        else:
            xticks_label.append(special_k_labels[i])
    ax.set_xticks(xticks_position)
    ax.set_xticklabels(xticks_label)
    ax.hlines(0.,kpath[0],kpath[-1],linewidth=1.0,edgecolor="gray",linestyles="--")
    for band in bands:
        ax.plot(kpath,band,color="blue")
    #plt.show()
    plt.savefig(outputfile,transparent=True)

    
        
def plot_gnuplot(KPT_file,BANDS_file,efermi=0.0,erange=[-1.0,1.0],outputfile=""):
    special_k_index,special_k_labels = read_KPT(KPT_file) 


if CONF_file != "":
    BANDS_file,KPT_file,E_fermi,E_range=read_CONF(CONF_file)
special_k_index,special_k_labels = read_KPT(KPT_file)
kpath,band = read_BANDS(BANDS_file,E_fermi)
plot_matplotlib(kpath,band,special_k_index,special_k_labels,E_range,OUTPUT_file)
