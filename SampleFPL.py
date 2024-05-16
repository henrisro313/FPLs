"""
Python3 code to generate FPL (Monte Carlo) and count its number of loops using NetworkX
author: Henrik Schou Roeising, 16/05/2024
"""

# Imports and packages:
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import sys, string, os, subprocess

plt.rcParams.update({'font.size':20})
plt.rcParams['mathtext.fontset']='stix'
plt.rcParams['font.family']='STIXGeneral'
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.preamble'] = r"\usepackage{bm} \usepackage{amsmath} \usepackage{amssymb}"

# System setup and NetworkX object:
mm, nn = 8,8
Nx, Ny = 2*mm, 2*nn
G = nx.grid_2d_graph(2*mm,2*nn) # square grid of 2m x 2n
edgeslist = list(G.edges()) # list with edges
faceliste = nx.minimum_cycle_basis(G) # list with faces

def Readfile_hor(filename, m, n):
    """
    Read horizontal bond variables from text file
    """
    infile = open(filename, 'r')
    lines = infile.readlines()
    hor = np.zeros(shape=(2*m,2*n-1),dtype=np.int32)
    i = 0
    for line in lines:
        vals = (line.strip()).split()
        for j in range(len(vals)):
            hor[j,i] = vals[j]
        i+= 1
    return hor

def Readfile_ver(filename, m, n):
    """
    Read vertical bond variables from text file
    """
    infile = open(filename, 'r')
    lines = infile.readlines()
    ver = np.zeros(shape=(2*m-1,2*n),dtype=np.int32)
    i = 0
    for line in lines:
        vals = (line.strip()).split()
        for j in range(len(vals)):
            ver[j,i] = vals[j]
        i+=1
    return ver

def count_loops(G):
     """
     Count the number of loops defined by the variables = +1 in graph G
     """
     eligible_edges = [(from_node,to_node,edge_attributes) for from_node,to_node,edge_attributes in G.edges(data=True) if edge_attributes['weight'] == 1]
     G_ups = nx.DiGraph()
     G_ups.add_edges_from(eligible_edges)
     no_loops = len(nx.cycle_basis(G_ups.to_undirected()))
     return no_loops

def count_length_loops(G):
     """
     Counting the average length of all loops (defined by the variables = +1) in G
     """
     L = 0.0
     eligible_edges = [(from_node,to_node,edge_attributes) for from_node,to_node,edge_attributes in G.edges(data=True) if edge_attributes['weight'] == 1]
     G_ups = nx.DiGraph()
     G_ups.add_edges_from(eligible_edges)
     LoopList = nx.cycle_basis(G_ups.to_undirected())
     NO = len(LoopList)
     for l in range(NO):
        L += len(LoopList[l])
     return L/NO

def plot_read(mm, nn, hor, ver):
    """
    Plot graph with horizontal and vertical bond variables hor and ver
    """
    hor = Readfile_hor('FinalState_hort.txt',m=mm,n=nn)
    ver = Readfile_ver('FinalState_vert.txt',m=mm,n=nn)
    s = 8
    cup = 'cornflowerblue'
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    for i in range(2*mm):
        for j in range(2*nn-1):
            if hor[i,j] == +1:
                plt.plot([i,i],[j,j+1],'-',color=cup,linewidth=4.0, markersize=s)
    for i in range(2*mm-1):
        for j in range(2*nn):
            if ver[i,j] == +1:
                plt.plot([i,i+1],[j,j],'-',color=cup,linewidth=4.0, markersize=s)
    for i in range(2*mm):
        for j in range(2*nn):
            plt.plot([i],[j],'.',color='mediumblue',linewidth=3.0, markersize=1.5*s)
    plt.axis('off')
    plt.tight_layout()
    ax.set_aspect('equal')
    #plt.savefig("LoopSample.pdf")
    plt.show()
    return None

up = {"weight":+1} # graph attributes associated with blue bonds
down = {"weight":-1} # graph attributes associated with white bonds

def ProduceGraph_countloops(mm, nn, hor, ver):
    """
    Produce a NetworkX graph of the variables contained in hor and ver, 
    and count the number of loops and average loop length
    """
    edattr_boxes = {}
    for ed in edgeslist:
        edattr_boxes[ed] = down # initialize all edges as white

    for i in range(2*mm):
        for j in range(2*nn-1):
            if hor[i,j] == +1:
                edattr_boxes[((i,j),(i,j+1))] = up
    for i in range(2*mm-1):
        for j in range(2*nn):
            if ver[i,j] == +1:
                edattr_boxes[((i,j),(i+1,j))] = up
    nx.set_edge_attributes(G, edattr_boxes) # associating the edge attributes in edattr_boxes to the graph G
    return G, count_loops(G), count_length_loops(G)

# Compile fortran90 code using subprocess:
subprocess.call(["gfortran","-c","MyDefs.f90"])#create
subprocess.call(["gfortran","FPL_MC_generate.f90","MyDefs.f90", "-o", "FPL_MC_generate.exe", "-ffree-line-length-none", "-llapack", "-O3"])

# Run fortran90 code to generate FPL output file:
os.system(r"./gfortran -c MyDefs.f90")
os.system(r"/Users/henrik/Documents/FPL_lengthdist/ForGitHub/FPL_MC_generate.exe "+str(mm)+" "+str(nn)+" 0.01 01")

# Read output files:
hor = Readfile_hor('FinalState_hort.txt',m=mm,n=nn)
ver = Readfile_ver('FinalState_vert.txt',m=mm,n=nn)

# Analyze result and create NetworkX object:
G, no, avlen = ProduceGraph_countloops(mm, nn, hor, ver)
print("Number of loops: ", no, " average length of loops: ", avlen)
plot_read(mm, nn, hor, ver)
