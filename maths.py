#!/usr/bin/env python
# encoding: utf-8
"""
maths.py

Created by Alfredo Di Napoli on 2011-03-29.
"""

from safe_io import *
from itertools import *


def get_connection_arcs(filename):
    """Returns the arcs connecting the atoms of the given compounds.
    Example: [[1,2], [1,7], [2,8]] where 1,2,7,8 are atoms id."""
    
    arcs = []
    
    pdb_file = fetch_and_open(filename, "r")
           
    
    extend = arcs.extend
    for line in pdb_file:
        
        splitted_tuple = line.split()
        if(splitted_tuple[0] == "CONECT"):
            extend([[int(splitted_tuple[1]), int(el)] for el in splitted_tuple[2:]])
    
    
    #Filtering the arcs list if (x1,x2) : x1 > x2
    arcs = list(ifilter(lambda x: x[0] < x[1], arcs))
    pdb_file.close()
    
    return arcs
    
    
def get_graph(filename):
    """Returns the tuple(node,arcs) representing the compound as a graph."""

    parser = SafePDBParser()
    structure = parser.get_structure('molecule', filename)
    model = structure[0]
    chain = model[' ']
    residue = chain[0]

    nodes = [atom.get_coord().tolist() for atom in residue]

    return nodes, get_connection_arcs(filename)


def get_trigraph(filename):
    """Same ad get_graph but contains even the atom code."""

    parser = SafePDBParser()
    structure = parser.get_structure('molecule', filename)
    model = structure[0]
    chain = model[' ']
    residue = chain[0]

    nodes = [atom.get_coord().tolist() for atom in residue]
    labels = [atom.get_id() for atom in residue]

    return nodes, get_connection_arcs(filename), labels


def get_complex(filename):
    """Questa funzione dovrebbe essere talmente generica da
    visualizzare qualsiasi file .PDB"""

    parser = SafePDBParser()
    structure = parser.get_structure('molecule', filename)

    nodes = []
    labels = []
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                nodes.extend([atom.get_coord().tolist() for atom in residue])
                labels.extend([atom.get_id() for atom in residue])

    return nodes, labels



def get_atoms(filename):

    parser = SafePDBParser()
    structure = parser.get_structure('molecule', filename)
    model = structure[0]
    chain = model[' ']
    residue = chain[0]

    atoms = [atom for atom in residue]

    return atoms
    
    
def link_graphs(g1, g2):
    
    v1, e1, l1 = g1
    v2, e2, l2 = g2
    
    #Attacco v1 a v2
    v = v1 + v2
    
    #Bisogna rimappare tutti gli archi di v2 rispetto alle nuove posizioni
    offset = len(v1)
    e = e1 + list(imap(lambda x: [x[0] + offset, x[1] + offset], e2))
    #e = e1 + [[el[0] + offset, el[1] + offset] for el in e2]
    
    l = l1 + l2
    return [v,e,l]
    
    
if __name__ == '__main__':
    g1 = get_trigraph("ALA.pdb")
    g2 = get_trigraph("ARG.pdb")
    
    print g1
    print g2
    print link_graphs(g1,g2) 