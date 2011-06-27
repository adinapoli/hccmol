#!/usr/bin/env python
# encoding: utf-8
"""
graphics.py

Created by Alfredo Di Napoli on 2011-03-29.
"""

from atomic_radius import *
from maths import *
from itertools import *
from hasselib.graphLib6 import *


atom_color = {
    'H': Color4f([0.8, 0.8, 0.8, 1.0]), # ligth gray
    'C': Color4f([0.3, 0.3, 0.3, 1.0]), # dark gray (quite black)
    'N': BLUE,
    'O': RED,
    'F': Color4f([0.0, 0.75, 1.0, 1.0]), # ligth blue
    'P': ORANGE,
    'S': YELLOW,
    'Cl': GREEN,
    'K': Color4f([200./255, 162./255, 200./255, 1.0]) # lilac
}



def create_primary_structure(graph):
    """Returns a Hpc representing the primary structure of a chemical compound."""
    
    v, e, l = graph 
    return MKPOL([v,e,None])
    

def mkatom(atom, radius_type = 3, scale = 100.0):
    """Returns a graphic representation of the given atom.
    
    The representation is a Plasm Batch object.
    """
    
    coords, name = atom
    
    batch = get_sphere().next()
    atom_code = name[0]
    sf = atomic_radius[atom_code][radius_type]/scale
    
    batch.matrix = Mat4f.translate(*coords) * Mat4f.scale(sf, sf, sf) * batch.matrix
    batch.diffuse = atom_color[atom_code]
    
    return batch
    
    
def mkligand(graph):
    """Costruct an OpenGL representation of the given ligand name."""
    
    primary = create_primary_structure(graph)
    atoms_list = [mkatom(i,scale = 1000.0) for i in izip(graph[0], graph[2])]
    primary = Plasm.getBatches(primary)
    
    return atoms_list + list(primary)
    
    
def link(g1, g2):

    v1,e1,l1 = g1
    v2,e2,l2 = g2
    
    #First step, make the last oxygen and the first nitrogen collide
    ligand1 = map(lambda x: mkatom(x, scale = 1000.0), izip(v1,l1))
    ligand2 = map(lambda x: mkatom(x, scale = 1000.0), izip(v2,l2))

    
    #IL SECONDO LIGANDO RUOTA, DAL PRIMO PRENDO C-OXT
    v1_reversed = v1[::-1]
    l1_reversed = l1[::-1]
    
    
    O_index = l1_reversed.index("OXT")
    O_atom =  v1_reversed[O_index]
    C_atom = v1[2]
    
    #Dal secondo prendo N-H
    N_atom = v2[l2.index("N")]
    H_atom = v2[l2.index("H")]
    
    
    diff1 = VECTDIFF([[0.,0.,0.], O_atom])
    diff2 = VECTDIFF([[0.,0.,0.], N_atom])
     
    
    #Porto il primo ligando all'origine traslando di -N_atom    
    for batch in ligand1:
        batch.matrix = Mat4f.translate(*diff1) * batch.matrix
    
    #Porto il secondo ligando all'origine traslando di -O_atom
    for batch in ligand2:
        batch.matrix = Mat4f.translate(*diff2) * batch.matrix
     
    
    #Rotazione del primo portando H sull'asse C-O.
    h_vect = VECTDIFF([H_atom, N_atom])
    c_vect = VECTDIFF([C_atom, O_atom])
    axis = UNITVECT(VECTPROD([h_vect, c_vect]))
    angle = ACOS(INNERPROD([UNITVECT(h_vect), UNITVECT(c_vect)]))
    
    
    #Porto tutto il polipeptide creato nel sistema di riferimento
    #del SECONDO ligando.
    new_v2 = []
    append = new_v2.append
    for batch in ligand2:
        batch.matrix = Mat4f.translate(*N_atom)*Mat4f.rotate(Vec3f(*axis), angle)*batch.matrix
        append([batch.matrix.a14(), batch.matrix.a24(), batch.matrix.a34()])
    
        
    new_v1 = []
    append = new_v1.append   
    for batch in ligand1:
        batch.matrix = Mat4f.translate(*N_atom) * batch.matrix
        append([batch.matrix.a14(), batch.matrix.a24(), batch.matrix.a34()])
     

    return link_graphs([new_v1, e1, l1], [new_v2, e2, l2])

  
def mkpolypeps(peptides_list):

    graphs = [get_trigraph(filename + ".pdb") for filename in peptides_list]    
    result = reduce(link, graphs)
    
    return mkligand(result)


def view_compound(batches):
    
    if not isinstance(batches,list):
        batches = [batches]
    
    octree = Octree(Batch.Optimize(batches))
    viewer = Viewer(octree)
    viewer.Run()


def smart_hccbatches(atom, base = None):

    g = atom.graph
    color = atom.color
    coords = atom.coords
    sf = atom.scale

    if not base:
        h = Hpc(g)
        batches = Plasm.getBatches(h)

        for b in batches:
            b.diffuse = color
            b.matrix = Mat4f.translate(*coords) * Mat4f.scale(sf, sf, sf) * b.matrix

        return batches

    else:
        b = Batch(base)
        b.diffuse = color
        b.matrix = Mat4f.translate(*coords) * Mat4f.scale(sf, sf, sf) * b.matrix
        return b


def MOLVIEW(batches):

    """Non faccio altro che visualizzare le batches. Ho
    separato la parte di visualizzazione da CHE COSA visualizzare."""

    if not isinstance(batches,list):
        batches = [batches]

    octree = Octree(Batch.Optimize(batches))
    viewer = Viewer(octree)
    viewer.Run()


def sectionize(atom):

    g = atom.graph
    all_cells = CAT([CELLSPERLEVEL(g)(i) for i in xrange(0,4)])

    to_display = []
    for cell in all_cells:
        centroid = CENTROID(g)(cell)
        if centroid[1] > 0.0:
            to_display.append(cell)

    return hccvolume(atom, to_display, [1.2, 1.2, 1.2])


def hccvolume(atom, to_display, expl=[1,1,1]):

    g = atom.graph
    coords = atom.coords
    sf = atom.scale
    color = atom.color

    def offset(point,expl=[1,1,1]):
        scaledpoint = [point[k]*expl[k] for k in range(3)]
        vect = VECTDIFF([scaledpoint,point])
        return vect

    def cells(batches,cellpoints,expl=[1,1,1]):
        for points in cellpoints:
            n = len(points)
            center = [coord/float(n) for coord in VECTSUM(points)]
            vect = offset(center,expl)
            points = [[point[k]+vect[k] for k in range(3)] for point in points]
            cell = MKPOL([points,[range(1,n+1)],None])
            cellBatch = Plasm.getBatches(cell)
            cellBatch[0].diffuse = color
            batches += cellBatch
            # view rotation
            rot = ROTN([ ACOS(INNERPROD([ [1,1,1],[0,0,1] ])), VECTPROD([ [1,1,1],[0,0,1] ]) ])
            batches += Plasm.getBatches(STRUCT([rot, MK([1,1,1])]))
        return batches

    m = g.getPointDim()
    d = g.getMaxDimCells()
    mapping = [dict(zip(CELLSPERLEVEL(g)(k),range(len(CELLSPERLEVEL(g)(k)))))
               for k in range(d+1)]

    chain = to_display
    chains = [[],[],[],[]]
    [chains[g.Level(node)].append(node) for node in chain[::-1]]

    nodepoints = [[g.getVecf(node)[i] for i in range(1,m+1)]
                  if m>2 else
                  [g.getVecf(node)[i] for i in range(1,m+1)]+[0.0]
                  for node in CELLSPERLEVEL(g)(0)]

    def translate(pointIds):

        return [nodepoints[mapping[0][vert[1]]] for vert in
                      enumerate(pointIds)]


    if chains[2]:
        facesAsEdges = [DOWNCELLS(g)(face) for face in chains[2]]
        facesAsVerts = [list(set(CAT(AA(DOWNCELLS(g))(face))))
                        for face in facesAsEdges]
        facepoints = AA(translate)(facesAsVerts)
    if d == 3:
        if chains[3]:
            solidsAsVerts = [UPCELLS(g)(cell) for cell in chains[3]]
            cellpoints = AA(translate)(solidsAsVerts)

    #this is the list of batches you want to display
    batches = []
    batches = cells(batches,cellpoints,expl)


    #Last transformations
    for b in batches:
        b.matrix = Mat4f.translate(*coords) * Mat4f.scale(sf, sf, sf) * b.matrix

    return batches


if __name__ == '__main__':
    
    peptides = ["ARG", "ASN","ASP","CYS", "GLN", "GLU", "GLY", 
                "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
                "SER", "THR", "TRP", "TYR", "VAL"]
    
    #peptides = ["ALA", "ARG", "ASN"]
    polipeps = mkpolypeps(peptides)

    view_compound(polipeps)