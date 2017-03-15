#!/usr/bin/python3

import Bio.PDB
import Bio.pairwise2
import numpy as np
import sys, getopt

def help_print():
      print('Multifunction RMSD calculator utility', '\r\n',
              '##################################', '\r\n',
              'This tool takes in a refence pdb and a search object,', '\r\n',
              'and calculates rmsd for each reference chain against', '\r\n',
              'the search object chain.')
      print('usage: alrmsd.py -r <reference-pdb> -s <search-pdb>')


def atoms_2_coords(atoms, coords=[]):
    for atom in atoms:
        if ca_only == False or atom.id == 'CA':
            # Append all atoms if ca_only is false
            # or just carbon alpha otherwise
            coords.append(atom.coord)
    return coords

def get_coords(chain):
    '''
    Create NumPy Array of Atom Coords from Chain.
    In: Bio.PDB chain generator, Bool of which atoms to return
    Return: numpy array of all/ca atom coordinates
    '''   
    atoms = chain.get_atoms()
    coords_l = atoms_2_coords(atoms, coords=[])
    coords_np = np.asarray(coords_l)
    return coords_np

def get_aligned_coords(chain, align):
    coords = []
    rn = 0
    residues = chain.get_residues()
    for res in residues:
        if align[rn] != '-':
             coords = atoms_2_coords(res, coords)
        rn += 1
    coords = np.asarray(coords)
    return coords

def get_chain_seq(chain):
    pep_codes = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    seq_1 = "" 
    seq_3 = []
    residues = chain.get_residues()
    for residue in residues: 
        try: 
            seq_1 = seq_1 + pep_codes[residue.get_resname()]
        except:
            seq_1 = seq_1 + 'x'
        seq_3.append(residue.get_resname())
    return seq_1       

    ref_struc = pdb_parser.get_structure("ref_struc", ref_file)
    c_rmsd = {} #dictionary of chain ids to rmsds

def get_multichain_coords_d(struc, search_seq=None):
    coords_d = {}
    for chain in struc.get_chains():
        #iterate over the refence chains
        c_id = chain.get_id()
        c_ob = struc[0][c_id]

        if search_seq is None:
            #grab whole chain coords, calculate rmsd and store in dict
            coords_d[c_id] = get_coords(c_ob)

        else:
            #get the ref chain sequence and generate alignments
            chain_seq = get_chain_seq(c_ob)
            #2 points for match, -100 mismatch penalty, -0.5 gap open penalty, 0 gap extension penalty
            alignments = Bio.pairwise2.align.globalms(chain_seq, search_seq, 2, -100, -.5, 0)

            if '-' in alignments[0][0]:
                #if the reference chain has gaps following alignments, disregard
                coords_d[c_id] = None
                if verbose == True: print("Reference Chain", c_id, "not gapplesly aligned to search chain",
                        file=sys.stderr)
                if verbose == True: print(Bio.pairwise2.format_alignment(*alignments[0]),
                        file=sys.stderr)
            else:
                #grab the coords of ref residues that align to the search
                coords_d[c_id] = get_aligned_coords(c_ob, alignments[0][1])
                if verbose == True: print(Bio.pairwise2.format_alignment(*alignments[0]))
    return coords_d

def calc_rmsd(ref_coords, search_coords):
    '''
    Calculate Root Mean Squared Deviation
    In: Numpy Array of Reference and Search structure atom coords
    Out: Float RMSD value
    '''
    if ref_coords is None or search_coords is None:
        return float("inf")
    elif search_coords.shape != ref_coords.shape:
        print("Atom ranges for reference and search chains don't match", "\r\n",
                "Ref:", ref_coords.shape,
                "Tar:", search_coords.shape,
                file=sys.stderr)
        return float("inf")
    else:
        dist = search_coords - ref_coords
        l = ref_coords.shape[0]
        rms = np.sqrt(np.sum(dist**2) / l)
        return rms

def main(argv):
    # Declaring default arg states
    ref_file = ''
    search_file = ''
    global verbose
    verbose = False
    global ca_only
    ca_only = False
    global align_mode
    align_mode = False
    try:
        opts, args = getopt.getopt(argv,"havcr:s:",["ref_file=","search_file="])
    except getopt.GetoptError:
        help_print()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            help_print()
            sys.exit()
        elif opt in ("-v", "--verbose"):
            verbose = True
        elif opt in ("-c", "--carbonalpha"):
            ca_only = True
        elif opt in ("-a", "--seq-align"):
            align_mode = True
        elif opt in ("-r", "--reference"):
            ref_file = arg
        elif opt in ("-s", "--search"):
            search_file = arg
    pdb_parser = Bio.PDB.PDBParser(QUIET = True)

    #### LOAD SEARCH STRUCTURE ####
    search_struc = pdb_parser.get_structure("search_struc", search_file)
 
    search_chains = []
    for chain in search_struc.get_chains():
        search_chains.append(chain.get_id())

    #Get first chain in the first search structure model
    search_chain = search_struc[0][search_chains[0]]
    search_coords = get_coords(search_chain)
    
    if align_mode == True:
        #get the search chain sequence
        search_seq = get_chain_seq(search_chain)

    if verbose == True:
        print('Search PDB file:', search_file)
        print('Search chain:', search_chains[0])
        print('Reference PDB file:', ref_file)

    #### LOAD REFERENCE STRUCTURE ####
    ref_struc = pdb_parser.get_structure("ref_struc", ref_file)
    c_rmsd = {} #dictionary of chain ids to rmsds

    # Get sequence aligned/whole reference structure coords in a dictionary, per chain
    if align_mode == True:
        ref_coords_d = get_multichain_coords_d(ref_struc, search_seq)
    else:
        ref_coords_d = get_multichain_coords_d(ref_struc)

    for c_id in sorted(ref_coords_d):
        c_rmsd[c_id] = calc_rmsd(ref_coords_d[c_id], search_coords)
        if verbose == True: print("Reference Chain", c_id, "RMSD: \t", c_rmsd[c_id])

    best_chain = min(c_rmsd, key=c_rmsd.get)
    if verbose == True: print("Chain \t RMSD")
    print(best_chain, "\t", c_rmsd[best_chain])

global verbose
verbose = False
global ca_only
ca_only = False
global align_mode
align_mode = False
if __name__ == "__main__":
   main(sys.argv[1:])
