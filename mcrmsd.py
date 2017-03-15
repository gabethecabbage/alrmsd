#!/usr/bin/python3

import Bio.PDB
import Bio.pairwise2
import numpy as np
import sys, getopt

import alrmsd

def main(argv):
    # Declaring default arg states
    ref_file = ''
    search_file = ''
    try:
        opts, args = getopt.getopt(argv,"havcr:s:",["ref_file=","search_file="])
    except getopt.GetoptError:
        #help_print()
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            #help_print()
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
    search_struc = pdb_parser.get_structure("search_struc", search_file)
    ref_struc = pdb_parser.get_structure("ref_struc", ref_file)
    search_coords_d = alrmsd.get_multichain_coords_d(search_struc)
    ref_coords_d = alrmsd.get_multichain_coords_d(ref_struc)
    chain_pairs = {}
    for sc_id in sorted(search_coords_d):
        best_rmsd = float("inf")
        best_chain = ""
        for rc_id in sorted(ref_coords_d):
            c_rmsd = alrmsd.calc_rmsd(ref_coords_d[rc_id], search_coords_d[sc_id])
            if c_rmsd < best_rmsd:
                best_rmsd = c_rmsd
                best_chain = rc_id
        chain_pairs[sc_id] = best_chain

    search_coords_assembled = []
    ref_coords_assembled = []
    for sc_id in sorted(chain_pairs):
        search_coords_assembled.append(search_coords_d[sc_id])
        ref_coords_assembled.append(ref_coords_d[chain_pairs[sc_id]])
    
    search_coords_assembled = np.concatenate(search_coords_assembled)
    ref_coords_assembled = np.concatenate(ref_coords_assembled)
    
    print(alrmsd.calc_rmsd(ref_coords_assembled, search_coords_assembled))
       

global verbose
verbose = False
global ca_only
ca_only = False
align_mode = False

if __name__ == "__main__":
   main(sys.argv[1:])


