import itertools
import numpy as np
import os
import sys
import time
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from phylodm import PhyloDM
sys.path.insert(0, os.path.join('scripts', 'programs'))
sys.path.insert(0, 'tools/families')
import launch_fastme
import fam

def get_key_order(keys):
    """
    [4, 2, 3, 1] -> [3, 1, 2, 0]
    """
    return np.argsort(keys)

def sparse_triu_to_sym(triu):
    # [[0], [1, 0]] -> [[0, 0], [1, 0]]
    triu =  np.array([np.pad(np.array(triu[i]), (0, len(triu) - len(triu[i])), 'constant') for i in range(len(triu))])
    # add upper triangle to lower 0s
    # `- np.diag(...)` not necessary in this case as biopython always sets diagonal to 0
    return triu + triu.T - np.diag(np.diag(triu))

def sort(matrix, order):
    # arr[order] reorders array, doing it twice reorders columns and rows
    return np.array([matrix[order][i][order] for i in range(len(order))])

def from_newick(file, norm = True):
    # calculate distance matrix
    pdm = PhyloDM.load_from_newick_path(file)
    # make labels same as in alignments (map function?)
    labels = pdm.taxa()
    # sort distance matrix and labels to add the correct values in later steps (standardized output)
    dist_matrix = sort(pdm.dm(norm), get_key_order(labels))
    labels = np.sort(labels)
    # write to file
    write_phylip(dist_matrix, labels, open(os.path.join(os.path.dirname(file), "speciesTree" + ".matrix" + ".phy"), "w"))

def from_fasta(file):
    # calculate distance matrix
    msa = AlignIO.read(file, 'fasta')
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(msa)
    # sort distance matrix and labels to add the correct values in later steps (standardized output)
    dist_matrix = sort(sparse_triu_to_sym(dm.matrix), get_key_order(dm.names))
    labels = np.sort(dm.names)
    # write to file
    write_phylip(dist_matrix, labels, open(os.path.join(file + ".matrix" + ".phy"), "w"))

def write_phylip(dist_matrix, labels, handle):
    """
    Modified version from Bio.Phylo.TreeConstruction.py to fit my needs
    (= strict phylip format and generalized to not only work with biopython)
    """
    handle.write(f"{len(labels)}\n")
    # Phylip needs space-separated, vertically aligned columns
    name_width = max(10, max(map(len, labels)) + 1)
    value_fmts = ("{" + str(x) + ":.19f}" for x in range(1, len(dist_matrix) + 1))
    row_fmt = "{0:" + str(name_width) + "s}" + "  ".join(value_fmts) + "\n"
    for i, (label, values) in enumerate(zip(labels, dist_matrix)):
        fields = itertools.chain([label], values)
        handle.write(row_fmt.format(*fields))

def convert_input(datadir, cores):
    # species tree
    start = time.time()
    from_newick(fam.get_true_species_tree(datadir))
    print("Converting the species trees takes {}s".format(time.time() - start))
    # gene alignments
    start = time.time()
    launch_fastme.run_fastme_matrix(datadir, cores=cores)
    print("Converting all gene trees takes {}s".format(time.time() - start))


""" # TIME TEST
datadir = "/home/fili/Desktop/2023/BA/code/output/families/ssim_DL_s100_f100_sites200_GTR_bl1.0_d1.0_l1.0_t0.0_gc0.0_p0.0_pop10_ms0.0_mf0.0_seed9724682"
s1 = time.time()
for family in fam.get_families_list(datadir):
    from_fasta(fam.get_alignment(datadir, family))
t1 = time.time() - s1
s2 = time.time()
for family in fam.get_families_list(datadir):
    launch_fastme.calculate_dist_matrix(datadir, family, True, os.path.join(fam.get_family_path(datadir, family), "alignment.msa.matrix.phy"))
t2 = time.time() - s2
s3 = time.time()
launch_fastme.run_fastme_matrix(datadir, cores=8)
t3 = time.time() - s3
print("AlignIO takes {}s ({}s per tree)".format(t1, t1 / 100))
print("FastME-single takes {}s ({}s per tree)".format(t2, t2 / 100))
print("FastME-sched takes {}s ({}s per tree)".format(t3, t3 / 100)) """
