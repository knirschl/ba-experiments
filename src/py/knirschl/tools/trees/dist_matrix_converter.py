import os
import itertools
import sys
from phylodm import PhyloDM
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
sys.path.insert(0, 'tools/families')
import fam

def from_newick(file, norm = True):
    pdm = PhyloDM.load_from_newick_path(file)
    dist_matrix = pdm.dm(norm)
    labels = pdm.taxa()
    max_len = -1
    for r in dist_matrix:
        for e in r:
            if (len(str(e)) > max_len):
                max_len = len(str(e))
    with open(os.path.join(os.path.dirname(file), "speciesTree" + ".phy"), "w") as writer:
        writer.write(str(len(labels)) + '\n')
        for i in range(len(labels)):
            writer.write(labels[i] + " " * (10 - len(labels[i])))
            for e in dist_matrix[i]:
                writer.write(str(e))
                diff = max_len - len(str(e)) 
                if (diff > 0):
                    writer.write("0" * diff)
                writer.write('\t')
            writer.write('\n')

def format_phylip(dm, handle):
    """
    Slightly modified version from Bio.Phylo.TreeConstruction.py to fit the strict phylip format
    """
    handle.write(f"{len(dm.names)}\n")
    # Phylip needs space-separated, vertically aligned columns
    name_width = max(10, max(map(len, dm.names)) + 1)
    value_fmts = ("{" + str(x) + ":.19f}" for x in range(1, len(dm.matrix) + 1))
    row_fmt = "{0:" + str(name_width) + "s}" + "  ".join(value_fmts) + "\n"
    for i, (name, values) in enumerate(zip(dm.names, dm.matrix)):
        # Mirror the matrix values across the diagonal
        mirror_values = (dm.matrix[j][i] for j in range(i + 1, len(dm.matrix)))
        fields = itertools.chain([name], values, mirror_values)
        handle.write(row_fmt.format(*fields))

def from_fasta(file):
    msa = AlignIO.read(file, 'fasta')
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(msa)
    format_phylip(dm, open(os.path.join(file + ".matrix" + ".phy"), "w"))

def convert_input(datadir):
    from_newick(fam.get_true_species_tree(datadir))
    for family in fam.get_families_list(datadir):
        from_fasta(fam.get_alignment(datadir, family))