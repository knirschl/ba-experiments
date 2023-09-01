import subprocess
import sys
from ete3 import Tree
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/trees')
import paths
import rescale_bl as ete3conv


def ete3_rf_files(tree1, tree2, unrooted = True):
    return ete3_rf(ete3conv.read_tree(tree1), ete3conv.read_tree(tree2))[:2]

def ete3_rf(tree1, tree2, unrooted = True):
    if (len(tree2.children) == 3):
        tree2.set_outgroup(tree2.children[0])
    if (len(tree1.children) == 3):
        tree1.set_outgroup(tree1.children[0])
    res = tree1.robinson_foulds(tree2, unrooted_trees=unrooted, skip_large_polytomies = True, correct_by_polytomy_size = True) # may throw TreeError
    return [res[0], res[0] / (len(res[3]) + len(res[4]) - 2 * len(res[2]) - 2) if res[0] != 0.0 else res[0], res[1]]
    
def raxmlng_rf(tree1, tree2):
    command = []
    command.append(paths.raxml_exec)
    command.append("--rf")
    command.append(tree1 + "," + tree2)
    out = subprocess.check_output(command).decode("utf-8") # may throw CalledProcessError
    lines = out.split("\n")
    rf_abs = lines[0].split(" ")[-1]
    rf_rel = lines[1].split(" ")[-1]
    #print(rf_abs, rf_rel)
    return [int(float(rf_abs)), float(rf_rel)]
