import subprocess
import sys
from ete3 import Tree
sys.path.insert(0, 'scripts')
import paths


def ete3_rf(tree1, tree2, unrooted = True):
    if (len(tree2.children) == 3):
        tree2.set_outgroup(tree2.children[0])
    if (len(tree1.children) == 3):
        tree1.set_outgroup(tree1.children[0])
    return tree1.robinson_foulds(tree2, unrooted_trees=unrooted, skip_large_polytomies = True, correct_by_polytomy_size = True)

def raxmlng_rf(tree1, tree2):
    command = []
    command.append(paths.raxml_exec)
    command.append("--rf")
    command.append(tree1 + "," + tree2)
    try:
        out = subprocess.check_output(command).decode("utf-8")
    except subprocess.CalledProcessError as error:
        return ["abort", error]
    lines = out.split("\n")
    rf_abs = lines[0].split(" ")[-1]
    rf_rel = lines[1].split(" ")[-1]
    #print(rf_abs, rf_rel)
    return [float(rf_abs), float(rf_rel)]
