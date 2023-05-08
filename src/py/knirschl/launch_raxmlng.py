import subprocess
import os
import benoitmorel.tools.families.fam as fam
import benoitmorel.scripts.experiments as exp

def run_raxml_all(families_dir, model):
    # number = number_of_families = max(map(lambda s: int(s.split("_")[1]), os.listdir(families_dir)))
    for family in os.listdir(families_dir):
        run_raxml(fam.get_alignment(families_dir, family), model)

# raxml-ng --msa output/families/ssim_.../families/family_XXX/alignment.msa --model GTR+G
def run_raxml(alignment, model):
    commands = []
    commands.append(exp.raxml_exec_no_mpi) # no_mpi version ?
    commands.append("--msa")
    commands.append(alignment)
    commands.append("--model")
    commands.append(model)
    subprocess.check_call(commands)
