import subprocess
import os
import sys
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
import paths
import fam

# raxml-ng --msa output/families/ssim_.../families/family_XXX/alignment.msa --model GTR+G
def get_raxml_command(alignment, model, with_mpi = False, cores = 1):
    command = []
    if (with_mpi):
        command.append(paths.raxml_mpi_exec)
        command.append("-np")
        command.append(str(cores))
    else:
        command.append(paths.raxml_exec)
    command.append("--msa")
    command.append(alignment)
    command.append("--model")
    command.append(model)
    return command

def get_run_specs(command):
    specs = []
    specs.append("RAxML-NG")
    specs.append(time.strftime("%d.%m.%Y-%H:%M:%S"))
    specs.append(command)
    return specs

def run_raxml_all(datadir, model, with_mpi = False, cores = 1):
    families_dir = fam.get_families_dir(datadir)
    for family in os.listdir(families_dir):
        run_raxml(datadir, fam.get_alignment(families_dir, family), model, cores, with_mpi)

def run_raxml(datadir, alignment, model, with_mpi = False, cores = 1):
    command = get_raxml_command(alignment, model, with_mpi, cores)
    run_name = get_run_specs(command)
    start = time.time()
    #subprocess.check_call(command)
    end = time.time() - start
    #saved_metrics.save_metrics
    #metrics.save_metrics(datadir, run_name, end, "runtimes")
    