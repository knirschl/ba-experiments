import subprocess
import os
import benoitmorel.tools.families.fam as fam
import benoitmorel.scripts.experiments as exp
import knirschl.msa.fasta_to_phy as convert


# fastme -i alignment.phylip -d -o <output-file>
def get_fastme_command(alignment, output_file):
    command = []
    command.append(exp.fastme_exec)
    command.append("-i")
    command.append(alignment)
    command.append("-d") # DNA data
    command.append("-o")
    command.append(output_file)
    return command

def get_run_specs(command):
    specs = []
    specs.append("FastME")
    specs.append(time.strftime("%d.%m.%Y-%H:%M:%S"))
    specs.append(command)
    return specs

def run_fastme_all(datadir, output_file):
    families_dir = fam.get_families_dir(datadir)
    for family in os.listdir(families_dir):
        # convert alignment.msa to alignment.phylip
        input_fasta = fam.get_alignment(families_dir, family)
        output_phy = input_fasta[:input_fasta.rfind(".") + 1] + "phy" # changes file ending
        convert.fasta_to_phy_file(input_fasta, output_phy) 
        run_fastme(datadir, output_phy, os.path.join(families_dir, family, output_file))

def run_fastme(datadir, alignment, output_file):
    command = get_fastme_command(alignment, output_file)
    run_name = get_run_specs(command)
    start = time.time()
    subprocess.check_call(command)
    end = time.time() - start
    metrics.save_metrics(datadir, run_name, end, "runtimes")