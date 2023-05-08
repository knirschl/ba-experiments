import subprocess
import os
import benoitmorel.tools.families.fam as fam
import benoitmorel.scripts.experiments as exp
import knirschl.msa.fasta_to_phy as convert

def run_fastme_all(families_dir, output_file):
    for family in os.listdir(families_dir):
        # convert alignment.msa to alignment.phylip
        input_fasta = fam.get_alignment(families_dir, family)
        output_phy = input_fasta[:input_fasta.rfind(".") + 1] + "phy" # changes file ending
        convert.fasta_to_phy_file(input_fasta, output_phy) 
        run_fastme(output_phy, os.path.join(families_dir, family, output_file))

# fastme -i alignment.phylip -d -o <output-file>
def run_fastme(alignment, output_file):
    commands = []
    commands.append(exp.fastme_exec)
    commands.append("-i")
    commands.append(alignment)
    commands.append("-d") # DNA data
    commands.append("-o")
    commands.append(output_file)
    subprocess.check_call(commands)