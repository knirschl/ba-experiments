import sys
import os
import ete3

def fasta_to_phy_file(input_fasta, output_phy):
    lines = open(input_fasta).readlines()
    fasta_str = ""
    for line in lines:
      if ("_" in line):
        fasta_str += (line[:-1] + " \n")
      else:
        fasta_str += line 
    msa = ete3.SeqGroup(fasta_str, format="fasta")
    msa.write("phylip", output_phy)

def fasta_to_phy_dir(input_fasta_dir, output_phy_dir):
  try:
    os.makedirs(output_phy_dir)
  except:
    print("WARNING: directory " + output_phy_dir + " already exists")
  for input_fasta_base in os.listdir(input_fasta_dir):
    input_fasta = os.path.join(input_fasta_dir, input_fasta_base)
    output_phy_base = ".".join(input_fasta_base.split(".")[:-1]) + ".phy"
    output_phy = os.path.join(output_phy_dir, output_phy_base)
    fasta_to_phy_file(input_fasta, output_phy)(input_fasta, output_phy)


if (__name__ == "__main__"):
  if (len(sys.argv) != 3):
    print("Syntax: python fasta_to_phy.py input_fasta_dir output_phy_dir")
    exit(1)
  fasta_to_phy_dir(sys.argv[1], sys.argv[2])
