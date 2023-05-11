import math
import os
import sys
import subprocess
import shutil
import copy
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/trees')
print(sys.path)
import paths
import utils
import fam
import rescale_bl
import analyze_tree
import discordance_rate
import sample_missing_data

class SimphyParameters():
  def __init__(self, tag = "dtl", prefix = "ssim", speciations_per_year = 0.000000005,
              extinction_per_year = 0.0000000049, species_taxa = 20 ,
              families_number = 100, bl = 1.0, loss_rate = 0.0, dup_rate = 0.0,
              transfer_rate = 1.0, gene_conversion_rate = 0.0, sites = 200,
              model = "GTR", seed = 42, distance_hgt = False, population = 10,
              miss_species = 0.0, miss_fam = 0.0):
    self.tag = tag
    self.prefix = prefix
    self.speciations_per_year = speciations_per_year
    self.extinction_per_year = extinction_per_year
    self.species_taxa = species_taxa
    self.families_number = families_number
    self.bl = bl     
    self.loss_rate = loss_rate 
    self.dup_rate = dup_rate
    self.transfer_rate = transfer_rate
    self.gene_conversion_rate = gene_conversion_rate
    self.sites = sites
    self.model = model
    self.seed = seed
    self.distance_hgt = distance_hgt
    self.population = population
    self.miss_species = miss_species 
    self.miss_fam = miss_fam

def build_config_file(parameters, output_dir):
  config_file = os.path.join(output_dir, "simphy_config.txt")
  with open(config_file, "w") as writer:
    writer.write("// SPECIES TREE\n")
    # number of replicates
    writer.write("-RS 1 // number of replicates\n")
    # speciation rates (speciation per yer)
    writer.write("-sb f:" + str(parameters.speciations_per_year) + "\n")
    writer.write("-sd f:" + str(parameters.extinction_per_year) + "\n")
    # number of species taxa 
    writer.write("-sl f:" + str(parameters.species_taxa) + "\n")
    # species tree height in years (I don't understand this)
    writer.write("-st ln:21.25,0.2\n")
    # substitution rate
    #writer.write("-su ln:-21.9," + str(0.1 * parameters.bl) + "\n")
    writer.write("-su ln:-21.9,0.1\n")
    # L, D, T global rates 
    lognormal_scale = 1.0
    lognormal_location = 0.0 #math.log(1.0 - 0.5 * pow(lognormal_scale, 2.0))
    #lognormal_mean = math.exp((lognormal_location + pow(lognormal_scale, 2.0)) / 2.0)
    loss_freq =     0.00000000049 * parameters.loss_rate 
    dup_freq =      0.00000000049 * parameters.dup_rate
    transfer_freq = 0.00000000049 * parameters.transfer_rate 
    gene_conversion_freq = 0.00000000049 * parameters.gene_conversion_rate 

    assert(loss_freq == dup_freq)
    writer.write("-gd f:" + str(loss_freq) + "\n")
    writer.write("-gb f:" + str(dup_freq) + "\n")
    writer.write("-gt f:" + str(transfer_freq) +"\n")
    writer.write("-gg f:" + str(gene_conversion_freq) +"\n")

    # L, D, T per family rates
    writer.write("-ld sl:" + str(lognormal_location) + "," + str(lognormal_scale) + ",gd\n")
    writer.write("-lb f:ld\n")
    writer.write("-lt sl:" + str(lognormal_location) + "," + str(lognormal_scale) + ",gt\n")
    writer.write("-lg f:gg\n")
    #writer.write("-lt f:ld\n")
    
    lk = 0
    if (parameters.distance_hgt):
        lk = 1
    writer.write("-lk " + str(lk) + "\n")
    
    #writer.write("-lb sl:" + str(lognormal_location) + "," + str(lognormal_scale) + ",gb\n")
    #writer.write("-lt sl:" + str(lognormal_location) + "," + str(lognormal_scale) + ",gt\n")

    writer.write("// POPULATION\n")
    writer.write("-SP f:" + str(parameters.population) + "\n")

    writer.write("// LOCUS\n")
    writer.write("-rl f:" + str(parameters.families_number) + " // locus (gene family) per replicate\n")

    writer.write("// Subsitution rates heterogeneity parameters\n")
    writer.write("-hs ln:1.5,1\n")
    writer.write("-hl ln:1.551533,0.6931472\n")
    writer.write("-hg ln:1.5,1\n")


    writer.write("// GENERAL\n")
    writer.write("-cs " + str(parameters.seed) + "\n") 
    writer.write("-O " + str(output_dir) + " // output directory\n")
    writer.write("-OM 1 // output the mappings\n")
    writer.write("-OC 1 // log the configuration file\n")
    writer.write("-OD 1 // log the configuration file\n")
    writer.write("-OP 1 // log the configuration file\n")

  return config_file

def build_indelible_config_file(parameters, output_dir):
  config_file = os.path.join(output_dir, "indelible_config.txt")
  with open(config_file, "w") as writer:
    
    sites_mean = parameters.sites
    sites_teta = 0.25
    sites_mu = 0
    sites_scale = sites_mean / math.exp(sites_mu + sites_teta * sites_teta / 2.0)
    #sites_min = str(20)
    #sites_max = str(2 * int(parameters.sites) - 20)
    writer.write("[TYPE] NUCLEOTIDE 1\n") # DNA using algorithm 1 
    writer.write("[SETTINGS] [fastaextension] fasta\n")
    writer.write("[SIMPHY-UNLINKED-MODEL] modelA \n")
    if ("GTR" == parameters.model):
      writer.write("  [submodel] GTR $(rd:16,3,5,5,6,15) // GTR with rates from a Dirichlet  \n")
      writer.write("  [statefreq] $(d:36,26,28,32)  // frequencies for T C A G sampled from a Dirichlet \n")
      writer.write("[rates] 0 $(e:2) 0 // Site-specific rate heterogeneities: 0 p-inv, alpha from an E(2) and using a continuous gamma distribution.\n")
    else:
      assert(False)

    #writer.write("[SIMPHY-PARTITIONS] simple [1.0 modelA $(u::" + sites_min + "," + sites_max + ")]\n")
    writer.write("[SIMPHY-PARTITIONS] simple [1.0 modelA $(sl:" + str(sites_mu) + "," + str(sites_teta) + "," + str(sites_scale) + ")]\n")

    writer.write("[SIMPHY-EVOLVE] 1 dataset \n")

  return config_file

def build_mapping(simphy_mapping, phyldog_mapping):
  lines = open(simphy_mapping).readlines()[1:]
  dico = {}
  for line in lines:
    if (not line.startswith("'")):
      continue
    split = line.replace("'", "").split("\t")
    gene = split[0]
    species = gene.split("_")[0]
    if (not species in dico):
      dico[species] = []
    dico[species].append(gene)
  with open(phyldog_mapping, "w") as phyldog_writer:
    for species, genes in dico.items():
      phyldog_writer.write(species + ":" + ";".join(genes) + "\n")

def run_simphy(output_dir, config_file):
  commands = []
  commands.append(paths.simphy_exec)
  commands.append("-I")
  commands.append(config_file)
  subprocess.check_call(commands)

  os.remove(config_file)

  simphy_output_dir = os.path.join(output_dir, "1")
  species_tree = os.path.join(simphy_output_dir, "s_tree.trees")
  generations = analyze_tree.check_ultrametric_and_get_length(species_tree)

def rescale_gene_tree_bl(output_dir, bl):
  simphy_output_dir = os.path.join(output_dir, "1")
  bl = float(bl)
  print("BEFORE RESCALE")
  if (bl == 1.0):
    return
  families = []
  print("outputdir " + simphy_output_dir)
  for f in os.listdir(simphy_output_dir):
    if (f.startswith("g_trees")):
      families.append("family_" + f.split("g_trees")[1].split(".")[0])
  for family in families:
    family_number = family.split("_")[1] 
    # true trees
    gene_tree = os.path.join(simphy_output_dir, "g_trees" + family_number + ".trees")
    rescale_bl.rescale_bl(gene_tree, gene_tree, bl)

def run_indelible(output_dir, config_file, cores):
  commands = []
  seed = "42"
  commands.append("perl")
  commands.append(paths.simphy_indelible_wrapper_exec)
  commands.append(output_dir)
  commands.append(config_file)
  commands.append(seed)
  commands.append(str(cores))
  subprocess.check_call(commands)

def copy_trim(input_file, output_file):
  s = open(input_file).read()
  open(output_file, "w").write(s.replace(" ", ""))

def export_to_family(output_dir, replicate = 1):
  print("Start exporting to families format...")
  fam.init_top_directories(output_dir)
  simphy_output_dir = os.path.join(output_dir, str(replicate))
  families = []
  for f in os.listdir(simphy_output_dir):
    if (f.startswith("g_trees")):
      families.append("family_" + f.split("g_trees")[1].split(".")[0])
  fam.init_families_directories(output_dir, families)
  # species tree
  species = os.path.join(simphy_output_dir, "s_tree.trees")
  shutil.copyfile(species, fam.get_species_tree(output_dir))
  for family in families:
    family_number = family.split("_")[1] 
    # true trees
    alignment = os.path.join(output_dir, "1", "dataset_" + family_number + ".fasta")
    
    # check that the alignment contain all characters (otherwise phyldog crashes)
    alignment_content = open(alignment).read()
    ok = True
    for c in ['A', 'C', 'G', 'T']:
      if (not c in alignment_content):
        ok = False
    if (not ok):
      shutil.rmtree(fam.get_family_path(output_dir, family))
      print("rm " + fam.get_family_path(output_dir, family))
      continue

    gene_tree = os.path.join(simphy_output_dir, "g_trees" + family_number + ".trees")
    shutil.copy(gene_tree, fam.get_true_tree(output_dir, family))
    
    simphy_mapping = os.path.join(simphy_output_dir, family_number + "l1g.maplg")
    phyldog_mapping = fam.get_mappings(output_dir, family)
    build_mapping(simphy_mapping,  phyldog_mapping)
    # alignment
    # true trees
    # alignment
    copy_trim(alignment, fam.get_alignment(output_dir, family))
    #copy_and_rename_alignment(alignment, fam.get_alignment(out, family), family)
  fam.postprocess_datadir(output_dir)

def get_output_dir(parameters, root_output):
  res = parameters.prefix
  res += "_" + parameters.tag
  res += "_s" + str(parameters.species_taxa)
  res += "_f" + str(parameters.families_number)
  res += "_sites" + str(parameters.sites)
  res += "_" + str(parameters.model).replace("+", "")
  res += "_bl" + str(parameters.bl)
  res += "_d" + str(parameters.dup_rate)
  res += "_l" + str(parameters.loss_rate)
  res += "_t" + str(parameters.transfer_rate)
  res += "_gc" + str(parameters.gene_conversion_rate)
  res += "_p0.0"
  res += "_pop" + str(parameters.population)
  res += "_ms" + str(parameters.miss_species)
  res += "_mf" + str(parameters.miss_fam)
  """
  res += "_hgt" 
  if (parameters.distance_hgt):
    res += "dist"
  else:
    res += "unif"
  """
  res += "_seed" + str(parameters.seed)
  return os.path.join(root_output, res)

def compute_and_write_discordance_rate(parameters, output_dir):
  d = 0.0
  if (int(parameters.population) > 20):
      if (parameters.dup_rate == 0.0 and parameters.transfer_rate == 0.0):
        d = discordance_rate.get_discordance_rate(output_dir)
      else:
        no_dtl_parameters = copy.deepcopy(parameters) 
        no_dtl_parameters.dup_rate = 0.0
        no_dtl_parameters.loss_rate = 0.0
        no_dtl_parameters.transfer_rate = 0.0
        temp_output_dir = generate_from_parameters(no_dtl_parameters, output_dir)
        d = fam.get_discordance_rate(temp_output_dir)
        shutil.rmtree(temp_output_dir)
    
  print("Discordance rate: " + str(d))
  fam.write_discordance_rate(output_dir, d)

def generate_from_parameters(parameters, root_output, cores = 1):
  cores = cores
  output_dir = get_output_dir(parameters, root_output)
  utils.reset_dir(output_dir)
  config_file = build_config_file(parameters, output_dir)
  run_simphy(output_dir, config_file)
  rescale_gene_tree_bl(output_dir, parameters.bl)
  indelible_config_file = build_indelible_config_file(parameters, output_dir)
  run_indelible(output_dir, indelible_config_file, cores)
  export_to_family(output_dir)
  compute_and_write_discordance_rate(parameters, output_dir)
  print("Done! output in " + output_dir) 
  return output_dir

def generate_simphy(tag, species, families, sites, model, bl_factor, dup_rate, loss_rate, transfer_rate, gene_conversion_rate, population, miss_species, miss_fam, root_output, seed, cores):
  p = SimphyParameters()
  p.tag = tag
  p.species_taxa = int(species)
  p.families_number = int(families)
  p.sites = int(sites)
  p.model = model
  p.bl = float(bl_factor)
  p.dup_rate = float(dup_rate)
  p.loss_rate = float(loss_rate)
  p.transfer_rate = float(transfer_rate)
  p.gene_conversion_rate = float(gene_conversion_rate)
  p.seed = int(seed)
  p.population = population
  p.miss_species = float(miss_species)
  print(miss_fam)
  p.miss_fam = float(miss_fam)
  if (p.miss_species != 0.0 or p.miss_fam != 0.0):
    p.prefix = p.prefix + "temp"
    temp_output_dir = generate_from_parameters(p, root_output)
    p.prefix = p.prefix[:-4]
    output_dir = get_output_dir(p, root_output) 
    print("Now move " + temp_output_dir + " to " + output_dir)
    sample_missing_data.sample_missing_data(temp_output_dir, output_dir, p.miss_species, p.miss_fam)
    shutil.rmtree(temp_output_dir)
  else:
    return generate_from_parameters(p, root_output, cores)

def get_param_from_dataset_name(parameter, dataset):
  print(parameter)
  print(dataset)
  if (parameter == "tag"):
    return dataset.split("_")[1]
  elif (parameter == "species"):
    return dataset.split("_")[2][1:]
  elif (parameter == "families"):
    return dataset.split("_")[3][1:]
  elif (parameter == "sites"):
    return dataset.split("_")[4][5:]
  elif (parameter == "model"):
    return dataset.split("_")[5]
  elif (parameter == "bl"):
    return dataset.split("_")[6][2:]
  elif (parameter == "dup_rate"):
    return dataset.split("_")[7][1:]
  elif (parameter == "loss_rate"):
    return dataset.split("_")[8][1:]
  elif (parameter == "transfer_rate"):
    return dataset.split("_")[9][1:]
  elif (parameter == "gene_conversion_rate"):
    return dataset.split("_")[10][2:]
  elif (parameter == "perturbation"):
    return dataset.split("_")[11][1:]
  elif (parameter == "population"):
    return dataset.split("_")[12][3:]
  elif (parameter == "ms"):
    return dataset.split("_")[13][2:]
  elif (parameter == "msmf"):
    return dataset.split("_")[13][2:]
  elif (parameter == "mf"):
    return dataset.split("_")[14][2:]
  elif (parameter == "av_miss"):
    ms = float(get_param_from_dataset_name("ms", dataset))
    mf = float(get_param_from_dataset_name("mf", dataset))
    return str(round(ms * mf, 2))
  elif (parameter == "seed"):
    return dataset.split("_")[15][4:]
  elif (parameter == "tl_ratio"):
    t = get_param_from_dataset_name("transfer_rate", dataset)
    l = get_param_from_dataset_name("loss_rate", dataset)
    if (float(l) == 0.0):
      return "-1.0"
    return  str(float(t)/float(l))
  elif (parameter == "dl_ratio"):
    d = get_param_from_dataset_name("dup_rate", dataset)
    l = get_param_from_dataset_name("loss_rate", dataset)
    if (float(l) == 0.0):
      return "-1.0"
    return str(float(d)/float(l))
  elif (parameter == "dt_ratio"):
    d = get_param_from_dataset_name("dup_rate", dataset)
    t = get_param_from_dataset_name("transfer_rate", dataset)
    if (float(t) == 0.0):
      return "-1.0"
    return str(float(d)/float(t))
  elif (parameter == "av_rate"):
    d = float(get_param_from_dataset_name("dup_rate", dataset))
    l = float(get_param_from_dataset_name("loss_rate", dataset))
    t = float(get_param_from_dataset_name("transfer_rate", dataset))
    if (float(t) == 0.0):
      return "-1.0"
    return str((d + t + l) / 2.0)
  elif (parameter == "discordance"):
    return float(fam.get_discordance_rate(fam.get_datadir(dataset)))
  else:
    return "invalid"
   
def generate_simphy(dataset):
  if (not dataset.startswith("ssim")):
    print("Unknown simulator tag for dataset " + dataset)
    sys.exit(1)
  tag = get_param_from_dataset_name("tag", dataset)
  species = int(get_param_from_dataset_name("species", dataset))
  families = int(get_param_from_dataset_name("families", dataset))
  sites = int(get_param_from_dataset_name("sites", dataset))
  model = get_param_from_dataset_name("model", dataset)
  bl_factor = float(get_param_from_dataset_name("bl", dataset))
  d = float(get_param_from_dataset_name("dup_rate", dataset))
  l = float(get_param_from_dataset_name("loss_rate", dataset))
  t = float(get_param_from_dataset_name("transfer_rate", dataset))
  gc = float(get_param_from_dataset_name("gene_conversion_rate", dataset))
  p = float(get_param_from_dataset_name("perturbation", dataset))
  seed = int(get_param_from_dataset_name("seed", dataset))
  population = get_param_from_dataset_name("population", dataset)
  miss_species = get_param_from_dataset_name("ms", dataset)
  miss_fam = get_param_from_dataset_name("mf", dataset)
  generate_simphy(tag, species, families, sites, model, bl_factor, d, l, t, gc, p, population, miss_species, miss_fam, output,  seed) 
  
if (__name__ == "__main__"):
  parameters = SimphyParameters()
  generate_from_parameters(parameters, paths.families_datasets_root)
