import os
import benoitmorel.scripts.experiments as exp
import benoitmorel.scripts.generax.launch_generax as  launch_generax
import benoitmorel.tools.families.generate_families_with_simphy as generate_families_with_simphy
import benoitmorel.tools.families.fam as fam
import knirschl.launch_raxmlng as launch_raxmlng
import knirschl.launch_fastme as launch_fastme


#parameters = generate_families_with_simphy.SimphyParameters()
#generate_families_with_simphy.generate_from_parameters(parameters, exp.families_datasets_root)


root_output = exp.families_datasets_root # output/families/
numer_of_families = 100

# -- simulation --
datadir = generate_families_with_simphy.generate_simphy(
    # tag, species, families, sites, model, bl_factor
    "test", 20, numer_of_families, 200, "GTR", 1.0,
    # dup_rate, loss_rate, transfer_rate, gene_conversion_rate
    0.0, 0.0, 1.0, 0.0,
    # perturbation, population, miss_species, miss_fam, root_output, seed
    perturbation, 10, 0.0, 0.0, root_output, 42,
    # cores
    1)

# -- infer gene trees --
launch_raxmlng.run_raxml_all(fam.get_families_dir(datadir), "GTR+G")

output = "misc/fastme-test.newick"
launch_fastme.run_fastme_all(fam.get_families_dir(datadir), output)
# maybe use benoitmorel/tools/families/run_fastme.py ?

# create families file
generax_families_file = os.path.join(datadir, "families.txt")
launch_generax.build_generax_families_file(datadir, "", "GTR+G", generax_families_file)
resultsdir = os.path.join(datadir, "misc")
launch_generax.run_generax(datadir, None, "SPR", "true", generax_families_file, "false", 1, "-r UndatedDL", resultsdir)