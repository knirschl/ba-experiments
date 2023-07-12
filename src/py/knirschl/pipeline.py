import sys
import os
import time
import random
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'scripts/programs')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/simulation')
sys.path.insert(0, 'tools/trees')
import utils
import paths
import fam
import generate_with_simphy as simphy
import launch_raxml
import launch_generax
import launch_fastme
import launch_ba
import rescale_bl
import dist_matrix_converter
import compare_trees
import metrics

class RunFilter():
    def __init__(self, generate = True):
        self.disable_all()
        self.generate = True
        self.raxml = True
        self.generax = True
        self.fastme = True
        self.ba = True

        #self.force_overwrite = True
        self.compare = True
    
    def disable_all(self):
        self.generate = False
        self.raxml = False
        self.generax = False
        self.fastme = False
        self.ba = False
        self.force_overwrite = False
        self.compare = False
    
    def script_ba(self):
        self.generate = False
        self.raxml = False
        self.generax = False
        self.fastme = False

    def run_compare(self):
        self.script_ba()
        self.ba = False

    def run_methods(self, datadir, subst_model, cores):
        if (self.generate):
            if (not os.path.isdir(datadir) or self.force_overwrite):
                utils.printFlush("Run simphy...")
                dataset = os.path.basename(datadir)
                simphy.generate_dataset(dataset, cores)
        print("**************************************************************************")
        print("Run tested gene tree inference tools for dataset " + datadir)
        print("**************************************************************************")
        if (len(datadir.split("/")) == 1):
          datadir = fam.get_datadir(datadir) 
        save_stdout = sys.stdout
        #redirected_file = os.path.join(datadir, "runs", "logs_run_all_genes." + subst_model + ".txt")
        #print("Redirected logs to " + redirected_file)
        sys.stdout.flush()
        # RUN
        if(self.raxml):
            utils.printFlush("Run raxml-ng...")
            try:
                launch_raxml.run_raxmlng_on_families(datadir, subst_model, cores)
            except Exception as exc:
                utils.printFlush("Failed running RAxML-NG\n" + str(exc))
        if (self.generax):
            utils.printFlush("Run generax...")
            try:
                species_tree = fam.get_species_tree(datadir)
                resultsdir = os.path.join(datadir, "runs", subst_model)
                launch_generax.run(datadir, subst_model, "SPR", species_tree, "random", cores, "", resultsdir)#, do_analyze=False)
            except Exception as exc:
                utils.printFlush("Failed running GeneRax\n" + str(exc))
        if (self.fastme):
            utils.printFlush("Run fastme...")
            try:
                launch_fastme.run_fastme_on_families(datadir, subst_model, 1, 1)
            except Exception as exc:
                utils.printFlush("Failed running FastME\n" + str(exc))
        if (self.ba):
            utils.printFlush("Run ba...")
            try:
                dist_matrix_converter.convert_input(datadir)
                species_tree = fam.get_true_species_tree_matrix_sorted(datadir)
                launch_ba.run_ba_on_families(datadir, "exp", species_tree, cores)
            except Exception as exc:
                utils.printFlush("Failed running bachelor thesis script\n" + str(exc))
        # COMPARE INFERRED TREES WITH TRUE TREE
        if (self.compare):
            utils.printFlush("Run compare...")
            try:
                compare_trees.compare_all(datadir)
            except Exception as exc:
                utils.printFlush("Failed running compare\n" + str(exc))


# TOGGLE PIPELINE ELEMENTS
# ====== ! CAREFUL ! ======
run_filter = RunFilter() # all enabled
#run_filter.force_overwrite = True # regenerate old dataset
run_filter.raxml = False
#run_filter.script_ba() # only ba script
#run_filter.compare = False
#run_filter.run_compare() # only compare inferred trees
# ====== ! CAREFUL ! ======

root_output = paths.families_datasets_root # output/families/
#seeds = [42, 1007, 19732311, 121873, 14976684177860080345]
#seeds = [1007, 1058026, 1091512, 1105070, 1143740, 11872, 121873, 125546, 1316419, 1320646, 14976684177860080345, 1570525, 1674005, 19732311, 2181913, 2366262, 2453027, 2525498, 2530855, 2545197, 2622038, 2650353, 2835503, 2862791, 2967039, 2967641, 3174338, 3219882, 3279056, 3389159, 340025, 3547108, 3614207, 3686221, 3939144, 3959336, 4108907, 4130562, 4144379, 4193573, 42, 4229, 4290505, 429510, 4354376, 436633, 4412646, 4602564, 4754615, 4771372, 4785291, 4824437, 491587, 5476033, 558959, 5597160, 5616680, 5749959, 6050422, 6077656, 6214982, 6382595, 6449414, 6645706, 6722758, 6734655, 7060210, 7265233, 7340707, 7342556, 7350931, 7538361, 7724757, 7927143, 7940684, 8039292, 8124350, 8386866, 8388955, 842036, 8438226, 8527944, 8563774, 8670475, 8710215, 8928866, 901464, 9070912, 9128753, 9179296, 9387015, 9442585, 9476398, 9488374, 9495297, 9510809, 9718279, 9919929, 9930269, 9955696]
seeds = []
while (len(seeds) != 50):
    seeds.append(random.randrange(0, 9999999))
tag = "DL"
d = l = 1.0
s = 100 # def = 20
replicates = []

# Run multiple replicates
for seed in seeds:
    # SET simphy PARAMETERS 
    simphy_parameters = simphy.SimphyParameters(tag=tag, seed=seed, dup_rate=d, loss_rate=l, species_taxa=s)
    datadir = simphy.get_output_dir(simphy_parameters, root_output)
    replicates.append(datadir)
    print(datadir)

    # RUN PIPELINE
    rep_start = time.time()
    try:
        run_filter.run_methods(datadir, "F81", 8)
    finally:
        elapsed = time.time() - rep_start
        print("End of single experiment. Elapsed time: " + str(elapsed) + "s")
        metrics.save_metrics(datadir, "pipeline_" + tag + str(seed), elapsed, "runtimes")

# AVERAGE OVER ALL REPLICATES
abs_name = "rf_distance_avg-abs"
rel_name = "rf_distance_avg-rel"
abs_avgs_dico = metrics.get_metrics(replicates[0], abs_name)
rel_avgs_dico = metrics.get_metrics(replicates[0], rel_name)
rep_counter = 1
for rep in replicates[1:]:
    cur_abs = metrics.get_metrics(rep, abs_name)
    cur_rel = metrics.get_metrics(rep, rel_name)
#    abs_avgs_dico = {x: (float(abs_avgs_dico[x]) * rep_counter + float(cur_abs[x])) / (rep_counter + 1) for x in set(abs_avgs_dico).union(cur_abs)}
#    rel_avgs_dico = {x: (float(rel_avgs_dico[x]) * rep_counter + float(cur_rel[x])) / (rep_counter + 1) for x in set(rel_avgs_dico).union(cur_rel)}
    for x in set(abs_avgs_dico).union(cur_abs):
        if not x in abs_avgs_dico:
            abs_avgs_dico[x] = 0
            rel_avgs_dico[x] = 0
        if not x in cur_abs:
            cur_abs[x] = 0
            cur_rel[x] = 0
        abs_avgs_dico[x] = (float(abs_avgs_dico[x]) * rep_counter + float(cur_abs[x])) / (rep_counter + 1)
        rel_avgs_dico[x] = (float(rel_avgs_dico[x]) * rep_counter + float(cur_rel[x])) / (rep_counter + 1)
    rep_counter += 1
metrics.save_dico(root_output, abs_avgs_dico, tag + "_global__rf_distance_avg-abs")
metrics.save_dico(root_output, rel_avgs_dico, tag + "_global__rf_distance_avg-rel")
print("seeds = ", seeds)