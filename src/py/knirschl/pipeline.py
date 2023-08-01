import os
import sys
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'scripts/programs')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/simulation')
sys.path.insert(0, 'tools/trees')
import evaluate
import compare_trees
import dist_matrix_converter
import fam
import generate_with_simphy as simphy
import launch_ba
import launch_fastme
import launch_generax
import launch_raxml
import metrics
import paths
import utils


class RunFilter():
    def __init__(self, generate = True):
        self.disable_all()
        self.generate = True
        self.raxml = True
        self.generax = True
        self.fastme = True
        self.ba = True
        self.ba_fastme = True
        self.generax_pick = True
        self.compare = True
    
    def disable_all(self):
        self.generate = False
        self.force_overwrite = False
        self.raxml = False
        self.generax = False
        self.fastme = False
        self.ba = False
        self.ba_fastme = False
        self.generax_pick = False
        self.compare = False
    
    def sim(self):
        self.generate = True
        self.force_overwrite = False
        self.raxml = False
        self.generax = False
        self.fastme = False
        self.ba = False
        self.compare = False

    def bacomp_full(self):
        self.generate = False
        self.raxml = False
        self.generax = False
        self.fastme = False
        self.ba = True
        self.ba_fastme = True
        self.generax_pick = True
        self.compare = True
    
    def pick_comp(self):
        self.generate = False
        self.raxml = False
        self.generax = False
        self.fastme = False
        self.ba = False
        self.ba_fastme = False
        self.generax_pick = True
        self.compare = True

    def comp(self):
        self.generate = False
        self.raxml = False
        self.generax = False
        self.fastme = False
        self.ba = False
        self.ba_fastme = False
        self.generax_pick = False
        self.compare = True

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
            utils.printFlush("Run raxml-ng...\n***************")
            try:
                launch_raxml.run_raxmlng_on_families(datadir, subst_model, cores)
            except Exception as exc:
                utils.printFlush("Failed running RAxML-NG\n" + str(exc))
        if (self.generax):
            utils.printFlush("Run generax...\n**************")
            try:
                species_tree = fam.get_species_tree(datadir)
                resultsdir = fam.get_run_dir(datadir, subst_model, "generax_run")
                launch_generax.run(datadir, subst_model, "SPR", species_tree, "random", cores, "--rec-model UndatedDL", resultsdir)
            except Exception as exc:
                utils.printFlush("Failed running GeneRax\n" + str(exc))
        if (self.fastme):
            utils.printFlush("Run fastme...\n*************")
            try:
                launch_fastme.run_fastme_on_families(datadir, subst_model, is_dna=True, algo="B", use_spr=True, only_mat=False, cores=cores)
            except Exception as exc:
                utils.printFlush("Failed running FastME\n" + str(exc))
        if (self.ba):
            utils.printFlush("Run ba...\n*********")
            try:
                start = time.time()
                # convert species tree and alignments to distance matrix
                dist_matrix_converter.convert_input(datadir, cores)
                species_tree = fam.get_true_species_tree_matrix(datadir)
                # run ba script
                inferred_trees = launch_ba.run_ba_on_families(datadir, "p", species_tree,
                                                              mat_out=int(self.ba_fastme),
                                                              cores=cores)
                print("=#=#= BA-Code took {}s per tree =#=#=".format((time.time() - start) / ((int)(simphy.get_param_from_dataset_name("families", datadir)) * inferred_trees)))
            except Exception as exc:
                utils.printFlush("Failed running bachelor thesis script\n" + str(exc))
        if (self.ba_fastme):
            utils.printFlush("Run ba matrix with fastme trees...\n**********************************")
            try:
                if (not self.ba):
                    # convert species tree and alignments to distance matrix
                    # dist_matrix_converter.convert_input(datadir, cores)
                    species_tree = fam.get_true_species_tree_matrix(datadir)
                    # run ba script
                    inferred_trees = launch_ba.run_ba_on_families(datadir, "p", species_tree,
                                                                  algo="APro", mat_out=2, cores=cores)
                launch_fastme.run_fastme_on_families_matrices(datadir, "ba.p", algo="B",
                                                              use_spr=True, cores=cores)
            except Exception as exc:
                utils.printFlush("Failed running bachelor thesis matrix correction with FastME\n" + str(exc))
        if (self.generax_pick):
            utils.printFlush("Run picking...\n**************")
            try:
                # run generax evaluation and select best tree
                species_tree = fam.get_species_tree(datadir)
                resultsdir = fam.get_run_dir(datadir, subst_model, "generax_eval_run")
                launch_generax.run(datadir, subst_model, "EVAL", species_tree, "ba", cores, ["--rec-model", "UndatedDL", "--per-family-rates"], resultsdir, False)
            except Exception as exc:
                utils.printFlush("Failed running pick\n" + str(exc))
        # COMPARE INFERRED TREES WITH TRUE TREE
        if (self.compare):
            utils.printFlush("Run compare...\n**************")
            try:
                compare_trees.compare_all(datadir)
            except Exception as exc:
                utils.printFlush("Failed running compare\n" + str(exc))

def pipeline(datadir, run_filter, seed, tag):
    # RUN PIPELINE
    rep_start = time.time()
    try:
        run_filter.run_methods(datadir, "F81", 8)
    finally:
        elapsed = time.time() - rep_start
        print("End of single experiment. Elapsed time: " + str(elapsed) + "s")
        metrics.save_metrics(datadir, "pipeline_" + tag + str(seed), elapsed, "runtimes")

def run_pipeline(enabled = True):
    # TOGGLE PIPELINE ELEMENTS
    # ====== ! CAREFUL ! ======
    run_filter = RunFilter()  # all enabled
    # run_filter.force_overwrite = True # regenerate old dataset
    # run_filter.raxml = False
    # run_filter.generax = False
    # run_filter.force_overwrite = False
    run_filter.bacomp_full()
    run_filter.ba = False  # keep this False
    # run_filter.pick_comp() # only compare inferred trees
    #run_filter.disable_all() # collect avgs
    # run_filter.generax_pick = True
    # ====== ! CAREFUL ! ======

    root_output = paths.families_datasets_root  # output/families/
    # seeds 
    seeds5 = [42, 1007, 19732311, 121873, 14976684177860080345]
    seeds50 = [45364, 254864, 422868, 710192, 1142238, 1158688, 1424996, 1472190, 1513197, 1650734, 1716696, 1994998,
               2190362, 2289044, 2436736, 2615935, 2620431, 2661376, 2722114, 2908180, 3044371, 3245068, 4040147,
               4156463, 4257503, 4683309, 5087218, 5198455, 5355281, 5536663, 5543415, 5755677, 5942622, 6089587,
               6329094, 6783683, 6947993, 7132734, 7308089, 7406226, 7781467, 8643158, 8726470, 8771421, 8846466,
               9133069, 9376940, 9387366, 9481423, 9724682]
    seeds100 = [42, 1007, 4229, 11872, 121873, 125546, 340025, 429510, 436633, 491587, 558959, 842036, 901464, 1058026,
                1091512, 1105070, 1143740, 1316419, 1320646, 1570525, 1674005, 2181913, 2366262, 2453027, 2525498,
                2530855, 2545197, 2622038, 2650353, 2835503, 2862791, 2967039, 2967641, 3174338, 3219882, 3279056,
                3389159, 3547108, 3614207, 3686221, 3939144, 3959336, 4108907, 4130562, 4144379, 4193573, 4290505,
                4354376, 4412646, 4602564, 4754615, 4771372, 4785291, 4824437, 5476033, 5597160, 5616680, 5749959,
                6050422, 6077656, 6214982, 6382595, 6449414, 6645706, 6722758, 6734655, 7060210, 7265233, 7340707,
                7342556, 7350931, 7538361, 7724757, 7927143, 7940684, 8039292, 8124350, 8386866, 8388955, 8438226,
                8527944, 8563774, 8670475, 8710215, 8928866, 9070912, 9128753, 9179296, 9387015, 9442585, 9476398,
                9488374, 9495297, 9510809, 9718279, 9919929, 9930269, 9955696, 19732311, 14976684177860080345]
    seeds = seeds50
    # while (len(seeds) != 50):
    #    seeds.append(random.randrange(0, 9999999))
    # simphy params
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
        if enabled:
            pipeline(datadir, run_filter, seed, tag)
    return root_output, seeds, tag, replicates

if (__name__ == "__main__"):
    start = time.time()
    root_output, seeds, tag, reps = run_pipeline(True)
    best_avg_tree, _ = evaluate.global_compare(root_output, reps, tag)
    evaluate.collect_generax_picks(root_output, reps, tag, True)
    evaluate.generax_likelihood_comp(root_output, reps, best_avg_tree, os.path.join("runs", "F81", "generax_eval_run"))
    print("seeds = ", seeds)
    print("End of pipeline. Elapsed time:", time.time() - start)
