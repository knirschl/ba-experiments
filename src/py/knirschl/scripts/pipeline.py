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
        self.disable_all()
        self.generate = True
        self.force_overwrite = False

    def bacomp_full(self):
        self.disable_all()
        self.ba = True
        self.ba_fastme = True
        self.generax_pick = True
        self.compare = True
    
    def pick_comp(self):
        self.disable_all()
        self.generax_pick = True
        self.compare = True

    def comp(self):
        self.disable_all()
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
                launch_generax.run(datadir, subst_model, "SPR", species_tree, "random", cores, ["--rec-model UndatedDL"], resultsdir)
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
                    dist_matrix_converter.convert_input(datadir, cores)
                    species_tree = fam.get_true_species_tree_matrix(datadir)
                    # run ba script
                    inferred_trees = launch_ba.run_ba_on_families(datadir, "p", species_tree, algo="APro", mat_out=2, cores=cores)
                launch_fastme.run_fastme_on_families_matrices(datadir, "ba.p", algo="B",        use_spr=True, cores=cores)
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

def run_pipeline(reps = 50, tag = "DL", val=0, run_filter_str = "ba-cp", enabled = True):
    # TOGGLE PIPELINE ELEMENTS
    # ====== ! CAREFUL ! ======
    run_filter = RunFilter()  # all enabled
    # run_filter.force_overwrite = True # regenerate old dataset
    if (run_filter_str == "sim"):
        run_filter.sim()
    elif (run_filter_str == "full"):
        run_filter.generate = False
    elif (run_filter_str == "fm-ba-cp"):
        run_filter.generate = False
        run_filter.generax = False
        run_filter.raxml = False
    elif (run_filter_str == "ba-cp"):
        run_filter.bacomp_full()
    elif (run_filter_str == "pc"):
        run_filter.pick_comp()
    elif (run_filter_str == "comp"):
        run_filter.comp()
    # DEBUG
    run_filter.ba = False  # keep this False
    # run_filter.disable_all() # run selected
    # run_filter.generax_pick = True

    # SEEDS
    seeds100 = [42, 1007, 39104, 45364, 121873, 178811, 254864, 364465, 422868, 592240, 710192, 733230, 785319, 1142238, 1158688, 1230166, 1381421, 1424996, 1472190, 1513197, 1650734, 1656898, 1690222, 1716696, 1886712, 1990093, 1994998, 2073911, 2133265, 2190362, 2289044, 2436736, 2615935, 2620431, 2661376, 2722114, 2895719, 2908180, 3044371, 3072898, 3245068, 3633308, 3799911, 4040147, 4156463, 4257503, 4413346, 4683309, 4788078, 5012216, 5019389, 5087218, 5182115, 5198455, 5287118, 5334140, 5355281, 5398355, 5536663, 5543415, 5574927, 5755677, 5942622, 5966994, 6063401, 6089587, 6317675, 6329094, 6396236, 6503199, 6783683, 6947993, 7132734, 7308089, 7406226, 7425579, 7555193, 7781467, 7966578, 8142353, 8197298, 8199273, 8534020, 8643158, 8726470, 8771421, 8821971, 8846466, 8850161, 9037043, 9133069, 9300400, 9316715, 9376940, 9387366, 9438631, 9481423, 9724682, 9824219, 19732311]
    seeds = seeds100[:reps]
    while (len(seeds) < reps):
       seeds.append(random.randrange(0, 9999999))
    
    # DATASET PARAMETERS
    # base
    s = 25
    f = 100
    sites = 100
    bl = 1.0
    d = 1.0
    l = 1.0
    t = 0.0
    pop = 10
    if (tag == "SPECIES"):
        s = int(val)
    elif (tag == "FAM"):
        f = int(val)
    elif (tag == "SITES"):
        sites = int(val)
    elif (tag == "BRALEN"):
        bl = float(val)
    elif (tag == "DUPLOS" or tag == "DL"):
        d = l = float(val)
    elif (tag == "TRA" or tag == "T"):
        t = float(val)
    elif (tag == "DUPLOSTRA" or tag == "DTL"):
        val = val.split(',')
        d = l = float(val[0])
        t = float(val[1])
    elif (tag == "POP"): # ils
        pop = int(val)
    
    root_output = paths.families_datasets_root  # output/families/
    replicates = []
    # Run multiple replicates
    for seed in seeds:
        # SET simphy PARAMETERS 
        simphy_parameters = simphy.SimphyParameters(tag=tag, species_taxa=s, families_number=f, 
            sites=sites, bl=bl, dup_rate=d, loss_rate=l, transfer_rate=t, population=pop, seed=seed)
        datadir = simphy.get_output_dir(simphy_parameters, root_output)
        replicates.append(datadir)
        if enabled:
            pipeline(datadir, run_filter, seed, tag)
    return root_output, seeds, tag, replicates

if (__name__ == "__main__"):
    if (len(sys.argv) != 7):
        print("Syntax: python scripts/pipeline.py reps tag tag_val run_filter enable_pip compare_picks")
        sys.exit(1)
    reps = int(sys.argv[1])
    tag = sys.argv[2]
    tag_val = sys.argv[3]
    run_filter_str = sys.argv[4]
    enable_pip = int(sys.argv[5])
    compare_picks = int(sys.argv[6])

    start = time.time()
    root_output, seeds, tag, reps = run_pipeline(reps, tag, tag_val, run_filter_str, enable_pip)
    if (run_filter_str != "sim"):
        best_avg_tree, _ = evaluate.global_compare(root_output, reps, tag)
        evaluate.collect_generax_picks(root_output, reps, tag, compare_picks)
        evaluate.generax_likelihood_comp(root_output, reps, best_avg_tree, os.path.join("runs", "F81", "generax_eval_run"))
    print("seeds = ", seeds)
    print("End of pipeline. Elapsed time:", time.time() - start)
