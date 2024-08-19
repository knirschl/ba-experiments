import os
import sys
import subprocess
import time
import random
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'scripts/programs')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/msa')
sys.path.insert(0, 'tools/simulation')
sys.path.insert(0, 'tools/trees')
import evaluate
import compare_trees
import fam
import generate_with_simphy as simphy
import launch_fastme
import launch_generax
import launch_raxml
import metrics
import paths
import utils

class RunFilter():
    def __init__(self):
        self.all(True)
        self.force_overwrite = False 
        self.is_dna = True
        self.subst_model = "F81"

    def all(self, t):
        self.generate = t
        self.raxml = t
        self.generax = t
        self.fastme = t
        self.ba = t
        self.ba_fastme = t
        self.spearfish = t
        self.generax_pick = t
        self.compare = t
    
    def disable_all(self):
        self.all(False)

    def run_methods(self, datadir, cores):
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
        #redirected_file = os.path.join(datadir, "runs", "logs_run_all_genes." + self.subst_model + ".txt")
        #print("Redirected logs to " + redirected_file)
        sys.stdout.flush()
        # RUN
        if(self.raxml):
            utils.printFlush("Run raxml-ng...\n***************")
            try:
                launch_raxml.run_raxmlng_on_families(datadir, self.subst_model, cores)
            except Exception as exc:
                utils.printFlush("Failed running RAxML-NG\n" + str(exc))
        if (self.generax):
            utils.printFlush("Run generax...\n**************")
            try:
                species_tree = fam.get_species_tree(datadir)
                resultsdir = fam.get_run_dir(datadir, self.subst_model, "generax_run")
                launch_generax.run(datadir, self.subst_model, "SPR", species_tree, "random", cores, ["--rec-model UndatedDL"], resultsdir)
            except Exception as exc:
                utils.printFlush("Failed running GeneRax\n" + str(exc))
        if (self.fastme):
            utils.printFlush("Run fastme...\n*************")
            try:
                launch_fastme.run_fastme_on_families(datadir, self.subst_model, is_dna=self.is_dna, algo="B", use_spr=True, only_mat=False, cores=cores)
            except Exception as exc:
                utils.printFlush("Failed running FastME\n" + str(exc))
        if (self.spearfish):
            utils.printFlush("Run Spearfish benchmarking...\n*****************************")
            try:
                # wrapper dataset dna subst_model cores compute algo
                subst_model = "p/" + self.subst_model # TODO for testing, later: self.subst_model
                command = [paths.python3(), paths.spearfish_wrapper_py, datadir, subst_model, str(int(self.is_dna)), str(cores), "test"]
                subprocess.check_call(command, stdout=sys.stdout, cwd=os.path.join(paths.root, "Spearfish", "src", "py"))
            except Exception as exc:
                utils.printFlush("Failed running Spearfish benchmarking\n" + str(exc))
        # COMPARE INFERRED TREES WITH TRUE TREE
        if (self.compare):
            utils.printFlush("Run compare...\n**************")
            try:
                compare_trees.compare_all(datadir)
            except Exception as exc:
                utils.printFlush("Failed running compare\n" + str(exc))


def get_run_filter(run_filter_str):
    run_filter_str = run_filter_str.replace("full", "rgfbpc")
    # TOGGLE PIPELINE ELEMENTS
    # ====== ! CAREFUL ! ======
    run_filter = RunFilter()  # all enabled
    # run_filter.force_overwrite = True # regenerate old dataset
    run_filter.disable_all()
    if ("s" in run_filter_str):
        run_filter.generate = True
    if ("r" in run_filter_str):
        run_filter.raxml = True
    if ("g" in run_filter_str):
        run_filter.generax = True
    if ("f" in run_filter_str):
        run_filter.fastme = True
    if ("b" in run_filter_str):
    #    run_filter.ba = False # keep this False as it's worse than ba+fm
    #    run_filter.ba_fastme = True
        run_filter.spearfish = True
    #if ("p" in run_filter_str):
    #    run_filter.generax_pick = True
    if ("c" in run_filter_str):
        run_filter.compare = True
    subst_model_start = run_filter_str.find("_")
    if (subst_model_start != -1):
        run_filter.is_dna = run_filter_str[subst_model_start + 1] == "D"
        run_filter.subst_model = run_filter_str[subst_model_start + 2:]
    return run_filter


def get_seeds(rank, n_families):
    seeds100 = [42, 1007, 39104, 45364, 121873, 178811, 254864, 364465, 422868,
                592240, 710192, 733230, 785319, 1142238, 1158688, 1230166,
                1381421, 1424996, 1472190, 1513197, 1650734, 1656898, 1690222,
                1716696, 1886712, 1990093, 1994998, 2073911, 2133265, 2190362,
                2289044, 2436736, 2615935, 2620431, 2661376, 2722114, 2895719,
                2908180, 3044371, 3072898, 3245068, 3633308, 3799911, 4040147,
                4156463, 4257503, 4413346, 4683309, 4788078, 5012216, 5019389,
                5087218, 5182115, 5198455, 5287118, 5334140, 5355281, 5398355,
                5536663, 5543415, 5574927, 5755677, 5942622, 5966994, 6063401,
                6089587, 6317675, 6329094, 6396236, 6503199, 6783683, 6947993,
                7132734, 7308089, 7406226, 7425579, 7555193, 7781467, 7966578,
                8142353, 8197298, 8199273, 8534020, 8643158, 8726470, 8771421,
                8821971, 8846466, 8850161, 9037043, 9133069, 9300400, 9316715,
                9376940, 9387366, 9438631, 9481423, 9724682, 9824219, 19732311]
    start_rep = rank * n_families
    seeds = seeds100[start_rep:start_rep + n_families]
    while (len(seeds) < n_families):
        seeds.append(random.randrange(0, 9999999))
    return seeds


def run(datadir, run_filter, tag, seed=-1):
    if (run_filter.spearfish):
        # Let GeneRax eval not skip if it already run
        try:
            os.remove(os.path.join(datadir, "runs", "F81", "generax_eval_run", "gene_optimization_0", "checkpoint_commands.txt"))
            os.remove(os.path.join(datadir, "runs", "F81", "generax_eval_run", "gene_optimization_0", "statistics.svg"))
            print("Removed generax checkpoints")
        except:
            print("Failed to rm in", datadir)
    # RUN PIPELINE
    rep_start = time.time()
    try:
        run_filter.run_methods(datadir, cores=16)
    finally:
        elapsed = time.time() - rep_start
        print("End of single experiment. Elapsed time: " + str(elapsed) + "s")
        metrics.save_metrics(datadir, "pipeline_" + tag + str(seed) * (seed > -1), elapsed, "runtimes")


def pipeline_sim_data(argv):
    # ((exec "sim")) tag tag_val rank n_families run_filter
    if (len(argv) != 5):
        print("Syntax: python scripts/pipeline.py \"sim\" tag tag_val rank n_families run_filter")
        sys.exit(1)
    tag = argv[0] # e.g. SPECIES
    tag_val = argv[1] # e.g. 75
    rank = int(argv[2]) # rank-th pipeline running on this tag-tag_val combination
    n_families = int(argv[3]) # families resolved in this pipeline

    start = time.time()
    # Prepare for run
    run_filter = get_run_filter(argv[4])
    seeds = get_seeds(rank, n_families)
    root_output = paths.families_datasets_root
    # SimPhy parameters
    # BASE dataset
    s = 25
    f = 100
    sites = 100
    bl = 1.0
    d = 1.0
    l = 1.0
    t = 0.0
    pop = 10
    # Specialized datasets
    match (tag):
        case "SPECIES":
            s = int(val)
        case "FAM":
            f = int(val)
        case "SITES":
            sites = int(val)
        case "BRALEN":
            bl = float(val)
        case "DUPLOS":
            d = l = float(val)
        case "DL":
            d = l = float(val)
        case "TRA":
            t = float(val)
        case "T":
            t = float(val)
        case "DUPLOSTRA":
            val = val.split(',')
            d = l = float(val[0])
            t = float(val[1])
        case "DTL":
            val = val.split(',')
            d = l = float(val[0])
            t = float(val[1])
        case"POP": # ils
            pop = int(val)
    # Run multiple replicates
    for seed in seeds:
        # Set SimPhy parameters
        simphy_parameters = simphy.SimphyParameters(tag=tag, species_taxa=s, families_number=f, 
            sites=sites, bl=bl, dup_rate=d, loss_rate=l, transfer_rate=t, population=pop, seed=seed)
        datadir = simphy.get_output_dir(simphy_parameters, root_output)
        run(datadir, run_filter, tag, seed)
    print("seeds =", seeds)
    print(f"End of pipeline ({tag + tag_val}). Elapsed time: {time.time() - start}")


def pipeline_real_data(argv):
    # ((exec "real")) tag rank n_families run_filter
    if (len(argv) != 4):
        print("Syntax: python scripts/pipeline.py \"real\" tag rank n_families run_filter")
        sys.exit(1)
    tag = argv[0] # full directory qualifier
    rank = int(argv[1]) # unused, could be used for family splitting
    n_families = int(argv[2]) # unused, could be used for family splitting

    start = time.time()
    # Prepare for run
    run_filter = get_run_filter(argv[3])
    datadir = fam.get_datadir(tag)
    run(datadir, run_filter, tag)
    print("dataset =", tag)
    print(f"End of pipeline ({tag}). Elapsed time: {time.time() - start}")


def tagval_in_string(string, tag, tag_val):
    match(tag):
        case "SPECIES":
            return "_s" + tag_val in string
        case "FAM":
            return "_f" + tag_val in string
        case "SITES":
            return "_sites" + tag_val in string
        case "BRALEN":
            return "_bl" + tag_val in string
        case "DUP":
            return "_d" + tag_val in string
        case "LOS":
            return "_l" + tag_val in string
        case "LOSS":
            return "_l" + tag_val in string
        case "DUPLOS":
            return "_d" + tag_val in string and "_l" + tag_val in string
        case "DL":
            return "_d" + tag_val in string and "_l" + tag_val in string
        case "TRA":
            return "_t" + tag_val in string
        case "T":
            return "_t" + tag_val in string
        case "DUPLOSTRA":
            tag_val = tag_val.split(",")
            return "_d" + tag_val[0] in string and "_l" + tag_val[0] in string and "_t" + tag_val[1] in string
        case "POP":
            return "pop" + tag_val in string


def pipeline_eval_data(argv):
    # ((exec "eval")) tag
    if (len(argv) < 1 or len(argv) > 2):
        print("Syntax: python scripts/pipeline.py \"eval\" tag [tag_val]")
        sys.exit(1)
    tag = argv[0]
    root_output = paths.families_datasets_root

    start = time.time()
    # Collect datasets
    replicates = [] # datasets with same tag & tag_val (probably differ only in seed)
    #replicates.append(os.path.join(root_output, "generax_data/cyano_empirical"))
    for dataset in os.listdir(root_output):
        if (not tag in dataset):
            continue
        if (len(argv) == 2):
            if (not (tag == "BASE" or tagval_in_string(dataset, tag, argv[1]))):
                continue
        replicates.append(os.path.join(root_output, dataset))
    if (len(argv) == 2):
        tag += argv[1]
    # Evaluate datasets
    best_avg_tree, best_avg_dist = evaluate.global_compare(root_output, replicates, tag)
    print("Best global results:", best_avg_tree, best_avg_dist)
    evaluate.collect_generax_picks(root_output, replicates, tag, compare=True)
    #evaluate.generax_likelihood_comp(root_output, replicates, tag, best_avg_tree, os.path.join("runs", "F81", "generax_eval_run"))
    print(f"End of evaluation of {tag} datasets. Elapsed time: {time.time() - start}") 

if (__name__ == "__main__"):
    match (sys.argv[1]):
        case "sim":
            pipeline_sim_data(sys.argv[2:])
        case "real":
            pipeline_real_data(sys.argv[2:])
        case "eval":
            pipeline_eval_data(sys.argv[2:])
        case _:
            print("Syntax: python scripts/pipeline.py \"<sim/real/eval>\" [args...]")
            sys.exit(1)