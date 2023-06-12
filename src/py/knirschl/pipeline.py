import sys
import os
import time
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

        self.force_overwrite = True
        self.compare = True
        #self.debug = False
    
    def disable_all(self):
        self.simphy = False
        self.raxml = False
        self.generax = False
        self.fastme = False
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
            # TODO
            try:
                dist_matrix_converter.convert_input(datadir)
                species_tree = fam.get_true_species_tree_matrix_sorted(datadir)
                launch_ba.run_ba_on_families(datadir, "experimental", species_tree, cores)
            except Exception as exc:
                utils.printFlush("Failed running bachelor thesis script\n" + str(exc))
        # COMPARE INFERRED TREES WITH TRUE TREE
        if (self.compare):
            utils.printFlush("Run compare...")
            try:
                compare_trees.compare_all(datadir)
            except Exception as exc:
                utils.printFlush("Failed running compare\n" + str(exc))

root_output = paths.families_datasets_root # output/families/
simphy_parameters = simphy.SimphyParameters()
datadir = simphy.get_output_dir(simphy_parameters, root_output)
print(datadir)
run_filter = RunFilter()
run_filter.generate = False
run_filter.force_overwrite = False
#run_filter.raxml = False
#run_filter.generax = False
#run_filter.fastme = False
run_filter.ba = False
run_filter.compare = False
start = time.time()
try:
    run_filter.run_methods(datadir, "F81", 1)
finally:
    elapsed = time.time() - start
    print("End of single experiment. Elapsed time: " + str(elapsed) + "s")
    metrics.save_metrics(datadir, "pipeline", elapsed, "runtimes")
