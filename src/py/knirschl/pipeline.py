import sys
import os
import time
sys.path.insert(0, 'scripts')
sys.path.insert(0, 'scripts/programs')
sys.path.insert(0, 'tools/families')
sys.path.insert(0, 'tools/simulation')
import utils
import paths
import fam
import generate_with_simphy as simphy
import launch_raxmlng__MY as launch_raxmlng
import launch_generax
import launch_fastme

class RunFilter():
    def __init__(self, generate = True):
        self.disable_all()
        self.generate = True
        self.raxml = True
        self.generax = True
        self.fastme = True
        self.ba = True

        self.analyze = False # TODO implement
        #self.debug = False
    
    def disable_all(self):
        self.simphy = False
        self.raxml = False
        self.generax = False
        self.fastme = False
        self.ba = False

    def run_methods(self, datadir, subst_model, cores):
        if (self.generate):
            if (not os.path.isdir(datadir)):
                dataset = os.path.basename(datadir)
                simphy.generate_dataset(dataset)
        print("*************************************")
        print("Run tested gene tree inference tools for dataset " + datadir)
        print("*************************************")
        if (len(datadir.split("/")) == 1):
          datadir = fam.get_datadir(datadir) 
        save_stdout = sys.stdout
        redirected_file = os.path.join(datadir, "runs", "logs_run_all_genes." + subst_model + ".txt")
        print("Redirected logs to " + redirected_file)
        sys.stdout.flush()
        # RUN
        if(self.raxml):
            utils.printFlush("Run raxml-ng...")
            try:
                # TODO
                # -- my version --
                launch_raxmlng.run_raxml_all(datadir, subst_model)
            except Exception as exc:
                utils.printFlush("Failed running RAxML-NG\n" + str(exc))
        if (self.generax):
            utils.printFlush("Run generax...")
            try:
                # TODO
                #run_generax.run_generax_on_families(datadir, subst_model, cores)
                # -- my version --
                # create families file
                generax_families_file = os.path.join(datadir, "families.txt")
                launch_generax.build_generax_families_file(datadir, "", subst_model, generax_families_file)
                resultsdir = os.path.join(datadir, "misc")
                launch_generax.run_generax(datadir, None, "SPR", "true", generax_families_file, "false", 1, resultsdir)
            except Exception as exc:
                utils.printFlush("Failed running GeneRax\n" + str(exc))
        if (self.fastme):
            utils.printFlush("Run fastme...")
            try:
                # TODO
                # -- my version --
                output = "misc/fastme-test.newick"
                launch_fastme.run_fastme_all(datadir, output)
                # maybe use benoitmorel/tools/families/run_fastme.py ?
                run_fastme.run_fastme_on_families(datadir, subst_model, 1, 1)
            except Exception as exc:
                utils.printFlush("Failed running FastME\n" + str(exc))
        if (self.ba):
            utils.printFlush("Run ba...")
            # TODO
        # ANALYZE
        if (self.analyze):
            utils.printFlush("Run analyze...")
            try:
                # TODO
                tmp = 0
                #fast_rf_cells.analyze(datadir, "all", cores)
            except Exception as exc:
                utils.printFlush("Failed running analyze\n" + str(exc))


root_output = paths.families_datasets_root # output/families/
simphy_parameters = simphy.SimphyParameters()
datadir = simphy.get_output_dir(simphy_parameters, root_output)
run_filter = RunFilter()
start = time.time()
try:
    run_filter.run_methods(datadir, "GTR+G", 1)
finally:
    elapsed = time.time() - start
    print("End of single experiment. Elapsed time: " + str(elapsed) + "s")
