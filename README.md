# Species tree aware gene tree inference with distance methods

This repository was created during my bachelore thesis. It contains all code
for conducting my experiments, generating graphs from their results, as
well as the implementation of Spearfish.

## TODO

- [ ] README
  - [X] Spearfish full name
  - [X] TODOs
  - [X] Usage
    - [X] Dependencies
- [ ] Spearfish
  - [ ] Extract Spearfish into own repository
  - [ ] Clean up
  - [ ] Optimize
  - [ ] Own README
    - [ ] Python dependencies
  - License
- [ ] src/
  - [ ] Clean up
    - [ ] Use less 3rd party 

### Usage

To run the experiments, either edit the settings of `src/py/knirschl/launcher.py`
or call `src/py/knirschl/scripts(pipeline.py)` directly for a single experiment.
The launcher was designed to work on clusters managed by slurm and creates a
slurm submit file accordingly.  
The `run_filter` can be used to run only a subset of tools:
  - `s`: Generate new synthetic datasets with SimPhy and INDELIble
  - `r`: Run RAxML-NG
  - `g`: Run GeneRax
  - `f`: Run FastME
  - `b`: Run code developed for the thesis (i.e. Spearfish)
  - `p`: Run GeneRax with (i.e. Spearfish)
  - `c`: Compare all existing trees against true tree
  - `full`: Run all of the above but don't generate new datasets (i.e. `rgfbpc`)
The `bm_set` value selects which datasets the tools specified in `run_filter`
should run on. The specific benchmarks (`benchmarks` and `benchmarks_ext`) are
defined in `launcher.py` and can be extended. If a new key is added, it has to set the SimPhy paramters accordingly in `pipeline.py#run_pipeline()`.

### Dependencies

All external programs (see "Executables" and "Tools that are compared") need to
be build first. Their paths are defined in `src/py/knirschl/scripts/paths.py` and should be updated accordingly.

#### Python

  - ETE3
    - Tree
    - SeqGroup
  - Bio
    - AlignIO
    - Phylo.TreeConstruction
  - PhyloDM
  - NumPy

For plotting the results:
  - Matplotlib

#### Executables

  - [SimPhy](https://github.com/adamallo/SimPhy)
  - [INDELIble](http://abacus.gene.ucl.ac.uk/software/indelible)
  - [MPIScheduler](https://github.com/BenoitMorel/MPIScheduler)

#### Tools that are compared

  - [RAxML-NG](https://github.com/amkozlov/raxml-ng)
  - [GeneRax](https://github.com/BenoitMorel/GeneRax)
  - [FastME](http://www.atgc-montpellier.fr/fastme)


# Spearfish

<p align="center">
  <img src="resources/spearfish.png" alt="Spearfish logo"/>
</p>
Spearfish (<b>SPE</b>cies tree <b>A</b>ware gene t<b>R</b>ee in<b>F</b>erence
with d<b>IS</b>tance met<b>H</b>ods) is a distance-based tree inference
algorithm that respects the species tree and can infer gene family trees. 

### Usage

### Dependencies

Spearfish uses multiple Python packages to pre-process the input data, as well
as some phylogenetic tools that need to be build in order to run Spearfish with
optimal settings.

#### Phylogenetic tools

  - [FastME](http://www.atgc-montpellier.fr/fastme)
  - [MADroot](https://github.com/davidjamesbryant/MADroot)
  - [GeneRax](https://github.com/BenoitMorel/GeneRax)

#### Others

  - [argparse](https://github.com/p-ranav/argparse)