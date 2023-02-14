# README

## Introduction

## Installation

To compile epics, epocs and epocs_mcmc, no modules are needed. Simply type `make`, which will build the three executables.

To be able to use evo-scope, you require `pastml` > v1.9.3, `epics`, `epocs` and `R` > v4.0.0 with `tidyverse` and `ape` installed as packages. These tools can be installed with bioconda.

## Usage

The details on how to use each separate tools, as well as the whole pipeline encapsulated in `evo-scope`, are presented hereafter.

### evo-scope

`evo-scope` works as a 5-step pipeline, combining the specificities of two of our previous tools, `epics` and `epocs`.

The only inputs required are a rooted phylogenetic tree and a tab-delimited file with discrete traits at the tips. The tip IDs need to correspond to the tips in the tree.

The first step is an ancestral character reconstruction performed by `pastml` with default parameters and the `JOINT` method (see the [github repository](https://github.com/evolbioinfo/pastml) and the [reference paper](https://academic.oup.com/mbe/article/36/9/2069/5498561?login=false)).

The second step in `evo-scope` reconstructs the mutations on the branches of the tree with a tailored Rscript `AnnotateATree_WithMasks_Clean.R`. Indeed, `epics` and `epocs` require as input a rooted phylogenetic tree with mutations formatted between square brackets and located at the bootstrap values in the newick string for internal branches (e.g., `)[0/1/0/1/1]`) or after the tip ids (e.g., `SPECIES1[0/0/0/1/1/1]`).

Next, `epics` is run with the identity matrix (see the [epics paper](https://academic.oup.com/sysbio/article/65/5/812/2223542)) to identify pairs of significantly co-occurring mutations within the tree. In `evo-scope`, epics works as a pre-filter to collect quickly evidence of correlated evolution, as the calculation speed allows it. We keep pairs that are significantly associated with a p-value < 0.05.

Following this step, we feed in `epocs` the significant pairs from the previous step to have an idea about the induction scenarios that describe the associations. We run four models within `epocs`, one of independence and three of dependence (see also the [epocs paper](https://doi.org/10.1093/sysbio/syab092)). The dependence models selected reflect either an induction of the first trait on the second one, vice-versa or reciprocal induction.

Finally, `evo-scope` uses a tailored R script `ParseEpocsFinal_WithLambda.R` to retrieve the results of all four runs of the previous steps. `evo-scope` then calculates the likelihood ratio tests between the encapsulated models of dependence and independence. The final output consists of three files: one `[prefix]_BestModel.tab` where for each pair whose LRT is significant, the best model is presented based on the LRT values; `[prefix]_LRT.tab` summarizing all LRT tests performed; `[prefix]_ConcatenatedResults.tab` containing the results of all runs of `epocs`.

The entire pipeline is fully automated and do not require any manual inputs.

Some options are available to tweak `evo-scope`, available after typing `./evoscope -h`:

```
Usage: ./evoscope [options]
Options:
    -t [FILENAME]: tab-delimited file of traits at the leaves [mandatory]
    -T [TREEFILE]: Phylogenetic tree in newick format [mandatory]
    -O [OUTGROUP]: Outgroup [Optional, default = no outgroup. format= "name1[:name2:name3:...]"]
    -o [PREFIX]  : Output prefix [Optional, default = "out"]
    -m [MATRIX]  : epics matrix [Optional, default=Identity]
    -g           : GWAS-like analysis [Optional, default runs all against all]
    -f           : full run [Optional, default runs until epics]
    -p           : Plotting the results [Optional, default does not output plot files]
    -s           : State-dependence analysis in epocs step [Optional, default not]
    -h           : Display this help message
```

### Separate tools

Any of the tools and scripts are usable outside of evo-scope. Arguments and usage are detailed hereafter.

#### Epics

Epics calculates very quickly a p-value describing the statistical association between pair of mutations placed on a rooted phylogenetic tree. Depending on the matrix choice, the user can either test whether pairs are statistically co-occurring or are genealogically ordered in the tree. We refer the user to the [original paper](https://academic.oup.com/sysbio/article/65/5/812/2223542) for details on the implementation.

type `./epics -h` to see the help and the main commands.

The standard format to run epics is as follows

`./epics [OPTIONS] "OUTGROUP" PHYLO_TREE`, where outgroup is either empty when no outgroup is present (`epics [OPTIONS] "" PHYLO_TREE`) or takes the form `"name1[:name2:name3:...]"`.

### List of options for Epics

```
[matrix choice]
    -I: inseparable pairs only will be computed (Identity matrix)
    -S: genealogically ordered pairs only will be computed (S matrix) [default]
    -B: (both) inseparable and genealogically ordered pairs will be computed (S+Id matrix)
    -A: (both) inseparable and adjacent pairs will be computed (A+Id matrix)
    -P: (both) inseparable and precedent pairs will be computed (P+Id matrix)
    -X: Exclusion matrix will be computed (X matrix)

 [filters on output pairs]
    -s <real value>: set p-value threshold to output on stdout a Newick tree with events in name (default: 0.05)
    -M: output also monomorphic event sites (i.e. couple of events that do not occur) [default not]
    -0: output also p-value of 0 inseparable/genealogically ordered pairs (p-value=1)
    -T <real value>: set p-value threshold to output for pairs (default: 1.0, i.e. show all)

 [output tree with events]
   -N: output on stdout a Newick tree with events (only pairs with pval < 0.05) [default not]

 [input and output operations]
   -E <events>: reads a 1-column file for displaying event IDs instead of default [default: e1,e2...]
   -O <output>: set a value for the output folder and files prefixes [default: out]
   -1 <event 1>: set the first event id (integer) to be tested
   -2 <event 2>: set the second event id (integer) to be tested

 [verbose options]
    -v: (verbose) outputs running progression of the program
    -V: (very verbose) outputs running progression of the program and more...
    -d: (distribution) output on stderr the inseparable and/or genealogically ordered pairs probabilities distribution
```

Epics outputs two files of the form `[prefix]_mat_epics_[matrix].tab` and `[prefix]_mat_signif_epics_[matrix]`. The second file is of the same organization of the first file but contain only the significant pairs. The information is organized as follows:

| IDe1 | IDe2 | Event1 | Event2 | Matrix | Pval | Nobs | Max |
| ---- | ---- | ------ | ------ | ------ | ---- | ---- | --- |
| first event in pair | second event in pair | user supplied e1 ID | user supplied e2 ID | pairs matrix | p-value | number of observed ordered pairs | maximum possible order of pairs |

### Epocs

Epocs maximizes the likelihood of independence and dependence models, with at most eight parameters. For details on the implementation and modelling, we refer the user to [the paper](https://doi.org/10.1093/sysbio/syab092). The standard approach is to run epocs multiple times on different models, then calculate the likelihood ratio tests between nested models and retrieve the model whose LRT is maximized while minimizing the number of parameters.

type `./epocs -h` to see the help and the main commands.

The standard format to run epocs is as follows

`./epocs [OPTIONS] "OUTGROUP" PHYLO_TREE`, where `OUTGROUP` is either empty when no outgroup is present (`epics [OPTIONS] "" PHYLO_TREE`) or takes the form `"name1[:name2:name3:...]"`.

```
[input options]
   -1 <event 1>: set the first event id (integer) to be tested
   -2 <event 2>: set the second event id (integer) to be tested
   -R <filename>: input filename containing significant events to run (epics format)

 [output options]
   -O: <output>: set a value for the output folder and files prefixes [default: out]
   -t: (tree) output the input tree on stderr
   -e: output on stderr the vector of states)
   -v: (verbose) outputs running progression of the program
   -V: (very verbose) outputs running progression of the program and more...

 [ML options]
   -C   : Maximum number of co-occurrences to consider in the analysis. Default is 10
   -I   : evaluate the Initial State (otherwise, assume <0,0>)
 [scenario selection]
   -S # : depending on the value of #, a different scenario is evaluated
      1 : 1 parameter: mu
      i : 2 parameters: mu1, mu2 (H0)
      x : 2 parameters: mu, lambda -- modified rates are mu*lambda --
      u : 3 parameters: mu1, nu1, mu2
      U : 3 parameters: mu1, mu2, nu2
      X : 3 parameters: mu, lambda, kappa --mu*lambda*kappa--
      r : 3 parameters: mu1, mu2, lambda (reciprocal induction: E1<->E2) -- modified rates are mu[12]*lambda --
      a : 3 parameters: mu1, mu2, mu2* (one-way induction: E1->E2)
      b : 3 parameters: mu1, mu1*, mu2 (one-way induction: E2->E1)
      l : 4 parameters: mu1, mu1*, mu2, mu2* (induction both ways)
      I : 4 parameters: mu1, mu2, nu1, nu2 (state dependance, no induction)
      R : 4 parameters: mu1, mu2, lambda, kappa (reciprocal induction: E1<->E2, and state dependance) --mu[12]*lambda*kappa--
      A : 5 parameters: mu1, mu2, mu2*, kappa1, kappa2 (one-way induction: E1->E2)
      B : 5 parameters: mu1, mu1*, mu2, kappa1, kappa2 (one-way induction: E2->E1)
      L : 6 parameters: mu1, mu1*, kappa1, mu2, mu2*, kappa2 (induction + state dependance)
 default: 8 parameters: mu1, mu1*, nu1, nu1*, mu2, mu2*, nu2, nu2*
```

epocs outputs one file per run, with the form `[prefix]_mat_epocs_[model][1-2].tab`. The number at the end reflects either a capitalized letter model (2) or not (1). This file naming is decided to work on case-insensitive filesystems such as MacOS. The content of the file is organized as follows"

| Event1 | Event2 | Tree | Model | lnML | ML | nu1 | nu1star | nu2 | nu2star | nu1 | nu1star | nu2 | nu2star |
| ------ | ------ | ---- | ----- | ---- | -- | --- | ------- | --- | ------- | --- | ------- | --- | ------- |
| first event in pair | second event in pair | tree ID | model | log-likelihood | likelihood | param | param | param | param | param | param | param | param |


### Epocs_mcmc

For a more precise description of the likelihood distribution, epocs_mcmc allows to perform likelihood surface exploration by a standard Metropolis-Hastings algorithm exploring parameter values bounded by the original tools (between 0 and 1000). In addition, epocs_mcmc explores the transition vectors as well. Every even sampling round, a new configuration is drawn and the new likelihood is tested as in the standard MH algorithm. The transition vectors are described in the original epocs paper (ref). In short, they summarize whether in co-occurrences, the event 1 is occurring before event 2 ('3') or the other way round ('4').

The input required is the same as `epics` and `epocs`. The options available are as follows, after typing `./epocs_mcmc -h`

```
[input options]
  -1 <event 1>: set the first event id (integer) to be tested
  -2 <event 2>: set the second event id (integer) to be tested
  -N <integer>: set the number of rounds for the MCMC chain. Default is 100,000
  -S <sampling>: set the intervals for sampling the MCMC chain. Default is sampling at every iteration.

[output options]
  -O: <output>: set a value for the output folder and files prefixes [default: out]
  -v: (verbose) outputs running progression of the program
  -V: (very verbose) outputs running progression of the program and more...
```

By default, without any selection, epocs_mcmc considers event 1 vs event 2 to analyze. You can select any length of the MCMC chain and the sampling intervals.

As output, epocs_mcmc writes a file of the form `[prefix]_epocs_mcmc.tab` containing the parameter values of all sampled chains. The format is compatible with Tracer (ref) if one is willing to plot the distributions and trace results.

| state | log-posterior | mu1 | mu1star | mu2 | mu2star | nu1 | nu1star | nu2 | nu2star | tvector |
| ----- | ------------- | --- | ------- | --- | ------- | --- | ------- | --- | ------- | ------- |
| sampling | double | param | param | param | param | param | param | param | param | transition vectors |

## Examples

### Epics and epocs

You can try epics and epocs on a test tree provided in TestData. The tree contains three WC pairs (see the [epocs paper](https://doi.org/10.1093/sysbio/syab092) for more details on the data set) for which the tests should show two significant associations.

```
./epics -O WCPairEpics -T 0.05 -I "AB682439" TestData/ThreeWCPairsOutput_ACR_tree.nwk
./epocs -O WCPairEpocs -Si "AB682439" TestData/ThreeWCPairsOutput_ACR_tree.nwk
```

Will give you two folders `WCPairEpics` and `WCPairEpocs`, containing the results of the two scripts.

For a complete epocs run and retrieval of the best models, the user must include at least two runs with at least the model of independence `i` and any model in `a, b, l, I, A, B, L`. Those models cover the basic inductions that can be retrieved between any two pairs on a tree (either 1 -> 2, 2 -> 1 or 1 <-> 2, with state-dependence as well with the capitalized models). Then, the Rscript can be run to parse the results.

```
./epocs -O WCPairEpocs -Si "AB682439" TestData/ThreeWCPairsOutput_ACR_tree.nwk
./epocs -O WCPairEpocs -Sa "AB682439" TestData/ThreeWCPairsOutput_ACR_tree.nwk
./epocs -O WCPairEpocs -Sb "AB682439" TestData/ThreeWCPairsOutput_ACR_tree.nwk
./epocs -O WCPairEpocs -Sl "AB682439" TestData/ThreeWCPairsOutput_ACR_tree.nwk
./epocs -O WCPairEpocs -SI "AB682439" TestData/ThreeWCPairsOutput_ACR_tree.nwk
./epocs -O WCPairEpocs -SA "AB682439" TestData/ThreeWCPairsOutput_ACR_tree.nwk
./epocs -O WCPairEpocs -SB "AB682439" TestData/ThreeWCPairsOutput_ACR_tree.nwk
./epocs -O WCPairEpocs -SL "AB682439" TestData/ThreeWCPairsOutput_ACR_tree.nwk
Rscript Scripts/ParseEpocsFinal_WithLambda.R WCPairEpocs
```
Outputting the final `[prefix]_BestModel.tab`

### epocs_mcmc

To run epocs_mcmc, select any pair of events and the mcmc options. It is recommended to select at least 10M runs and sampling every 50, which should give decent results. To verify convergence, autocorrelations and other statistics, we advise to use Tracer (ref).

```
./epocs_mcmc -O epocsMCMC_1v2 -1 1 -2 -2 -N 10000000 -S 50 "AB682439" TestData/ThreeWCPairsOutput_ACR_tree.nwk
```

Afterwards, we include a small Rscript to plot the parameter distributions, as well as a barplot summarizing the transition vectors at co-occurrences. The script takes in argument the epocs_mcmc output file and a prefix. For example:

```
Rscript Scripts/PlotMCMC.R epocsMCMC_1v2/epocsMCMC_1v2_epocs_mcmc.tab 1v2
```

which will output two pdf: `[prefix]_barplot_tvector_MCMC.pdf` and `[prefix]_histogram_parameter_distribution.pdf`.

### evo-scope

To run a full evo-scope analysis, you can try the test data in the TestData folder like:

```
./evoscope -t TestData/ThreeWCPairs_forTests.tab -T TestData/align_cured_UbyT_2_tree_rooted -f
```

Where here we want to have the full run.

## Questions?

If you have any questions, bugs to report or suggestions, feel free to post an issue.

## Authors

Maxime Godfroid

Guillaume Achaz

## References

Behdenna A et al. 2016. Testing for Independence between Evolutionary Processes. Syst Biol. 65:812–823. doi: 10.1093/sysbio/syw004.

Behdenna A, Godfroid M et al., 2022. A Minimal yet Flexible Likelihood Framework to Assess Correlated Evolution, Syst Biol, 71:823–838, doi: 10.1093/sysbio/syab092

Godfroid, M et al., 2022. Evo-Scope: Fully automated assessment of correlated evolution on phylogenetic trees, bioRxiv, 2022.12.08.519595; doi: https://doi.org/10.1101/2022.12.08.519595 

Ishikawa S et al., 2019. A Fast Likelihood Method to Reconstruct and Visualize Ancestral Scenarios, Mol Biol Evol, 36:2069–2085, doi:10.1093/molbev/msz131
