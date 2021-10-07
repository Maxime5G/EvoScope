# README

## How to install

To install both epics and epocs, simply type `make`. This will build two executables `epics` and `epocs`. You do not need any external libraries.

## Input

The only input required to run either epics or epocs is a rooted phylogenetic tree, in newick format, with mutations mapped on every branches. The mutations are are encapsulated between square brackets and located at the bootstrap values in the newick string for internal branches (e.g., `)[0/1/0/1/1]`) or after the tip ids (e.g., `SPECIES1[0/0/0/1/1/1]`). In this example, we have six traits, three of which mutate in SPECIES1.

## How to run

type `./epics -h` or `./epocs -h` to see the help and the main commands.

The standard format to run either epics or epocs with your data is as follows

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

### List of options for Epocs

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

## Example run

You can try epics and epocs on a test tree provided in TestData. The tree contains three WC pairs (see our paper: Behdenna, Godfroid et al., 2021, BioRxiv) for which the tests should show three significant associations.

```
./epics -O WCPairEpics -T 0.05 -I "" TestData/ThreeWCPairsOutput_ACR_tree.nwk
./epocs -O WCPairEpocs -Si "" TestData/ThreeWCPairsOutput_ACR_tree.nwk
```

Will give you two folders `WCPairEpics` and `WCPairEpocs`, containing the results of the two scripts.

## Output files

epics output files consist of 8 fields: the number of trait 1 (IDe1), the number of trait 2 (IDe2), the id of trait 1 (Event1), the id of trait 2 (Event2), the matrix used, the p-value, the number of associations (Nobs) (wrt the matrix selected) and the maximum possible number of observations (Max).

epocs output files consist of 14 fields: the number of trait 1 (IDe1), the number of trait 2 (IDe2), the tree, the model selected, the lnML, the ML, and then the values for the eight parameters of the model.

## Questions?

If you have any questions, bugs to report or suggestions, feel free to post an issue.

## Authors

Maxime Godfroid

Guillaume Achaz

## References

Behdenna A et al. 2016. Testing for Independence between Evolutionary Processes. Syst Biol. 65:812–823. doi: 10.1093/sysbio/syw004.

Behdenna A, Godfroid M et al. 2021. A minimal yet flexible likelihood framework to assess correlated evolution. bioRxiv. 2020.09.04.282954. doi: 10.1101/2020.09.04.282954.
