#!/bin/bash
# Writing a standard procedure from tip states to correlated evolution inference

# Steps are as follows:
# - states at the leaves in 1 file + phylogenetic tree w/ matching tip ids
# - pastml - JOINT
# - RSCRIPT: AnnotateATree_WithMasks_Clean.R
# - Epics I
# - Epocs i, a, b, l
# - RSCRIPT: ParseEpocsFinal_WithLambda.R

# FUNCTIONS

help() {
    cat << EOF
Usage: ./Pipeline.sh [options]
Options:
    -t: tab-delimited file of traits at the leaves [mandatory]
    -T: Phylogenetic tree in newick format [mandatory]
    -O: Outgroup [Optional, default = no outgroup. format= "name1[:name2:name3:...]"]
    -o: Output prefix [Optional, default = "out"]
    -m: epics matrix (default=Identity)
    -f: full run? (pastml-epics-epocs - default only until epics)
    -p: Plotting the results
    -h: Display this help message
EOF
}

verify_command() {
    command -v "$1" > /dev/null 2>&1
}

verify_rpackage() {
    R --slave -e "packageVersion('$1')" > /dev/null 2>&1
}

verify_all_dependencies() {
    if verify_command epics; then
        printf 'ok'
    else
        printf 'epics not installed'
    fi

    if verify_command epocs; then
        printf 'ok'
    else
        printf 'epocs not installed'
    fi

    if verify_command pastml; then
        printf 'ok'
    else
        printf 'pastml not installed'
    fi

    if verify_rpackage tidyverse; then
        printf 'ok'
    else
        printf 'tidyverse not installed'
    fi

    if verify_rpackage gplots; then
        printf 'ok'
    else
        printf 'gplots not installed'
    fi
}

# OPTIONS

tabfile=""
treefile=""
outputprefix=""
full=0
plot=0
scriptdir=`dirname "$0"`
epicsmatrix="I" # Matrix of chronology - by default is identity

while getopts t:T:O:o:m:fph opt; do
    case $opt in
        t)
            tabfile=$OPTARG
            ;;
        T)
            treefile=$OPTARG
            ;;
        o)
            outputprefix=$OPTARG
            ;;
        O)
            outgroup=$OPTARG
            ;;
        f)
            full=1
            ;;
        m)
            epicsmatrix=$OPTARG
            ;;
        p)
            plot=1
            ;;
        h)
            help; exit 0;
        esac
    done

# --- CHECKING COMMAND LINE ARGUMENTS ---

checking_command_line_args() {

    if [[ $tabfile == "" ]]; then
        printf "WARNING: Missing character table, exiting...\n\n";
        help;
        exit 1;
    fi

    if [[ $treefile == "" ]]; then
        printf "WARNING: Missing phylogenetic tree, exiting...\n\n";
        help;
        exit 1;
    fi

    if [ ! -e "$treefile" ]; then
        printf "Tree '$treefile' does not exist! Exiting...\n\n"
        help;
        exit 1;
    fi

    if [ ! -e "$tabfile" ]; then
        printf "File '$tabfile' does not exist! Exiting...\n\n"
        help;
        exit 1;
    fi

    if [[ $outputprefix == "" ]]; then
        printf "Setting the output prefix to default (out)\n\n";
        outputprefix="out";
    fi

}

# --- MAIN WRAPPERS ---

pastml_wrapper() {

    # Here: checking input tree - removing dots (better because easier to deal with the output name from pastml)
    # Then running pastml

    treefilenodot=`echo $treefile | sed 's#\.#_#g'`
    cp $treefile $treefilenodot
    basenametreefile=$(echo ${treefilenodot##*/})

    # --- ANCESTRAL CHARACTER RECONSTRUCTION ---
    # --- Tool: PastML, Model: Joint ---

    nthread=1       # Number of threads to use - default: 1

    outfolder=$outputprefix"_ACRFOLDER"
    outtab=$outputprefix"_ACR.tab"
    namedtree=$outfolder"/named.tree_"$basenametreefile".nwk"

    # Verbose
    printf 'Running pastml...\n'
    printf 'Command: pastml -t '$treefilenodot' -d '$tabfile' -o '$outtab' --prediction_method JOINT --forced_joint --work_dir '$outfolder' --threads '$nthread' --offline\n\n'

    # Actual command
    pastml -t $treefilenodot -d $tabfile -o $outtab --prediction_method JOINT --forced_joint --work_dir $outfolder --threads $nthread --offline

}

rscript1_wrapper(){

    # --- RECONSTRUCTING THE MUTATIONS ON THE TREE ---
    # --- Rscript: AnnotateATreeForEpoics.R ---
    # --- Arguments: - Named tree from pastml
    #                - Tab-delimited file from pastml
    #                - Output phylogenetic tree with mutations reconstructed
    #                - Output tab-delimited file with mutations reconstructed

    # Extracting the event IDs for epics afterwards

    cleanACR=$outputprefix"_ACR_clean.tab"
    matACR=$outputprefix"_ACR_clean.mat"
    epoicsTree=$outputprefix"_ACR_tree.nwk"

    # Verbose
    printf 'Reconstructing the mutations from the ACR...\n'
    printf 'Command: Rscript '$scriptdir'/Scripts/AnnotateATree_WithMasks_Clean.R '$namedtree' '$outtab' '$epoicsTree' '$matACR'\n\n'

    # Actual command
    Rscript $scriptdir/Scripts/AnnotateATree_WithMasks_Clean.R $namedtree $outtab $epoicsTree $matACR

    # --- EXTRACTING THE TRAIT IDS ---
    # --- Instead of displaying e1, e2 etc in epics
    # --- Display the trait ids that are the column ids in the input table
    # --- Saving the ids in a 1-column file "eventIDs.in"

    numbevents=$(awk '{print NF}' $matACR | head -n1)
    eventID=$outputprefix"eventIDs.in"
    cat $matACR | head -n1 | tr '\t' '\n' | tail -n `awk '{print NF-1}' $matACR | head -n1` > $eventID

}

epics_wrapper(){

    # --- EPICS OPTIONS ---
    # --- VERSION 1.0: - Matrix identity
    #                  - All-vs-all

    epicsevents=0   # 0: all vs all, N [1-nevt]: GWAS-like analysis (1 vs all - WIP)

    epoicsOutput=$outputprefix"_correlations"

    # Verbose
    printf 'Running epics...\n'
    if [[ $epicsevents -eq 0 ]]
    then
        # Verbose
        printf 'Command: '$scriptdir'/EpoicsMain/epics -E '$eventID' -T 0.05 -O '$epoicsOutput' -'$epicsmatrix' "'$outgroup'" '$epoicsTree'\n\n'

        # Actual command
        epics -1 $epicsevents -E $eventID -T 0.05 -O $epoicsOutput -$epicsmatrix "$outgroup" $epoicsTree
    else
        # Verbose
        printf 'Command: '$scriptdir'/EpoicsMain/epics -1 '$epicsevents' -E '$eventID' -T 0.05 -O '$epoicsOutput' -'$epicsmatrix' "'$outgroup'" '$epoicsTree'\n\n'

        # Actual command
        epics -1 $epicsevents -E $eventID -T 0.05 -O $epoicsOutput -$epicsmatrix "$outgroup" $epoicsTree
    fi

    # --- PLOT 1: QQPLOT ---
    if [[ $plot -eq 1 ]]
    then
        Rscript $scriptdir/Scripts/PlotUtilities.R $epoicsOutput QQ
    fi
}

epocs_wrapper(){

    # --- EPOCS OPTIONS ---
    # --- INPUT FILE: epics significant pairs
    # --- ARGUMENTS: Models to test are: -i [independence model - no state-dependent rates]
    #                                    -a [induction 1 -> 2 - no state-dependent rates]
    #                                    -b [induction 2 -> 1 - no state-dependent rates]
    #                                    -l [reciprocal induction 1 <-> 2 - no state-dependent rates]

    basenameepicsfile=$(echo ${epoicsOutput##*/})
    epicsSignif=$epoicsOutput"/"$basenameepicsfile"_mat_signif_epics_"$epicsmatrix".tab"

    # Verbose
    printf 'Running epocs...\n'
    printf 'Commands: epocs -R '$epicsSignif' -O '$epoicsOutput' -Si "'$outgroup'" '$epoicsTree'\n'
    printf '          epocs -R '$epicsSignif' -O '$epoicsOutput' -Sa "'$outgroup'" '$epoicsTree'\n'
    printf '          epocs -R '$epicsSignif' -O '$epoicsOutput' -Sb "'$outgroup'" '$epoicsTree'\n'
    printf '          epocs -R '$epicsSignif' -O '$epoicsOutput' -Sl "'$outgroup'" '$epoicsTree'\n\n'

    # Actual commands
    epocs -R $epicsSignif -O $epoicsOutput -Si "$outgroup" $epoicsTree
    epocs -R $epicsSignif -O $epoicsOutput -Sa "$outgroup" $epoicsTree
    epocs -R $epicsSignif -O $epoicsOutput -Sb "$outgroup" $epoicsTree
    epocs -R $epicsSignif -O $epoicsOutput -Sl "$outgroup" $epoicsTree
}

rscript2_wrapper(){
    # --- RETRIEVING THE BEST MODEL ---
    # --- RSCRIPT: ParseEpocsFinal_WithLambda.R ---
    # --- Arguments: The output prefix (i.e., folder where all files are saved)
    #                The event ids (optional - but in the pipeline should work)

    # Verbose
    printf 'Retrieving the best model from epocs\n'
    printf 'Command: Rscript '$scriptdir'/Scripts/ParseEpocsFinal_WithLambda.R '$epoicsOutput' '$eventID'\n\n'

    # Actual command
    Rscript $scriptdir/Scripts/ParseEpocsFinal_WithLambda.R $epoicsOutput $eventID

    # --- PLOT 2: HEATMAP ---

    if [[ $plot -eq 1 ]]
    then
        Rscript $scriptdir/Scripts/PlotUtilities.R $epoicsOutput HEATMAP
    fi

    # --- OUTPUT 3: CYTOSCAPE ---

    if [[ $plot -eq 1 ]]
    then
        epocsBestModel=$epoicsOutput"/"$basenameepicsfile"_BestModel.tab"
        cytoscapeFile=$epoicsOutput"/"$basenameepicsfile"_Cytoscape.tab"
        printf "source\ttarget\tinteraction\tWeight\n" > $cytoscapeFile

        if [ ! -e "$eventID" ]; then
            cut -f 1,2,7,9 $epocsBestModel | tail -n `awk 'END{print NR-1}' $epocsBestModel` >> $cytoscapeFile
        fi

        if [ -e "$eventID" ]; then
            cut -f 1,2,7,11 $epocsBestModel | tail -n `awk 'END{print NR-1}' $epocsBestModel` >> $cytoscapeFile
        fi
    fi
}

# MAIN

main(){

    # verifying the inputs from the user
    checking_command_line_args

    # Running pastml
    pastml_wrapper

    # Reconstructing the mutations on the branches of the tree
    rscript1_wrapper

    # Running epics
    epics_wrapper

    if [[ $full -eq 0 ]]
    then
        printf 'End of pipeline at epics!\n'
        exit 1;
    fi

    # Running epocs
    epocs_wrapper

    # Parsing the epocs results
    rscript2_wrapper

}

main