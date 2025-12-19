#!/bin/bash

TEMPLATE="crabConfig_full.py"

DATASETS=(
"QCD_Bin-PT-1000to1500_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-120to170_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-1500to2000_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-15to20_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-170to300_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-2000to2500_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-20to30_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-2500to3000_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-3000_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-300to470_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-30to50_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-470to600_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-50to80_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-600to800_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-800to1000_TuneCP5_13p6TeV_pythia8"
"QCD_Bin-PT-80to120_TuneCP5_13p6TeV_pythia8"
)

CFG_DIR="crab_cfgs"
rm -rf $CFG_DIR
mkdir $CFG_DIR

for DS in "${DATASETS[@]}"; do

    FULL_DATASET="/${DS}/RunIII2024Summer24MiniAODv6-150X_mcRun3_2024_realistic_v2-v2/MINIAODSIM"

    REQ="SVAnalysis_${DS}"
    OUTTAG="SVAnalysis_${DS}"
    CFG="$CFG_DIR/crab_cfg_${DS}.py"

    echo "Creating $CFG"

    sed \
        -e "s|REQUESTNAME|$REQ|" \
        -e "s|INPUTDATASET|$FULL_DATASET|" \
        -e "s|OUTPUTDATASETTAG|$OUTTAG|" \
        $TEMPLATE > $CFG

    echo "Submitting CRAB job: $REQ"
    crab submit -c $CFG

done
