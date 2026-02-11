# SVAnalyzer
NanoAOD samples have some ambiguity on SV information.
Therefore, SVAnalyzer makes Ntuples from MiniAOD samples, which can be used for the further analysis.

## Initial setup
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_X_Y_Z
# 15_0_4: CMSSW used for the 2024 MiniAOD sample (Dec 2025)
# 15_1_0: has function KalmanVertexFitter::KalmanVertexFitter(bool useSmoothing, bool useMuonSystemBounds) {
cd CMSSW_X_Y_Z/src
cmsenv
mkdir MyAnalyzer && cd MyAnalyzer
git clone <this_repository>
```
## Build the analyzer
```
cd CMSSW_X_Y_Z/src
vi MyAnaylzer/SVAnalyzer/plugins/SVAnalyzer.cc #Edit the analyzer here
scram b -j 4

# The message below shows up when the build is successful.
# >> Entering Package MyAnaylzer/SVAnalyzer
# >> Creating project symlinks
# >> Leaving Package MyAnaylzer/SVAnalyzer

# Test the analzyer locally before submitting the CRAB jobs
cd MyAnaylzer/SVAnalyzer/test
cmsRun test_cfg.py
```

## Submit the CRAB jobs
`test/crab_cfg.py` is the CRAB configuration file where you should specify the object to use.
```
voms-proxy-init --voms cms
cd MyAnaylzer/SVAnalyzer/test
vi test/crabConfig_test.py #Specify the object that you use

# Test the CRAB jobs before the full submission
crab submit -c crabConfig_test.py

# Full CRAB submission utilizing crabConfig_full.py
# The dataset should be declared here
./submit_crab.sh
```

## Some useful commmands for CRAB jobs
```
# Check the status
crab status -d crab_projects/<requestName>
# Resubmit failed jobs
crab resubmit -d crab_projects/<requestName>
# Kill the job
crab kill -d crab_projects/<requestName>
# Get log files for debugging
crab getlog -d crab_projects/<requestName> --jobids=1
# Check specific job output
crab getoutput -d crab_projects/<requestName> --jobids=1,2,3
```

## Plot the results from Ntuple files
```
cd MyAnaylzer/SVAnalyzer/plot

# Make a root file containing histograms for each Ntuple file
vi drawHisto.py

# Test a condor job for drawing histograms
./submit_condor.sh --pilot

# Full condor submission
./submit_condor.sh

# The results will be saved in condor_<date>_<time> directory
# hadd all the files by using hadd.sh
```
