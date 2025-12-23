#!/bin/bash

# Simple condor submission script for histogram making

BASE_DIR="/xrootd/store/user/taehee"
PILOT=0

# QCD directories
QCD_DIRS=(
    "QCD_Bin-PT-15to20_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-20to30_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-30to50_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-50to80_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-80to120_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-120to170_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-170to300_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-300to470_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-470to600_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-600to800_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-800to1000_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-1000to1500_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-1500to2000_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-2000to2500_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-2500to3000_TuneCP5_13p6TeV_pythia8"
    "QCD_Bin-PT-3000_TuneCP5_13p6TeV_pythia8"
)

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --pilot)
            PILOT=1
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--pilot]"
            exit 1
            ;;
    esac
done

echo "Creating condor submission files..."

# Create submission directory
SUBMIT_DIR="condor_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$SUBMIT_DIR"
mkdir -p "$SUBMIT_DIR/logs"

# Copy python script
cp drawHisto.py "$SUBMIT_DIR/"

# Create condor executable wrapper
cat > "$SUBMIT_DIR/run_job.sh" << 'EOF'
#!/bin/bash
INPUT_FILE=$1
XRD_FILE="${INPUT_FILE/\/xrootd/root://cms-xrdr.private.lo:2094//xrd}"
QCD_DIR=$2
OUTPUT_NUM=$3

echo "Processing: $XRD_FILE"
echo "QCD Dir: $QCD_DIR"
echo "Output number: $OUTPUT_NUM"

source /cvmfs/sft.cern.ch/lcg/views/LCG_106/x86_64-el9-gcc13-dbg/setup.sh

python3 drawHisto.py "$XRD_FILE" --dirname "$QCD_DIR" --output-number "$OUTPUT_NUM"
EOF

chmod +x "$SUBMIT_DIR/run_job.sh"

TOTAL_JOBS=0

# Scan directories and create separate submit file for each QCD directory
for qcd_dir in "${QCD_DIRS[@]}"; do
    full_dir="$BASE_DIR/$qcd_dir"
    
    if [ ! -d "$full_dir" ]; then
        continue
    fi
    
    echo "Scanning: $qcd_dir"
    
    # Create job list for this QCD directory
    JOBLIST="$SUBMIT_DIR/jobs_${qcd_dir}.txt"
    > "$JOBLIST"
    
    counter=1
    dir_jobs=0
    for input_file in $(find "$full_dir" -name "output_*.root" | sort); do
        echo "$input_file $qcd_dir $counter" >> "$JOBLIST"
        counter=$((counter + 1))
        dir_jobs=$((dir_jobs + 1))
        
        if [ $PILOT -eq 1 ]; then
            break
        fi
    done
    
    # Create separate submit file for this directory
    if [ $dir_jobs -gt 0 ]; then
        cat > "$SUBMIT_DIR/submit_${qcd_dir}.jds" << EOF
universe                = vanilla
executable              = run_job.sh
arguments               = \$(input_file) \$(qcd_dir) \$(output_num)

output                  = logs/\$(qcd_dir)_\$(output_num).out
error                   = logs/\$(qcd_dir)_\$(output_num).err
log                     = logs/\$(qcd_dir)_\$(output_num).log

should_transfer_files   = NO

accounting_group = group_cms

queue input_file,qcd_dir,output_num from jobs_${qcd_dir}.txt
EOF
        
        TOTAL_JOBS=$((TOTAL_JOBS + dir_jobs))
        echo "  -> $dir_jobs jobs"
    fi
    
    if [ $PILOT -eq 1 ]; then
        break
    fi
done

echo "========================================="
echo "Total jobs: $TOTAL_JOBS"
echo "Submit directory: $SUBMIT_DIR"
echo "========================================="

cd "$SUBMIT_DIR"

if [ $PILOT -eq 1 ]; then
    echo "PILOT MODE - Submitting 1 job"
    for jds_file in submit_*.jds; do
        condor_submit "$jds_file"
        break
    done
else
    echo "Submitting jobs for each QCD directory..."
    for jds_file in submit_*.jds; do
        qcd_name=$(echo "$jds_file" | sed 's/submit_//; s/.jds//')
        echo "Submitting: $qcd_name"
        condor_submit "$jds_file"
    done
    echo "========================================="
    echo "All jobs submitted! Monitor with: condor_q"
fi
