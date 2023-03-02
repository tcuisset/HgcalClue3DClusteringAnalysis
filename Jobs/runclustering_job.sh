#!/bin/bash

export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

# if ! cd "$ANALYSIS_DIR" && `scram runtime -sh` && cd - ; then
#     >&2 echo "ERROR : could not cd to analysis path to run cmsenv"
# fi

cd "$ANALYSIS_DIR"
eval `scram runtime -sh`
cd -

#Create output file folder
mkdir -p -- $(dirname "$2")


./runclustering "$@"
