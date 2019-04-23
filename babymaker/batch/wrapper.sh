#!/bin/bash

# --> Arguments <-- #
# Input
FILEID=$1
SAMPLE=$2
COPYDIR=$3
# Log input arguments
echo "[wrapper] FILEID    = " ${FILEID}
echo "[wrapper] COPYDIR   = " ${COPYDIR}

# --> Environment <-- #
# Version using cvmfs install of CMSSW
CMSSW_VERSION=CMSSW_9_4_1
echo "[wrapper] setting env"
export SCRAM_ARCH=slc6_amd64_gcc630
source /cvmfs/cms.cern.ch/cmsset_default.sh
OLDDIR=`pwd`
cd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/$CMSSW_VERSION/src
# Set up cmssw env
eval `scramv1 runtime -sh`
cd $OLDDIR
# Log basic job info
echo "[wrapper] hostname  = " `hostname`
echo "[wrapper] date      = " `date`
echo "[wrapper] linux timestamp = " `date +%s`

# --> Input Decompression <-- #
# Log cwd
echo "[wrapper] current directory is:"
pwd
# Un-tar input tarball
echo "[wrapper] extracting input sandbox"
TARFILE=`ls *.tar.xz`
tar -xJf $TARFILE
# Log input contents
echo "[wrapper] input contents are"
ls -a

# --> Main Process <-- #
# Input file
INFILE="${SAMPLE}/merged_ntuple_${FILEID}.root"
echo "[wrapper] running: python doAll.py $INFILE" 
# Run main process
python doAll.py $INFILE

# --> Output <-- #
# Log output
echo "[wrapper] output is"
ls -lh
# Skim output
echo "[wrapper] skimming output"
python copyTree.py output.root skimmed.root -1 0 "tree" "recoMeson_nCands > 0"
# Mark output to be copied
echo "[wrapper] copying file"
OUTPUT=skimmed.root
echo "[wrapper] OUTPUT = " ${OUTPUT}

# --> Clean-up <-- #
# Set up destination directory
if [ ! -d "${COPYDIR}" ]; then
    echo "creating output directory " ${COPYDIR}
    mkdir ${COPYDIR}
fi
# Copy output to destination
export LD_PRELOAD=/usr/lib64/gfal2-plugins/libgfal_plugin_xrootd.so
gfal-copy -p -f -t 4200 --verbose file://`pwd`/${OUTPUT} gsiftp://gftp.t2.ucsd.edu${COPYDIR}/output_${FILEID}.root
# Clean up remaining files
echo "[wrapper] cleaning up"
for FILE in `find . -not -name "*stderr" -not -name "*stdout"`; do rm -rf $FILE; done
echo "[wrapper] cleaned up"
pwd
ls
