export SCRAM_ARCH=slc5_amd64_gcc462
source /osg/app/cmssoft/cms/cmsset_default.sh
#source /apps/02/cmssoft/cms/cmsset_default.sh
source /osg/app/glite/etc/profile.d/grid_env.sh

echo ""
echo "----------------------------------------------------"
echo "Job started on `date` at WN: `hostname` "
echo "Job is running on `uname -a`"

cd /net/hisrv0001/home/rkunnawa/WORK/CMSSW_5_3_8_HI/
#eval `scram list CMSSW`
eval `scramv1 runtime -sh`
cd /net/hisrv0001/home/rkunnawa/WORK/CMSSW_5_3_8_HI/src/

echo "root directory: $ROOTSYS"

gcc --version

startfile=$1
endfile=$2
destination=$3
process=$4

echo "Processing..."

root -b -q merge_kurt_files_V3.C\+\($startfile,$endfile\)

hadoop dfs -moveFromLocal *.root /cms/store/user/rkunnawa/rootfiles/pPb/2013/data/purdue_forests/

echo "Done!"

echo "Copying output files to " $destination
