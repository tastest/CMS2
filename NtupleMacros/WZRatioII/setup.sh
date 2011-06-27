

DIR=`pwd`
DIR2=$PWD/../../
export CMS2_LOCATION=$DIR2

LOC=/Users/dlevans/store/
if [ `hostname` = "cms-tas03.cern.ch" ]; then
	LOC=/store/disk01/
fi

export CMS2_NTUPLE_LOCATION=$LOC

