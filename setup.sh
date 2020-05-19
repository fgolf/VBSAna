therelease=CMSSW_9_4_9
export SCRAM_ARCH=slc6_amd64_gcc700
mkdir -p common
if [ ! -d common/$therelease ]; then 
    cd common/ ;
    cmsrel $therelease;
    cd ~-
fi
cd common/$therelease/src
eval `scram ru -sh`
cd -
#

export ANABASE="$( cd "$(dirname "$BASH_SOURCE")" ; pwd -P )"

#echo $ANABASE
#export LD_LIBRARY_PATH=${ANABASE}/babymaking/batch/:$LD_LIBRARY_PATH
export PYTHONPATH=${ANABASE}/:$PYTHONPATH
export PYTHONPATH=${ANABASE}/common/matplottery:$PYTHONPATH

# bash PS1 gets mangled because CMSSW overrides this to LANG=C :(, so switch it back.
export LANG=en_US.UTF-8

# export PYTHONPATH=$PWD/analysis/bdt/xgboost/python-package/lib/python2.7/site-packages/:$PYTHONPATH
# export PYTHONPATH=$PWD/analysis/bdt/xgboost/python-package/:$PYTHONPATH
# export PYTHONPATH=$PWD/analysis/bdt/root_numpy-4.7.2/lib/python2.7/site-packages/:$PYTHONPATH

[[ -d ${ANABASE}/common/matplottery/ ]] || {
    git clone git@github.com:shchauha/matplottery.git ${ANABASE}/common/matplottery/;
    pip install --user matplotlib
    pip install --user uproot
    export PATH=/home/users/${USER}/.local/bin:$PATH
}

# checkout baby making code
#[[ -d ${ANABASE}/FTAnalysis/ ]] || git clone git@github.com:cmstas/FTAnalysis.git ${ANABASE}/FTAnalysis

[[ -d ${ANABASE}/babymaking/batch/NtupleTools/ ]] || git clone git@github.com:cmstas/NtupleTools.git ${ANABASE}/babymaking/batch/NtupleTools/

if [ ! -d ${ANABASE}/babymaking/batch/ProjectMetis/ ]; then
    git clone https://github.com/aminnj/ProjectMetis/ ${ANABASE}/babymaking/batch/ProjectMetis/
else
    source ${ANABASE}/babymaking/batch/ProjectMetis/setup.sh
fi
## [[ -d ${ANABASE}/common/Software/ ]] || git clone git@github.com:cmstas/Software.git ${ANABASE}/common/Software/
[[ -d ${ANABASE}/common/CORE/ ]] || {
    git clone git@github.com:cmstas/CORE.git ${ANABASE}/common/CORE/;
    cd ${ANABASE}/common/CORE; 
    make -j10 >&/dev/null &
    cd ${ANABASE};
}

[[ -d ${ANABASE}/common/${therelease}/src/HiggsAnalysis/CombinedLimit/ ]] || {
    pushd ${ANABASE}/
    cd $CMSSW_BASE/src
    cmsenv
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cd HiggsAnalysis/CombinedLimit
    git checkout 94x
    scramv1 b vclean
    scramv1 b -j 15
    popd
}
#
