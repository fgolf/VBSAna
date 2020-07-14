#!/usr/bin/env bash

outdir=$1
mkdir -p $outdir

for infile in `find /hadoop/cms/store/user/namin/run2_mc2017/VBSWWHLep_M125_13TeV_amcatnloFXFX_madspin_pythia8_PRIVATE_RunIIFall17MiniAODv3-PU2017_2Apr2020_94X_mc2017_realistic_v14-ext1-v2_MINIAODSIM_CMS4_V10-02-05/ -regex ^.*\.root$` ;
    do outfile=${infile#${infile%merged*}} ;
    nice -n 15 ./main.exe $infile -o $outdir/$outfile
done

