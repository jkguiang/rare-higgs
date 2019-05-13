#!/usr/bin/env sh

basedir="/hadoop/cms/store/user/jguiang/rare-higgs"
outdir="/nfs-7/userdata/jguiang/rare-higgs"
ver="v3-0-0"
years=(2016 2017 2018)

mkdir -p outputs/
for year in "${years[@]}"; do
    for vname in `ls $basedir/$year`; do
        for sname in `ls $basedir/$year/$vname`; do
            if [[ "$vname" == "$ver" ]]; then
                echo "$sname $vname"
                copyTree.py "$basedir/$year/$ver/$sname/*.root" $outdir/$year/$ver/output_${sname}.root -1 0 tree "recoGamma_pt > 0 && recoWLepton_pt > 0" &> logs.txt &
            fi
        done
    done
done
