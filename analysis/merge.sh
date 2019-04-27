#!/usr/bin/env sh

basedir="/hadoop/cms/store/user/jguiang/rare-higgs/"

mkdir -p outputs/
for sname in `ls $basedir`; do
    echo $sname
    hadd -j 4 -k -f outputs/output_${sname}.root $basedir/$sname/*.root
done
