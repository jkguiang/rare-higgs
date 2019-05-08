#! /bin/bash

if [ -f package.tar.xz ]; then
    rm ./package.tar.xz
fi

tar -hcJf package.tar.xz \
    ../CORE/*.h ../CORE/Tools/*.h ../CORE/Tools/jetcorr/*.h ../CORE/Tools/datasetinfo/*.h ../CMS3_CORE.so \
    ../ScanChain.C ../mcTree.h ../mcTree.C ../magicAngles.h \
    ../copyTree.py ../doAll.py \
    ../jetCorrections/* \
