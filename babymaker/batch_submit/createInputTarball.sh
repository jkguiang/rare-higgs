#! /bin/bash

# if [ -f package.tar.xz ]; then
#     rm ./package.tar.xz
# fi

tar -hcJf package.tar.xz ../CORE/*.h ../CORE/Tools/*.h ../CORE/jetcorr/*.h ../CORE/datasetinfo/*.h ../ScanChain_C.so ../mcTree_C.so ../doAll.py


