#!/bin/bash

set -e
set -u

# Fetches all of the relevant header files. We vendor it inside the package
# so that downstream packages can simply use LinkingTo to get access to them.

##########################################################

if [ ! -e source-byteme ]
then 
    git clone https://github.com/LTLA/byteme source-byteme
else 
    cd source-tatami
    git pull
    cd -
fi

cd source-byteme
git checkout ad206e3
rm -rf ../byteme
cp -r include/byteme/ ../byteme
git checkout master
cd -

##########################################################

if [ ! -e source-tatami ]
then 
    git clone https://github.com/LTLA/tatami source-tatami
else 
    cd source-tatami
    git pull
    cd -
fi

cd source-tatami
git checkout 80f7ba4
rm -rf ../tatami
cp -r include/tatami/ ../tatami
git checkout master
cd -

##########################################################

if [ ! -e source-irlba ]
then 
    git clone https://github.com/LTLA/CppIrlba source-irlba
else 
    cd source-irlba
    git pull
    cd -
fi

cd source-irlba
git checkout 89e6c2f
rm -rf ../irlba
cp -r include/irlba/ ../irlba
git checkout master
cd -

##########################################################

if [ ! -e source-aarand ]
then 
    git clone https://github.com/LTLA/aarand source-aarand
else 
    cd source-aarand
    git pull
    cd -
fi

cd source-aarand
git checkout 4e41649
rm -rf ../aarand
cp -r include/aarand/ ../aarand
git checkout master
cd -

##########################################################
#
#if [ ! -e source-tatami_mult ]
#then 
#    git clone https://github.com/tatami-inc/tatami_mult source-tatami_mult
#else 
#    cd source-tatami_mult
#    git pull
#    cd -
#fi
#
#cd source-tatami_mult
#git checkout d924315
#rm -rf ../tatami_mult
#cp -r include/tatami_mult/ ../tatami_mult
#git checkout master
#cd -