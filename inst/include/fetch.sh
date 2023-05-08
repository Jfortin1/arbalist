#!/bin/bash

set -e
set -u

# Fetches all of the relevant header files. We vendor it inside the package
# so that downstream packages can simply use LinkingTo to get access to them.

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
git checkout 0a19ec6d0b3b460e2e069cdcec6c7b7e8fda09eb
rm -rf ../tatami
cp -r include/tatami/ ../tatami
git checkout master
cd -

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
git checkout c638112d307eb873ec0c719046cf202d061dc4d0
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
git checkout 0a19ec6d0b3b460e2e069cdcec6c7b7e8fda09eb
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
git checkout 0eee555
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
git checkout 84d48b6
rm -rf ../aarand
cp -r include/aarand/ ../aarand
git checkout master
cd -

