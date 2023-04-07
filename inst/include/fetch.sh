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
git checkout 2243b5ff46a8947adf4c95cb48a1c15e2226ccb6 
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

