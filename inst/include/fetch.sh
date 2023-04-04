#!/bin/bash

set -e
set -u

# Fetches all of the header files from tatami. We vendor it inside the package
# so that downstream packages can simply use LinkingTo to get access to them.

if [ ! -e source-tatami ]
then 
    git clone https://github.com/LTLA/tatami source-tatami
else 
    cd source-tatami
    git pull
    cd -
fi

cd source-tatami
git checkout 488a14ad9fabc62c6cf88d79aab8d34bb5778a4a
rm -rf ../tatami
cp -r include/tatami/ ../tatami
git checkout master
cd -
