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