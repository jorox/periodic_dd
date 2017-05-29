#!/bin/bash

fname=screw.data
currdir=$(pwd)

fdata=$currdir/$fname               #location of data file 
potfile=../potentials/Zr_3.eam.fs  #location of potential file (relative to input script)
ename=Zr                            #name of element in EAM file
emass=91.224
tmin=10000                          #number of minimization steps <100000
resdir=$currdir/_$fname              #results directory

mkdir -p $resdir

mpirun -np 20 lmp_mpi \
-var fdata $fdata \
-var potfile $potfile \
-var ename $ename \
-var tmin $tmin \
-var res $resdir \
-var emass $emass \
-log $resdir/log.$fname \
-in ../scripts/in.main_md 

