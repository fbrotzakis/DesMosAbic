#!/bin/bash
i=$1
rm -r r$i
echo $i
mkdir r$i
cd r$i
cp ../TOPO/plumed.dat .
cp ../TOPO/*pdb .
gmx_mpi grompp -f ../TOPO/run.mdp -c ../TOPO/centered.gro -p ../TOPO/smog.top -o run.tpr -maxwarn 2
gmx_mpi mdrun -deffnm run -v -g --plumed plumed.dat -ntomp 1
cd -
