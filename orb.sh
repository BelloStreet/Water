#/bin/bash

######print the coordinates for each center and how many are them
awk 'NF' "$1">temp.molden

awk '/Atoms/{flag=1; next} /GTO/{flag=0} flag' temp.molden>center.inp
ncenter=$(wc -l < center.inp)
echo ${ncenter}>read.inp

#######print the gaussian functions and the order of appearence of each center
awk '/GTO/{flag=1; next} /MO/{flag=0} flag' temp.molden>gaussian.inp
for ((i=1; i <=ncenter; i++))
do
    n=$(awk -v x="${i} 0" '$0~x {print NR}' gaussian.inp)
    echo ${n}>>read.inp
done

nlines=$(wc -l < gaussian.inp)
echo ${nlines}>>read.inp

#######printf coefficients of each molecular orbital
#######counts how many orbitals are in the file and how many of them are occupied
awk '/Sym=/,/$0/' temp.molden>molecule.inp
i=$(grep 'Occup' molecule.inp | wc -l)
j=$(grep 'Occup=    0.000000' molecule.inp | wc -l)
num=$((i - j))
echo ${i}>>read.inp
echo ${num}>>read.inp

