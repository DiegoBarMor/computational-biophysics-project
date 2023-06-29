#!/bin/bash

mkdir ../data_analysis
cd ../data_analysis
mkdir _coords _trajectories bse clustering pca pyinteraph rama rgyr rmsd rmsf sasa vmd wald

cd ../data_md
for run in mt1_rep0 mt1_rep1 mt2_rep0 mt2_rep1 wt1_rep0 wt1_rep1 wt2_rep0 wt2_rep1
do
    cd $run

    ############################################################################ POST-PROCESSING
    echo 1 0 | gmx trjconv -s 4_MD/md_plain.tpr -f 4_MD/md_plain.xtc -o md-nojump.xtc -pbc nojump -center

    echo 3 0 | gmx trjconv -s 4_MD/md_plain.tpr -f md-nojump.xtc -o md-center.xtc -pbc mol -center
    rm md-nojump.xtc

    echo 1 0 | gmx trjconv -s 4_MD/md_plain.tpr -f md-center.xtc -o md-rottrans.xtc -fit rot+trans
    rm md-center.xtc

    mv md-rottrans.xtc ../../data_analysis/_trajectories/$run.xtc
    cp 3_NPT/npt.gro ../../data_analysis/_trajectories/$run.gro

    cd ..
done
