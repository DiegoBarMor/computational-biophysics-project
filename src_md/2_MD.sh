#!/bin/bash

for run in mt1_rep0 mt1_rep1 mt2_rep0 mt2_rep1 wt1_rep0 wt1_rep1 wt2_rep0 wt2_rep1
do
    qsub "./"$run".pbs"
done
