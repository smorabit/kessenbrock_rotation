#!/bin/bash
#$ -N transition_prob
#$ -q free64,pub64
#$ -m beas
#$ -ckpt restart
#$ -pe openmp 8

module load anaconda/3.6-5.0.1
python3 ~/bin/transition_prob.py
