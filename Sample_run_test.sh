#!/bin/bash

PT_ID="JH-2-002 JH-2-009 JH-2-055 JH-2-023 WU-356 WU-368 WU-436 JH-2-079 WU-487 WU-386 WU-561 WU-225" 
Working_path="/Volumes/lyu.yang/MAC_1/R/Project/14_MPNST_tumor_evolution"

# Iterate over the PT_IDs
for i in $PT_ID; do
  bash ./batch_running_2.sh "$i" "$Working_path"
done