#!/bin/bash

# Check if Patient_ID is passed as an argument
if [ -z "$1" ]; then
  echo "Error: Please provide a Patient_ID as the first argument."
  exit 1
fi

if [ -z "$2" ]; then
  echo "Error: Please provide a working path as the first argument."
  exit 1
fi


# Set the Patient_ID from the command-line argument
Patient_ID="$1"
echo "$Patient_ID"

# Set the Working dir from the command-line argument
Path="$2"
echo "$Path"


mkdir -p $Path/Rscript
# define your github downloaded folder
Github_path="/Volumes/lyu.yang/MAC_1/Github/MPNST_tumor_evolution/Rcode"
cp -rf $Github_path/* $Path/Rscript

cd $Path/Rscript
TMP_FILE=$Path/Rscript/Step2_config.yml
echo $TMP_FILE
 
 # Update the Step2_config.yml file with the Patient_ID, Path, and Password
echo "Updating PT_ID in $TMP_FILE to: $Patient_ID"
sed -i.bak "s/PT_ID:.*/PT_ID: \"$Patient_ID\"/" "$TMP_FILE"



echo "Updating Path in $TMP_FILE to: $Path"
sed -i.bak "s/Path:.*|Path: \"$Path\"|" "$TMP_FILE"

# Verify that PT_ID, Path, and Password have been updated
grep -E "PT_ID|Path" "$TMP_FILE"



Processed_path=$Path/Processed/$Patient_ID
pSubTitle="${Patient_ID}_data_curation.pdf" 
base_output_dir="$Path/Report/$Patient_ID"

Rscript -e "rmarkdown::render('$Path/Rscript/1_step1_RD_BAF_combination.Rmd', output_file='$pSubTitle', output_dir='$base_output_dir')"

echo "step1 finished"

# Find the total number of *.clonality.rdata files in the current working directory
#Total_samples=$(ls "$Processed_path"/*clonality.rdata 2>/dev/null | wc -l)
Total_samples=$(cat $Processed_path/n_sample.txt)
echo $Total_samples


# Check if there are any .clonality.rdata files
if [ "$Total_samples" -eq 0 ]; then
  echo "Error: no sample numbers was provided."
  exit 1
fi


echo "Total number of samples : $Total_samples"

# Loop through each Sample_input value from 1 to Total_samples
for Sample_input in $(seq 1 $Total_samples); do
    echo "Updating Sample_input in $TMP_FILE to: $Sample_input"
    
    # Use sed to update the Sample_input in the temporary YML file
    sed -i.bak "s/Sample_input:.*/Sample_input: \"$Sample_input\"/" "$TMP_FILE"
    
    # Print the updated Sample_input for verification
    grep "Sample_input" "$TMP_FILE"
  
    
    # Run the R Markdown file using rmarkdown::render
    
    pSubTitle="${Patient_ID}_${Sample_input}_CNV.pdf"
    echo "$pSubTitle"
    echo "$base_output_dir"
    Rscript -e "rmarkdown::render('$Path/Rscript/1_step2_CNV.Rmd', params=list(Sample_input='$Sample_input', PT_ID='$Patient_ID'), output_file='$pSubTitle', output_dir='$base_output_dir')"
  
    # Optional: Wait for the current iteration to finish before continuing
    wait
done


# MTAP_CDKN2A
for Sample_input in $(seq 1 $Total_samples); do
  pSubTitle="${Patient_ID}_${Sample_input}__MTAP_CDKN2A_Loss.pdf"
  output_dir="$Path/Report/MTAP_CDKA2A/$Patient_ID"
  echo "Updating Sample_input in $TMP_FILE to: $Sample_input"
    
    # Use sed to update the Sample_input in the temporary YML file
    sed -i.bak "s/Sample_input:.*/Sample_input: \"$Sample_input\"/" "$TMP_FILE"
    
    # Print the updated Sample_input for verification
    grep "Sample_input" "$TMP_FILE"
    Rscript -e "rmarkdown::render('$Path/Rscript/3_MTAP_CDKN2A_loss.Rmd', params=list(Sample_input='$Sample_input', PT_ID='$Patient_ID'), output_file='$pSubTitle', output_dir='$output_dir')"

    wait
  done