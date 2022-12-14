# PRIMITI-miRNA-target-mRNA-interactions

## Dependency

Conda is needed to be installed in order to create two conda environments that contain dependencies.
  1) PRIMITI environments - (PRIMITI_requirements.yml)
  2) PRIMITI_iLearn environments - (iLearn_requirements.yml)
  
Optional:
        If conda is not installed in the main root, please modify the variable MAIN_ENV_PATH and SUB_ENV_PATH in a config.py file to change to               environments path.

## Usage

python PRIMITI.py input_miRNA_csv input_transcript_csv job_id

positional arguments:
  1. input_miRNA_csv             choose csv file contain a list of miRNAs
  2. input_transcript_csv        choose csv file contain a list of transcripts
  3. job_id                      select job_id for outputs - the results with job_id as header can be found in the result folder
  
optional arguments:
  -h, --help                  show a help message and exit

## Instruction

Predict miRNA-target mRNA repression activities 

### Step 1: Create environments

Install two environments using conda - PRIMITI and PRIMITI_iLearn - with dependency installed.

```
conda env create -f PRIMITI_requirements.yml
conda env create -f iLearn_requirements.yml
```

### Step 2: Download required code data

Due to the excessive size of data required in PRIMITI.py, we are unable to upload them to github. 
Please download the code folder from google drive (https://drive.google.com/drive/folders/1UZ8mFEX1BwNUTKR1nOkaezpedAFd4nOQ?usp=sharing).

Put the code folder onto the same directoty as PRIMITI.py
If you have any inquery, please contact me through this email - k.uthayopas@uq.edu.au.

### Step 3: Run PRIMITI.py

run PRIMITI.py
```
python PRIMITI.py input_miRNA_csv input_transcript_csv job_id
```

### Step 4: Investigate the Result

The result can be found in a result folder. 
They include
  1) miRNA-target mRNA interaction table
  2) miRNA-target site binding table
  3) graphical representation for binding interaction (found in a subfolder with job_id header)
  4) error report (found in errors folder)

## Dataset 

Training dataset can be found in a dataset folder.

https://drive.google.com/drive/folders/1UZ8mFEX1BwNUTKR1nOkaezpedAFd4nOQ?usp=sharing
