# tcp-follow_up
This code builds an HMM based on a training set and then classifies the test sets. The class can belong to one of the following seven categories:
```bash
  1: Amino Acid   2: Anion    3: Cation   4:Electron    5:Protein   6:Sugar   7:Other
```
## FOLDERS

### metrix

Contains the metrix for the xdet

### tools
The packages needed to be installed in order to run the code.

## How to use

```bash
  training_set test_set database MSA_name SDS_name
```

## Result
the result of each step will be saved in the ```  MSA_name_result ```and ``` sds_name_result ``` folder.
  
