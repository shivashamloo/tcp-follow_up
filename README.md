# tcp-follow_up
This code builds an HMM based on a training set and then classifies the test sets. The class can belong to one of the following seven categories:
```bash
  1: Amino Acid   2: Anion    3: Cation   4:Electron    5:Protein   6:Sugar   7:Other
```
## FOLDERS

### matrix

Contains the matrix for use by xdet

### tools
The packages needed to be installed in order to run the code.

## Mechanism
The tool begins by separating individual sequences from the test and training sets, placing each sequence in its own fasta-formatted intermediate file. An intermediate BLAST database is then created by removing all sequences found in the test data from the provided database. The provided training sequences are then BLASTed against the generated database. The BLAST results are then converted to fasta-format and passed through the chosen MSA, whose results are in turn passed through the SDS and finally the results from the SDS are used to generated an HMM. This HMM is used to classify the sequences from the test set.

## How to use

```bash
  python main.py training_set test_set database MSA_name SDS_name
```

## Result
The final results, including the computed specificity, sensitivity, accuracy and mcc, are outpit in a pdf file called ```finalresult.pdf```.

In addition, the results from each of intermediate steps will be saved in their respective folders, with the individual sequeneces stored in the ```test_sequences``` and ```train_sequences``` folders, the initial blast results stored in the ```blastxml``` folder, the derived fasta-results in the ```blastfasta``` folder, the MSA and SDS results in ```[MSA NAME]_result ```, ```[SDS NAME]_result``` and for xdet the final result will be in ```SDS_XDet_result``` folders specific to the type of MSA/SDS used and the HMM stored in the ```HMM\HMM_build``` folder and the final result will be in the ```output``` folder.

  
