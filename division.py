from Bio import SeqIO
import os
from pathlib import Path

class DIVISION:

    def __init__(self,data, type):
        self.data = data
        self.type = type

    def divide(self):

        file_type = "fasta"
        if(self.type == "train"):
            # windows path
            # original_output = str(Path().absolute())+"\train_sequences\"
            original_output = str(Path().absolute())+"/train_sequences/"
        elif (self.type == "test"):
            # windows path
            # original_output = str(Path().absolute())+"\test_sequences\"
            original_output = str(Path().absolute()) + "/test_sequences/"

        if not os.path.exists(original_output):
            os.makedirs(original_output)

        for f, seq_record in enumerate(SeqIO.parse(self.data, file_type)):
            file_number = str(f)
            o_output = original_output + file_number + ".fasta"
            SeqIO.write(seq_record, o_output, "fasta")
