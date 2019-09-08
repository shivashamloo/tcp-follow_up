import os, subprocess
import re
from pathlib import Path
from Bio.Align.Applications import MuscleCommandline

class MSA:

    def __init__(self,name,data):
        self.name=name
        self.data=data

    def compute(self):
        if self.name == "muscle":
            output=self.compute_muscle()
        elif self.name == "mafft":
            output=self.compute_mafft()
        elif self.name == "tm-coffee":
            output = self.compute_tmcoffee()
        # elif self.name == "AQUA":
        #     output=aqua(self.data)
        return output

    def compute_muscle(self):
        # window path
        # muscle_exe=str(Path().absolute())+"\tools\+muscle3.8.31_i86darwin64"
        # path_out=str(Path().absolute())+"\MUSCLE_result\"
        muscle_exe = str(Path().absolute())+"/tools/muscle3.8.31_i86darwin64"
        path_out = str(Path().absolute())+"/MUSCLE_result/"
        if not os.path.exists(path_out):
            os.mkdir(path_out)
        blast_fasta = os.walk(self.data, topdown=True)
        list = []
        for path,dir_list,file_list in blast_fasta:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print ('Invalid File Name :' + file_name)
            for file_name in list:
                input_muscle = path + file_name
                fasta_file_name = file_name[0:(len(file_name)-6)]
                output_muscle = path_out + fasta_file_name + '.fasta'
                muscle_cline = MuscleCommandline(muscle_exe, input=input_muscle, out=output_muscle)
                subprocess.call(str(muscle_cline), stdout=subprocess.PIPE, shell=True)

        return(path_out)

    def compute_mafft(self):
        # windows path
        # path_out = str(Path().absolute())+"\MAFFT_result\"
        # mafft_exe = str(Path().absolute())+"\tools\mafft-mac\mafft.bat"
        path_out = str(Path().absolute())+"/MAFFT_result/"
        mafft_exe = str(Path().absolute())+"/tools/mafft-mac/mafft.bat"
        if not os.path.exists(path_out):
            os.mkdir(path_out)
        c = os.walk(self.data, topdown=True)
        list = []
        for path,dir_list,file_list in c:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print ('Invalid File Name :' + file_name)
            for file_name in list:
                input_path = path + file_name
                fasta_file_name = file_name[0:(len(file_name)-6)]
                out_put_path = path_out + fasta_file_name + '.fasta'
                with open(out_put_path, 'w') as outputFile:
                    subprocess.call([mafft_exe, input_path], stdout=outputFile)
        return path_out

    def compute_tmcoffee(self):
        # windows path
        # tm_coffee_out = str(Path().absolute())+"\TM-COFFEE_result\"
        tm_coffee_out = str(Path().absolute())+"/TM-COFFEE_result/"
        if not os.path.exists(tm_coffee_out):
            os.mkdir(tm_coffee_out)
        # windows path
        # make_db_output= str(Path().absolute())+"\database\database"
        make_db_output = str(Path().absolute())+"/database/database"
        d = os.walk(self.data, topdown=True)

        tm_coffee_valid_file_list = []
        for path,dir_list,file_list in d:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    tm_coffee_valid_file_list.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)
            for file_name in tm_coffee_valid_file_list:
                tm_coffee_input_path = path + file_name
                tm_coffee_fasta_file_name = file_name[0:(len(file_name)-6)]
                tm_coffee_out_put_path = tm_coffee_out + tm_coffee_fasta_file_name + '.fasta'
                subprocess.call(["t_coffee", tm_coffee_input_path, "-outfile", tm_coffee_out_put_path, "-output", "fasta"
                                , "-mode", "psicoffee", "-blast_server", "LOCAL", "-protein_db", make_db_output])

        return tm_coffee_out