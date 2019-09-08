import subprocess
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import SeqIO
from Bio import SearchIO
import re
from os.path import isfile, join
import os
import sys
from pathlib import Path

class BlastRun:

    blast_run = str(Path().absolute())+"/blastxml/"
    blast_fasta_folder = str(Path().absolute())+"/blastfasta/"
    queries = str(Path().absolute())+"/train_sequences/"
    # windows path
    # blast_run = str(Path().absolute()) + "\blastxml\"
    # blast_fasta_folder = str(Path().absolute()) + "\blastfasta\"
    # queries = str(Path().absolute()) + "\train_sequences\"

    # blast run
    def run(self):
        # windows path
        # blast_exe = str(Path().absolute()) + "\tools\ncbi-blast-2.9.0+\bin\blastp"
        # make_db_output = str(Path().absolute()) + "\database\blast_database"
        blast_exe = str(Path().absolute())+"/tools/ncbi-blast-2.9.0+/bin/blastp"
        make_db_output = str(Path().absolute())+"/database/blast_database"

        blast_out_path = self.blast_run
        a = os.walk(self.queries)
        blast_valid_file_list_first = []
        for path,dir_list,file_list in a:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    blast_valid_file_list_first.append(file_name)
                else:
                    print ('Invalid File Name :' + file_name)
            for file_name in blast_valid_file_list_first:
                blast_query_path = path + file_name
                blast_name = file_name[0:(len(file_name)-6)]
                blast_xml = blast_out_path + blast_name + '.xml'
                blastp_cline = NcbiblastpCommandline(blast_exe, query=blast_query_path, db=make_db_output, evalue=0.001,
                                                     outfmt=5, out=blast_xml, max_target_seqs=10)
                subprocess.call(str(blastp_cline), stdout=subprocess.PIPE, shell=True)

    # blast xml to fasta
    def xml(self):
        xml_files = self.blast_run
        b = os.walk(xml_files)
        blast_valid_file_list = []
        for path, dir_list, file_list in b:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if result:
                    blast_valid_file_list.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)
            for file_name in blast_valid_file_list:
                blast_xml_file = path + file_name
                fasta_file_name = file_name[0:(len(file_name) - 4)]
                blast_fasta_file = self.blast_fasta_folder + fasta_file_name + '.fasta'
                blast_qresult = SearchIO.read(blast_xml_file, 'blast-xml')
                records = []
                for hit in blast_qresult:
                    records.append(hit[0].hit)
                SeqIO.write(records, blast_fasta_file, "fasta")
                del records[:]

    def compute(self):
        if not os.path.exists(self.blast_run):
            os.makedirs(self.blast_run)
        if not os.path.exists(self.blast_fasta_folder):
            os.makedirs(self.blast_fasta_folder)
        self.run()
        self.xml()
        return(self.blast_fasta_folder)




