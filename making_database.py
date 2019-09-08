import os,subprocess
from pathlib import Path

class DATABASE:

    def __init__(self, testset, dataset):
        self.testset = testset
        self.dataset = dataset

    def removedup(self):
        small_content = []
        larget_content = {}
        with open(self.testset, "r") as smallFile:
            oneline = smallFile.readline()
            temp = []
            while oneline:
                if not oneline.strip("\n").strip():
                    oneline = smallFile.readline()
                    continue
                if oneline.startswith(">"):
                    if temp:
                        small_content.append("".join(temp))
                        temp = []
                else:
                    temp.append(oneline.strip('\n'))
                oneline = smallFile.readline()
        small_content.append("".join(temp))
        with open(self.dataset, "r") as largeFile:
            oneline = largeFile.readline()
            temp = []
            entire = []
            while oneline:
                if not oneline.strip("\n").strip():
                    oneline = largeFile.readline()
                    continue
                if oneline.startswith(">"):
                    if entire:
                        larget_content.update({"".join(temp): "".join(entire)})
                        temp = []
                        entire = []
                else:
                    temp.append(oneline.strip("\n"))
                entire.append(oneline)
                oneline = largeFile.readline()
        larget_content.update({"".join(temp): "".join(entire)})

        print("Before removing, the number of larger dataset is "+ str(len(larget_content)))

        for each in small_content:
            if each in larget_content:
                del larget_content[each]

        print("After removing duplicate records, the number of larger dataset is " + str(len(larget_content)))
        # windows path
        # clean_dataset = str(Path().absolute()) + "\database\blast_database"
        # if not os.path.exists(str(Path().absolute()) + "\database\"):
        #     os.makedirs(str(Path().absolute()) + "\database\")
        clean_dataset = str(Path().absolute())+"/database/blast_database"
        if not os.path.exists(str(Path().absolute()) + "/database/"):
            os.makedirs(str(Path().absolute()) + "/database/")
        with open(clean_dataset, "w") as result:
            result.writelines(larget_content.values())

        print("done!")
        return clean_dataset

    def makedb(self, input):
        # windows path
        # make_db_output = str(Path().absolute()) + "\database\blast_database"
        # blast_database_exe = str(Path().absolute()) + "\tools\ncbi-blast-2.9.0+\bin\makeblastdb"
        make_db_output= str(Path().absolute())+ "/database/blast_database"
        blast_database_exe=str(Path().absolute())+"/tools/ncbi-blast-2.9.0+/bin/makeblastdb"
        subprocess.call([blast_database_exe,"-in",input,"-title",make_db_output,"-dbtype","prot"])

    def compute(self):
        output = self.removedup()
        self.makedb(output)
        os.remove(output)

