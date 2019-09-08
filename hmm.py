import subprocess
import re
from os.path import join, isfile
import os
from pathlib import Path


class HMM:
    def __init__(self, data):
        self.data=data
    # windows path
    # hmm_out = str(Path().absolute()) + "\hmm\"
    hmm_out= str(Path().absolute())+"/hmm/"

    # hmmbuild
    def hmm_build(self):
        # windows path
        # path_out = self.hmm_out + "\hmm_build\"
        path_out = self.hmm_out+"/hmm_build/"
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
                file_name = file_name[0:(len(file_name)-6)]
                out_put_path = path_out + file_name
                # windows path
                # hmmbuild_exe = str(Path().absolute()) + "\tools\hmmer-3.2.1\src\hmmbuild"
                hmmbuild_exe = str(Path().absolute())+"/tools/hmmer-3.2.1/src/hmmbuild"
                subprocess.call([hmmbuild_exe, out_put_path, input_path])

    #hmmpress
    def hmm_press(self):
        # windows path
        # c = os.walk(self.hmm_out+"\hmm_build\", topdown=True)

        c = os.walk(self.hmm_out+"/hmm_build/", topdown=True)
        list = []
        for path,dir_list,file_list in c:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if result:
                    list.append(file_name)
                else:
                    print ('Invalid File Name :' + file_name)
            for file_name in list:
                input_path = path + file_name
                # windows path
                # hmmpress_exe=str(Path().absolute())+"\tools\hmmer-3.2.1\src\hmmpress"
                hmmpress_exe = str(Path().absolute())+"/tools/hmmer-3.2.1/src/hmmpress"
                subprocess.call([hmmpress_exe, input_path])

    def hmm_scan(self):
        # windows path
        # path_out = self.hmm_out + "\hmm_build\"
        # path_ou2 = str(Path().absolute()) + "\output\"
        # if not os.path.exists(path_ou2):
        #     os.mkdir(path_ou2)
        # a = os.walk(str(Path().absolute()) + "\test_sequences\")
        path_out = self.hmm_out+"/hmm_build/"
        path_ou2 = str(Path().absolute())+"/output/"
        if not os.path.exists(path_ou2):
            os.mkdir(path_ou2)
        a = os.walk(str(Path().absolute())+"/test_sequences/")

        list = []

        for path,dir_list,file_list in a:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if (result):
                    list.append(file_name)
                else:
                    print ('Invalid File Name :' + file_name)
            for file_name_2 in list:
                file = path + file_name_2
                path_out_2 = path_ou2 + file_name_2
                os.mkdir(path_out_2)
                pathout2 = path_out_2 + "/"
                for i in range(1,781,1):
                    # windows path
                    # hmmscan_exe=str(Path().absolute())+"\tools\hmmer-3.2.1\src\hmmscan"
                    hmmscan_exe = str(Path().absolute())+"/tools/hmmer-3.2.1/src/hmmscan"
                    path_3 = path_out + str(i)
                    pathout3 = pathout2 + str(i)
                    subprocess.call([hmmscan_exe, "-o", pathout3, "-E", "0.1", path_3, file,])

    def get_folders(self):
        # windows path
        # return[join(str(Path().absolute()) + "\output\", each_file) for each_file in os.listdir(str(Path().absolute()) + "\output\") if not isfile(join(str(Path().absolute()) + "\output\", each_file)) and not each_file.startswith('.')]
        return [join(str(Path().absolute())+"/output/", each_file) for each_file in os.listdir(str(Path().absolute())+"/output/")
                if not isfile(join(str(Path().absolute())+"/output/", each_file)) and not each_file.startswith('.')]

    def get_files(self,basic_path):
        return [join(basic_path, each_file) for each_file in os.listdir(basic_path)
                if isfile(join(basic_path, each_file)) and not each_file.startswith('.')]

    def start(self):
        result = []
        # windows path
        # folders = sorted(self.get_folders(), key=lambda x: int(x.split("\")[-1].split(".")[0]))
        folders = sorted(self.get_folders(), key=lambda x: int(x.split("/")[-1].split(".")[0]))
        for folder_number in range(len(folders)):
            print("Processing folder %s" % folders[folder_number])
            index = [-1.0]
            # windows path
            # files = sorted(self.get_files(folders[folder_number]), key=lambda x: int(x.split("\")[-1]))
            files = sorted(self.get_files(folders[folder_number]), key=lambda x: int(x.split("/")[-1]))
            for i in range(len(files)):
                content = ''
                contentSpare = ''
                with open(files[i], "r") as fileread:
                    oneline = fileread.readline()
                    while oneline:
                        if oneline.strip().startswith("--- full sequence ---"):
                            contentList = fileread.readlines()
                            content = contentList[2]
                            contentSpare = contentList[3]
                            break
                        oneline = fileread.readline()

                if not content.strip('\n'):
                    index.append(0.0)
                    continue

                try:
                    index.append(float(content.split()[1]))
                except ValueError:
                    index.append(float(contentSpare.split()[1]))

            # index[folder_number + 1] = -1.0
            max_number = max(index)
            print(index)
            if max_number == 0.0:
                print ("0")
                result.append(': '.join([str(folder_number + 1), "0"]) + '\n')
                continue

            for i in range(len(index)):
                if index[i] == max_number:
                    result.append(': '.join([str(folder_number + 1), str(i)]) + '\n')
        # windows path
        # with open(join(str(Path().absolute())+"\output\", "score.txt"), 'w') as writefile:
        #     writefile.writelines(result)
        with open(join(str(Path().absolute())+"/output/", "score.txt"), 'w') as writefile:
            writefile.writelines(result)

    def compute(self):
        if not os.path.exists(self.hmm_out):
            os.mkdir(self.hmm_out)
        self.hmm_build()
        self.hmm_press()
        self.hmm_scan()
        self.start()

