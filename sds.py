import os, subprocess,sys
import re
from pathlib import Path
from os.path import isfile, join


class SDS:

    def __init__(self,name,data):
        self.name = name
        self.data = data

    def compute(self):
        if self.name == "tcs":
            output = self.compute_tcs()
        elif self.name == "xdet":
            output = self.compute_xdet()
        return output

    def compute_tcs(self):
        # windows path
        # tm_coffee_out = str(Path().absolute())+"\TCS_result\"
        tm_coffee_out = str(Path().absolute()) + "/TCS_result/"
        if not os.path.exists(tm_coffee_out):
            os.mkdir(tm_coffee_out)
        # windows path
        # tcs_exe=str(Path().absolute())+"\tools\tcoffee-master\compile\t_coffee"
        tcs_exe = str(Path().absolute()) + "/tools/tcoffee-master/compile/t_coffee"
        d = os.walk(self.data, topdown=True)
        tm_coffee_valid_file_list = []
        for path, dir_list, file_list in d:
            for file_name in file_list:
                prog = re.compile('^\d')
                result = prog.match(file_name)
                if result:
                    tm_coffee_valid_file_list.append(file_name)
                else:
                    print('Invalid File Name :' + file_name)
            for file_name in tm_coffee_valid_file_list:
                tm_coffee_input_path = path + file_name
                tm_coffee_fasta_file_name = file_name[0:(len(file_name) - 6)]
                tm_coffee_out_put_path = tm_coffee_out + tm_coffee_fasta_file_name + '.fasta'
                subprocess.call([tcs_exe, "-infile", tm_coffee_input_path, "-outfile", tm_coffee_out_put_path,
                                 "-evaluate", "-output", "fasta"])
        return tm_coffee_out

    def compute_xdet(self):
        self.detec_pos()
        # windows path
        # sds_xdet_result=str(Path().absolute())+"\SDS_XDet_result\"
        sds_xdet_result = str(Path().absolute()) + "/SDS_XDet_result/"
        if os.path.exists(sds_xdet_result):
            os.mkdir(sds_xdet_result)
        self.startPoint()
        return sds_xdet_result

    def detec_pos(self):
        # windows path
        # metrix_path = str(Path().absolute())+"\metrix\blosum45.bla"
        # path_out = str(Path().absolute())+"\XDet_result\"
        metrix_path = str(Path().absolute())+"/metrix/blosum45.bla"
        path_out = str(Path().absolute())+"/XDet_result/"
        if not os.path.exists(path_out):
            os.mkdir(path_out)
        c = os.walk(self.data, topdown=True)
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
                filename = file_name[0:(len(file_name)-11)]
                out_put_path = path_out + filename + '.fasta'
                with open(out_put_path, 'w') as outputFile:
                    # windows path
                    # xdet_exe=str(Path().absolute())+"\tools\JDet\programs\xdet_osx"
                    xdet_exe = str(Path().absolute())+"/tools/JDet/programs/xdet_osx"
                    subprocess.call([xdet_exe, input_path, metrix_path, out_put_path, "-S", "10"], stdout=outputFile)

    def readFileList(self):
        # get original input files for the speer server
        org_files = [join(self.data, each_file) for each_file in os.listdir(self.data)
                     if isfile(join(self.data, each_file)) and not each_file.startswith('.')]

        # get Xdet output files
        # windows path
        # out_files = [join(str(Path().absolute()) + "\XDet_result\", each_file) for each_file in
        #              os.listdir(str(Path().absolute()) + "\XDet_result\")
        #              if isfile(join(str(Path().absolute()) + "\XDet_result\", each_file)) and not each_file.startswith(
        #         '.')]
        out_files = [join(str(Path().absolute())+"/XDet_result/", each_file) for each_file in os.listdir(str(Path().absolute())+"/XDet_result/")
                     if isfile(join(str(Path().absolute())+"/XDet_result/", each_file)) and not each_file.startswith('.')]

        if len(org_files) != len(out_files):
            print("The number of input files is not same as the number of output files")
            sys.exit(0)

        return org_files, out_files

    def readOriInput(self,filename):
        recordList = []
        temp = []
        key = ''
        with open(filename, 'r') as inputfile:
            oneline = inputfile.readline()
            while oneline:
                if oneline.startswith('>'):
                    if not temp:
                        key = oneline.replace('\n', '')
                    else:
                        recordList.append({key: ''.join(temp)})
                        key = oneline.replace('\n', '')
                        temp = []
                else:
                    temp.append(oneline.replace('\n', ''))
                oneline = inputfile.readline()

        recordList.append({key: ''.join(temp)})
        return recordList

    def readOutputFile(self,filename):
        positions = []
        with open(filename, 'r') as outputfile:
            oneline = outputfile.readline()
            while oneline:
                parts = oneline.replace('\n', '').split()
                positions.append(int(parts[0]))
                oneline = outputfile.readline()

        positions = set(positions)
        return positions

    def writeToFile(self,output, filename):
        # windows path
        # filename = str(Path().absolute()) + "\SDS_XDet_result\" + filename
        filename = str(Path().absolute())+"/SDS_XDet_result/" + filename
        with open(filename, 'w') as outputfile:
            for eachData in output:
                outputfile.writelines(eachData)
            outputfile.flush()

    def startPoint(self):
        path_splitor = "/"
        output = []
        dashStat = []
        original_files, output_files = self.readFileList()
        for i in range(len(original_files)):
            positions = self.readOutputFile(output_files[i])
            recordList = self.readOriInput(original_files[i])
            counter = 0
            length = 0
            counting = True
            for eachRecord in recordList:
                original = ['dummy']
                item = list(eachRecord)[0]
                key = item[0:len(item)]
                value = list(eachRecord.values())[0]
                original += list(value)
                target = []

                for j in range(len(original)):
                    if j in positions:
                        target.append(original[j])
                        if counting:
                            counter += 1
                    else:
                        target.append('-')

                output.append({key: ''.join(target)})

                if counting:
                    length = (len(original) - 1)
                    counting = False

            toFile = []
            for element in output:
                item = list(element)[0]
                toFile.append('\n'.join([item[0:len(item)], list(element.values())[0], '']))

            self.writeToFile(toFile, original_files[i].split(path_splitor)[-1])
            dashStat.append({original_files[i].split(path_splitor)[-1]: '/'.join([str(counter), str(length)])})

            print("Replace characters to '-' is done for %s" % original_files[i].split(path_splitor)[-1])
            output = []

        statistics = []
        for element in dashStat:
            item = list(element)[0]
            statistics.append(': '.join([item[0], list(element.values())[0]]) + '\n')
        # windows path
        # self.writeToFile(statistics, str(Path().absolute())+"\SDS_XDet_result\dash_statistics.txt")
        self.writeToFile(statistics, str(Path().absolute())+"/SDS_XDet_result/dash_statistics.txt")