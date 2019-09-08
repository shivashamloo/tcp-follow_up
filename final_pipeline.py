from division import DIVISION
from making_database import DATABASE
from blastrun import BlastRun
from msa import MSA
from sds import SDS
from hmm import HMM
from latex import ReadResult


class FinalPipeline:

    def __init__(self,train,test,database,msa,sds):
        self.train=train
        self.test=test
        self.database=database
        self.msa=msa
        self.sds=sds

    def result(self):
        training_set = DIVISION(self.train, "train")
        training_set.divide()
        testing_set = DIVISION(self.test, "test")
        testing_set.divide()
        database = DATABASE(self.test, self.database)
        database.compute()
        blastrun = BlastRun()
        msa_input = blastrun.compute()
        msa = MSA(self.msa, msa_input)
        sds_input = msa.compute()
        sds = SDS(self.sds, sds_input)
        hmm_input = sds.compute()
        hmm = HMM(hmm_input)
        hmm.compute()
        final_result = ReadResult
        final_result.read_result()



