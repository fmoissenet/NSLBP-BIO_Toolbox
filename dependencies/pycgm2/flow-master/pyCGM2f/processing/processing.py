# coding: utf-8
import logging
from  collections import OrderedDict


class ProcessingFilter(object):

    def __init__(self,procedure,manager,plotFlag, jsonAnalyses=OrderedDict()):
        self.m_procedure = procedure
        self.m_manager = manager
        self.m_plotFlag = plotFlag

        self.m_jsonAnalysisInstances=jsonAnalyses

    def getJsonAnalyses(self):
        return self.m_jsonAnalysisInstances


    def process(self,condition):
        jsonContent =  self.m_procedure.process(self.m_manager, condition, self.m_plotFlag)
        self.m_jsonAnalysisInstances[condition["ConditionID"]] = jsonContent
