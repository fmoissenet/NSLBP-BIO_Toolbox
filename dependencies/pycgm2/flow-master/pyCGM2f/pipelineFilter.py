# -*- coding: utf-8 -*-
import os
import collections
from  collections import OrderedDict
import logging
import copy



from pyCGM2.Utils import files



class Pipeline():

    def __init__(self,document=None):

        # if "AQM" not in document.keys():
        #     document.update({"AQM":OrderedDict()})
        #     document["AQM"].update({"Conditions":document["Conditions"]})
        #     document.pop("Conditions")

        self.Document=document


    def __getConditionItem(self,conditionID):

        flag = False
        index = 0
        for conditionIt in self.Document["Conditions"]:
            if conditionIt["ConditionID"] == conditionID:
                flag = True
                break
            index+=1

        if flag:
            return index
        else:
            raise Exception("[pyCGM2]  Condition Id not  find")


    def _makeDocument(self):

        document = OrderedDict()

        # document.update(self.Subject)
        # document["Visits"] = list()
        # document["Visits"].append(self.Visit)
        # document["Visits"][0]["Exams"]=OrderedDict()
        # document["Visits"][0]["Exams"].update({"AQM":list()})
        # document["Visits"][0]["Exams"].append(self.Session)
        # document["Visits"][0]["Exams"]["AQM"][0].update({"Conditions":self.Conditions})

        document.update(self.Subject)
        document.update(self.Visit)
        document.update(self.Session)
        document.update({"Conditions":OrderedDict()})
        document["Conditions"] = self.Conditions


        # document["AQM"].update(self.Session)
        # document["AQM"].update({"Conditions":self.Conditions})

        return document


    def getFittingFiles(self,conditionID):
        """Short summary.

        :param conditionItem: Description of parameter `conditionItem`.
        :type conditionItem: type
        :param index: Description of parameter `index`.
        :type index: type
        :return: Description of returned object.
        :rtype: type

        """
        index = self.__getConditionItem(conditionID)

        conditionItem = self.Document["Conditions"][index]

        out = list()
        for fittingIt in conditionItem["FittingTrials"]:
            out.append(fittingIt["File"])

        return out

    def getEmgFiles(self,conditionID):

        index = self.__getConditionItem(conditionID)
        conditionItem = self.Document["Conditions"][index]

        out = list()
        for emgTrialIt in conditionItem["EmgTrials"]:
            out.append(emgTrialIt)

        return out

    def getEmgConfiguration(self,conditionID,outputType = "list"):

        index = self.__getConditionItem(conditionID)
        conditionItem = self.Document["Conditions"][index]

        labels = []
        contexts =[]
        normalActivities = []
        muscles =[]
        for emg in  conditionItem["EmgSettings"]["CHANNELS"].keys():
            if emg !="None":
                if conditionItem["EmgSettings"]["CHANNELS"][emg]["Muscle"] != "None":
                    labels.append((emg))
                    muscles.append((conditionItem["EmgSettings"]["CHANNELS"][emg]["Muscle"]))
                    contexts.append((conditionItem["EmgSettings"]["CHANNELS"][emg]["Context"])) if conditionItem["EmgSettings"]["CHANNELS"][emg] != "None" else contexts.append("NA")
                    normalActivities.append((conditionItem["EmgSettings"]["CHANNELS"][emg]["NormalActivity"])) if conditionItem["EmgSettings"]["CHANNELS"][emg]["NormalActivity"] != "None" else normalActivities.append("NA")



        if outputType == "list":
            return labels,muscles,contexts,normalActivities
        else:
            return {"Labels": labels, "Muscles": muscles, "Contexts": contexts, "NormalActivity": normalActivities}


    def getConditionFirstLevelInfos(self,conditionID):

        index = self.__getConditionItem(conditionID)
        conditionItem = self.Document["Conditions"][index]

        out = OrderedDict()
        for key in conditionItem.keys():
            if type(conditionItem[key]) not in [list, collections.OrderedDict, dict]:
                out[key] = conditionItem[key]

        return out

    def getSubjectFirstLevelInfos(self):
        out = OrderedDict()

        for key in self.Document.keys():
            if type(self.Document[key]) not in [list, collections.OrderedDict, dict]:
                out[key] = self.Document[key]
        return out


    def getNestedCalibrations(self):

        calibrations = []
        ids =[]
        count=0
        for conditionIt in self.Document["Conditions"]:
            for fittingIt in conditionIt["FittingTrials"]:
                if count==0:
                    ids.append(fittingIt["Calibration"]["ID"])
                    calibrations.append(fittingIt["Calibration"])
                else:
                    if fittingIt["Calibration"]["ID"] not in ids:
                        calibrations.append(fittingIt["Calibration"])
            count +=1
        return calibrations

    def getCalibrationTranslators(self,calibrationID):

        for conditionIt in self.Document["Conditions"]:
            for fittingIt in conditionIt["FittingTrials"]:
                if fittingIt["Calibration"]["ID"] == calibrationID:
                    return conditionIt["Translators"]

    def getCalibrationGlobal(self,calibrationID):

        for conditionIt in self.Document["Conditions"]:
            for fittingIt in conditionIt["FittingTrials"]:
                if fittingIt["Calibration"]["ID"] == calibrationID:
                    return conditionIt["Global"]

    def getHJC(self,calibrationID):

        for conditionIt in self.Document["Conditions"]:
            for fittingIt in conditionIt["FittingTrials"]:
                if fittingIt["Calibration"]["ID"] == calibrationID:
                    return fittingIt["Calibration"]["HJC"]

    def getStaticWeightSet(self,calibrationID):

        for conditionIt in self.Document["Conditions"]:
            for fittingIt in conditionIt["FittingTrials"]:
                if fittingIt["Calibration"]["ID"] == calibrationID:
                    return fittingIt["Calibration"]["Weights"]


    def getRepresentativeEMGTrial(self,conditionID):

        index = self.__getConditionItem(conditionID)
        conditionItem = self.Document["Conditions"][index]

        return conditionItem["EmgRepresentativeTrial"]



class PipelineFilter(object):
    """
         Filter building an Analysis instance.
    """

    def __init__(self):
        self.__concretePipelineBuilder = None
        self.pipeline = Pipeline()


    def setBuilder(self,concreteBuilder):
        """
        """
        self.__concretePipelineBuilder = concreteBuilder



    def build (self) :
        """
            build member analysis from a concrete builder

        """
        self.__concretePipelineBuilder._internalSettingsManipulation()
        self.__concretePipelineBuilder._check()

        setattr(self.pipeline, 'Subject', self.__concretePipelineBuilder.getSubjectInfo())
        setattr(self.pipeline, 'Visit', self.__concretePipelineBuilder.getVisitInfo())
        setattr(self.pipeline, 'Session', self.__concretePipelineBuilder.getExamInfo())

        setattr(self.pipeline, 'Conditions', self.__concretePipelineBuilder.getConditions())


        # setattr(self.pipeline, 'Calibration', self.__concretePipelineBuilder.getCalibration())
        # setattr(self.pipeline, 'Translators', self.__concretePipelineBuilder.getTranslators())

        setattr(self.pipeline, 'Document', self.pipeline._makeDocument())

    def exportAQMinfos(self, DATA_PATH,modelVersion, name="AQM-exam.info",jsonAnalyses=None):

        files.saveJson(DATA_PATH, "Subject.info", self.pipeline.Subject)

        document = OrderedDict()
        document["Ipp"] = self.pipeline.Subject["Ipp"]
        document["Model"] = modelVersion
        document.update(self.pipeline.Visit)
        document.update(self.pipeline.Session)
        document.update({"Conditions":self.pipeline.Conditions})


        if jsonAnalyses is not None:
            for index in range(0,len(document["Conditions"])):
                conditionID = document["Conditions"][index]["ConditionID"]
                if conditionID in jsonAnalyses.keys():
                    document["Conditions"][index]["DATA"] = jsonAnalyses[conditionID]

        files.saveJson(DATA_PATH, name, document)
