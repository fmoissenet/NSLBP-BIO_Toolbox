# -*- coding: utf-8 -*-
import os
import collections
from  collections import OrderedDict
import logging
import copy

import pyCGM2
from pyCGM2.Utils import files
from pyCGM2 import enums
from pyCGM2.Eclipse import vskTools,eclipse


# --- BUILDERS-----
class AbstractPipelineBuilder(object):
    def __init__(self,userSettings):

        self._localInternalSettings = None
        self._localTranslators = None

        self._userSettings = userSettings
        self._vsk = None

        if os.path.isfile(pyCGM2.PYCGM2_APPDATA_PATH + "emg.settings"):
            self._emgSettings = files.openFile(pyCGM2.PYCGM2_APPDATA_PATH,"emg.settings")
        else:
            self._emgSettings = files.openFile(pyCGM2.PYCGM2_SETTINGS_FOLDER,"emg.settings")

    def setVsk(self,vsk):
        self._vsk = vsk

    def setLocalEmgSettings(self,localEmgSettings):
        self._emgSettings = localEmgSettings


    def getSubjectInfo(self):
        return self._userSettings["SubjectInfo"]

    def getVisitInfo(self):
        return self._userSettings["VisitInfo"]

    def getExamInfo(self):
        out = self._userSettings["ExamInfo"]


        return out

    def setLocalInternalSettings(self,localInternalSettings):
        self._localInternalSettings = localInternalSettings


    def setLocalTranslators(self,translators):
        self._localTranslators = translators


    def _MP(self):

        if self._vsk is not None:
            required_mp,optional_mp = vskTools.getFromVskSubjectMp(self._vsk, resetFlag=True)
            self._userSettings["MP"]["Required"].update(required_mp)
            self._userSettings["MP"]["Optional"].update(optional_mp)

    def _check(self):


        # duplication of fitting trial
        fittingFiles=list()
        if self._userSettings["Fitting"]["Trials"] is not None:
            for fittingTrialIt in self._userSettings["Fitting"]["Trials"]:
                if fittingTrialIt["File"]  in fittingFiles:
                    raise Exception("[pipeCGM2] Check your user settings. Fitting trial name duplicated !")
                else:
                    fittingFiles.append(fittingTrialIt["File"])

        # duplication of emg trial
        EmgFiles=list()
        if self._userSettings["Emg"]["Trials"] is not None:
            for emgTrialIt in self._userSettings["Emg"]["Trials"]:
                if emgTrialIt["File"]  in EmgFiles:
                    raise Exception("[pipeCGM2] Check your user settings. emg trial name duplicated !")
                else:
                    EmgFiles.append(emgTrialIt["File"])

        # duplication of calibration ID
        calibrationsIds=list()
        if self._userSettings["Calibration"] is not None:
            for it in self._userSettings["Calibration"]:
                if it["ID"]  in calibrationsIds:
                    raise Exception("[pipeCGM2] Check your user settings. calibration ID not unique !")
                else:
                    calibrationsIds.append(it["ID"])

        # duplication of condition ID
        conditionsIds=list()
        if self._userSettings["Protocol"]["Conditions"] is not None:
            for it in self._userSettings["Protocol"]["Conditions"]:
                if it["ConditionID"]  in conditionsIds:
                    raise Exception("[pipeCGM2] Check your user settings. Condition ID not unique !")
                else:
                    conditionsIds.append(it["ConditionID"])

        # check if EMG Reference condition point a valid condition

        for it in self._userSettings["Protocol"]["Conditions"]:
            if it["EmgReferenceConditionID"] is not None:
                if it["EmgReferenceConditionID"] not in conditionsIds:
                    raise Exception("[pipeCGM2] Check your user settings. Emg reference conditionID not pointing a valid ConditionID !")






class CGM1PipelineBuilder(AbstractPipelineBuilder):
    """
        **Description** :
    """

    def __init__(self,userSettings):

        super(CGM1PipelineBuilder, self).__init__(userSettings)



    def _internalSettingsManipulation(self):

        if self._localInternalSettings is None:
            if os.path.isfile(pyCGM2.PYCGM2_APPDATA_PATH + "CGM1-pyCGM2.settings"):
                self._internSettings = files.openFile(pyCGM2.PYCGM2_APPDATA_PATH,"CGM1-pyCGM2.settings")
            else:
                self._internSettings = files.openFile(pyCGM2.PYCGM2_SETTINGS_FOLDER,"CGM1-pyCGM2.settings")
        else:
            logging.info("Local internal setting found")
            self._internSettings = self._localInternalSettings

        if self._localTranslators is None:
            self._localTranslators = self._internSettings["Translators"]
        else:
            logging.info("Local translators found")


        for key in self._userSettings["Global"].keys():
            if self._userSettings["Global"][key] is None  and self._internSettings["Global"].has_key(key):
                self._userSettings["Global"][key] = self._internSettings["Global"][key] if self._internSettings["Global"][key] !="None" else None

        self._momentProjection = self._internSettings["Fitting"]["Moment Projection"]


    def _Calibration(self):
        self._MP()

        finalCalibration = list()

        calibrationIds_fromFitting=list()

        for fittingIt in self._userSettings["Fitting"]["Trials"]:
            id =  fittingIt["CalibrationID"]
            if id not in calibrationIds_fromFitting:
                calibrationIds_fromFitting.append(id)


        for i in range(0,len(self._userSettings["Calibration"])):

            if self._userSettings["Calibration"][i]["ID"]  in calibrationIds_fromFitting:

                self._userSettings["Calibration"][i].update( self._userSettings["Global"])

                copiedTranslators = copy.deepcopy(self._localTranslators)
                if "Translators" in self._userSettings["Calibration"][i].keys():
                    if self._userSettings["Calibration"][i]["Translators"] != "Inherit":
                        for translatorKey in self._userSettings["Calibration"][i]["Translators"].keys():
                            if translatorKey in copiedTranslators.keys():
                                copiedTranslators[translatorKey] = self._userSettings["Calibration"][i]["Translators"][translatorKey]
                self._userSettings["Calibration"][i]["Translators"] = copiedTranslators


                copiedMp = copy.deepcopy(self._userSettings["MP"])
                if "MP" in self._userSettings["Calibration"][i].keys():
                    if self._userSettings["Calibration"][i]["MP"] !="Inherit":
                        for mpKey in self._userSettings["Calibration"][i]["MP"].keys():
                            valueModified = self._userSettings["Calibration"][i]["MP"][mpKey]
                            if mpKey in copiedMp["Required"].keys():
                                copiedMp["Required"][mpKey] = valueModified
                            elif mpKey in copiedMp["Optional"].keys():
                                copiedMp["Optional"][mpKey] = valueModified
                            else:
                                logging.warning ("[pyCGM2] mp key not detected")
                self._userSettings["Calibration"][i]["MP"] =copiedMp

                finalCalibration.append(self._userSettings["Calibration"][i])


        self._userSettings["Calibration"] = finalCalibration

        return self._userSettings["Calibration"]

    def getConditions(self):

        self._Calibration()

        for index in range(0,len(self._userSettings["Protocol"]["Conditions"])):

            logging.info("Traitement de la condition [%i]"%(index))

            # protocol level
            self._userSettings["Protocol"]["Conditions"][index] .update({"ResearchProtocol": self._userSettings["Protocol"]["ResearchProtocol"]})

            # Global
            self._userSettings["Protocol"]["Conditions"][index] .update({"Global": self._userSettings["Global"]})


            self._userSettings["Protocol"]["Conditions"][index] .update({"FittingTrials":list()})
            self._userSettings["Protocol"]["Conditions"][index] .update({"EmgTrials":list()})

            # EMG trials
            copiedEMGsettings = copy.deepcopy(self._emgSettings)
            if "EmgSettings"  in self._userSettings["Protocol"]["Conditions"][index].keys():
                if self._userSettings["Protocol"]["Conditions"][index] ["EmgSettings"] != "Inherit":
                    if self._userSettings["Protocol"]["Conditions"][index] ["EmgSettings"].has_key("CHANNELS"):
                        copiedEMGsettings["CHANNELS"].update(self._userSettings["Protocol"]["Conditions"][index] ["EmgSettings"]["CHANNELS"])
                    if self._userSettings["Protocol"]["Conditions"][index] ["EmgSettings"].has_key("Processing"):
                        copiedEMGsettings["Processing"].update(self._userSettings["Protocol"]["Conditions"][index] ["EmgSettings"]["Processing"])
            self._userSettings["Protocol"]["Conditions"][index] ["EmgSettings"] = copiedEMGsettings

            # ----Moment projection----
            momentProjection = self._momentProjection
            if "MomentProjection" in self._userSettings["Protocol"]["Conditions"][index].keys():
                 if self._userSettings["Protocol"]["Conditions"][index] ["MomentProjection"] != "Inherit":
                    momentProjection = self._userSettings["Protocol"]["Conditions"][index] ["MomentProjection"]
                    if not enums.isValueBelongToEnum (momentProjection,enums.MomentProjection):
                        raise Exception("[pipeCGM2] Moment projection not known" )
            self._userSettings["Protocol"]["Conditions"][index] ["MomentProjection"] = momentProjection


            #---Fitting trials---
            for fittingIt in self._userSettings["Fitting"]["Trials"]:

                #global
                fittingIt .update( self._userSettings["Global"])


                copiedTranslators = copy.deepcopy(self._localTranslators)
                if "Translators" not in fittingIt.keys() or fittingIt["Translators"] == "Inherit":
                    fittingIt["Translators"] = copiedTranslators
                else:
                    for translatorKey in fittingIt["Translators"].keys():
                        if translatorKey in copiedTranslators.keys():
                            copiedTranslators[translatorKey] = fittingIt["Translators"][translatorKey]
                    fittingIt["Translators"] = copiedTranslators


                if fittingIt["ConditionId"] == self._userSettings["Protocol"]["Conditions"][index]["ConditionID"]:
                    # ajout de la calibration associé au fitting
                    calibrationID = fittingIt["CalibrationID"]
                    for calibrationIt in self._userSettings["Calibration"]:
                        if calibrationIt["ID"] == calibrationID:
                            fittingIt["Calibration"] = calibrationIt
                            pass
                    # ajout de fitting au protocole
                    self._userSettings["Protocol"]["Conditions"][index] ["FittingTrials"].append(fittingIt)


            #--- emg trials---
            emgFiles = list()
            for fittingIt in self._userSettings["Fitting"]["Trials"]:
                if fittingIt["ConditionId"] == self._userSettings["Protocol"]["Conditions"][index]["ConditionID"]:
                    if fittingIt["Emg"]:
                        if fittingIt["File"] not in emgFiles:
                            emgFiles.append(fittingIt["File"])


            if self._userSettings["Emg"]["Trials"] is not None:
                for emgIt in self._userSettings["Emg"]["Trials"]:
                    if emgIt["ConditionId"] == self._userSettings["Protocol"]["Conditions"][index]["ConditionID"]:
                        if emgIt["File"] not in emgFiles:
                            emgFiles.append(emgIt["File"])

            self._userSettings["Protocol"]["Conditions"][index] ["EmgTrials"] = emgFiles


        return self._userSettings["Protocol"]["Conditions"]



class CGM11PipelineBuilder(CGM1PipelineBuilder):
    """
        **Description** :
    """

    def __init__(self,userSettings):

        super(CGM11PipelineBuilder, self).__init__(userSettings)



    def _internalSettingsManipulation(self):

        if self._localInternalSettings is None:
            if os.path.isfile(pyCGM2.PYCGM2_APPDATA_PATH + "CGM1_1-pyCGM2.settings"):
                self._internSettings = files.openFile(pyCGM2.PYCGM2_APPDATA_PATH,"CGM1_1-pyCGM2.settings")
            else:
                self._internSettings = files.openFile(pyCGM2.PYCGM2_SETTINGS_FOLDER,"CGM1_1-pyCGM2.settings")
        else:
            logging.info("Local internal setting found")
            self._internSettings = self._localInternalSettings

        if self._localTranslators is None:
            self._localTranslators = self._internSettings["Translators"]
        else:
            logging.info("Local translators found")


        for key in self._userSettings["Global"].keys():
            if self._userSettings["Global"][key] is None  and self._internSettings["Global"].has_key(key):
                self._userSettings["Global"][key] = self._internSettings["Global"][key] if self._internSettings["Global"][key] !="None" else None

        self._momentProjection = self._internSettings["Fitting"]["Moment Projection"]



class CGM21PipelineBuilder(CGM1PipelineBuilder):
    """
        **Description** :
    """

    def __init__(self,userSettings):

        super(CGM21PipelineBuilder, self).__init__(userSettings)


    def _internalSettingsManipulation(self):

        if self._localInternalSettings is None:
            if os.path.isfile(pyCGM2.PYCGM2_APPDATA_PATH + "CGM2_1-pyCGM2.settings"):
                self._internSettings = files.openFile(pyCGM2.PYCGM2_APPDATA_PATH,"CGM2_1-pyCGM2.settings")
            else:
                self._internSettings = files.openFile(pyCGM2.PYCGM2_SETTINGS_FOLDER,"CGM2_1-pyCGM2.settings")
        else:
            logging.info("Local internal setting found")
            self._internSettings = self._localInternalSettings

        if self._localTranslators is None:
            self._localTranslators = self._internSettings["Translators"]
        else:
            logging.info("Local translators found")

        for key in self._userSettings["Global"].keys():
            if self._userSettings["Global"][key] is None  and self._internSettings["Global"].has_key(key):
                self._userSettings["Global"][key] = self._internSettings["Global"][key] if self._internSettings["Global"][key] !="None" else None

        self._momentProjection = self._internSettings["Fitting"]["Moment Projection"]

        self._HJC = self._internSettings["Calibration"]["HJC"]

    def _Calibration(self):
        self._MP()
        finalCalibration = list()

        calibrationIds_fromFitting=list()

        for fittingIt in self._userSettings["Fitting"]["Trials"]:
            id =  fittingIt["CalibrationID"]
            if id not in calibrationIds_fromFitting:
                calibrationIds_fromFitting.append(id)


        for i in range(0,len(self._userSettings["Calibration"])):

            if self._userSettings["Calibration"][i]["ID"]  in calibrationIds_fromFitting:

                self._userSettings["Calibration"][i].update( self._userSettings["Global"])

                copiedHJC = copy.deepcopy(self._HJC)
                if "HJC" in self._userSettings["Calibration"][i].keys():
                    if self._userSettings["Calibration"][i]["HJC"] != "Inherit":
                        for hjcKey in self._userSettings["Calibration"][i]["HJC"].keys():
                            valueModified = self._userSettings["Calibration"][i]["HJC"][hjcKey]
                            if hjcKey in copiedHJC.keys():
                                copiedHJC[hjcKey] = valueModified
                self._userSettings["Calibration"][i]["HJC"] =copiedHJC



                copiedTranslators = copy.deepcopy(self._localTranslators)
                if "Translators" in self._userSettings["Calibration"][i].keys():
                    if self._userSettings["Calibration"][i]["Translators"] != "Inherit":
                        for translatorKey in self._userSettings["Calibration"][i]["Translators"].keys():
                            if translatorKey in copiedTranslators.keys():
                                copiedTranslators[translatorKey] = self._userSettings["Calibration"][i]["Translators"][translatorKey]
                self._userSettings["Calibration"][i]["Translators"] = copiedTranslators


                copiedMp = copy.deepcopy(self._userSettings["MP"])
                if "MP" in self._userSettings["Calibration"][i].keys():
                    if self._userSettings["Calibration"][i]["MP"] !="Inherit":
                        for mpKey in self._userSettings["Calibration"][i]["MP"].keys():
                            valueModified = self._userSettings["Calibration"][i]["MP"][mpKey]
                            if mpKey in copiedMp["Required"].keys():
                                copiedMp["Required"][mpKey] = valueModified
                            elif mpKey in copiedMp["Optional"].keys():
                                copiedMp["Optional"][mpKey] = valueModified
                            else:
                                logging.warning ("[pyCGM2] mp key not detected")
                self._userSettings["Calibration"][i]["MP"] =copiedMp


                finalCalibration.append(self._userSettings["Calibration"][i])





        self._userSettings["Calibration"] = finalCalibration

        return self._userSettings["Calibration"]



class CGM22PipelineBuilder(CGM21PipelineBuilder):
    """
        **Description** :
    """

    def __init__(self,userSettings):

        super(CGM22PipelineBuilder, self).__init__(userSettings)


    def _internalSettingsManipulation(self):

        if self._localInternalSettings is None:
            if os.path.isfile(pyCGM2.PYCGM2_APPDATA_PATH + "CGM2_2-pyCGM2.settings"):
                self._internSettings = files.openFile(pyCGM2.PYCGM2_APPDATA_PATH,"CGM2_2-pyCGM2.settings")
            else:
                self._internSettings = files.openFile(pyCGM2.PYCGM2_SETTINGS_FOLDER,"CGM2_2-pyCGM2.settings")
        else:
            logging.info("Local internal setting found")
            self._internSettings = self._localInternalSettings

        if self._localTranslators is None:
            self._localTranslators = self._internSettings["Translators"]
        else:
            logging.info("Local translators found")


        for key in self._userSettings["Global"].keys():
            if self._userSettings["Global"][key] is None  and self._internSettings["Global"].has_key(key):
                self._userSettings["Global"][key] = self._internSettings["Global"][key] if self._internSettings["Global"][key] !="None" else None

        self._momentProjection = self._internSettings["Fitting"]["Moment Projection"]

        self._HJC = self._internSettings["Calibration"]["HJC"]
        self._weights = self._internSettings["Fitting"]["Weight"]


    def _Calibration(self):

        self._MP()
        finalCalibration = list()

        calibrationIds_fromFitting=list()

        for fittingIt in self._userSettings["Fitting"]["Trials"]:
            id =  fittingIt["CalibrationID"]
            if id not in calibrationIds_fromFitting:
                calibrationIds_fromFitting.append(id)


        for i in range(0,len(self._userSettings["Calibration"])):

            if self._userSettings["Calibration"][i]["ID"]  in calibrationIds_fromFitting:

                self._userSettings["Calibration"][i].update( self._userSettings["Global"])


                copiedHJC = copy.deepcopy(self._HJC)
                if "HJC" in self._userSettings["Calibration"][i].keys():
                    if self._userSettings["Calibration"][i]["HJC"] != "Inherit":
                        for hjcKey in self._userSettings["Calibration"][i]["HJC"].keys():
                            valueModified = self._userSettings["Calibration"][i]["HJC"][hjcKey]
                            if hjcKey in copiedHJC.keys():
                                copiedHJC[hjcKey] = valueModified

                self._userSettings["Calibration"][i]["HJC"] =copiedHJC


                copiedWeights = copy.deepcopy(self._weights)
                if "Weights" in self._userSettings["Calibration"][i].keys():
                    if self._userSettings["Calibration"][i]["Weights"] != "Inherit":
                        for weightKey in self._userSettings["Calibration"][i]["Weights"].keys():
                            if weightKey in copiedWeights:
                                copiedWeights[weightKey] = self._userSettings["Calibration"][i]["Weights"][weightKey]
                self._userSettings["Calibration"][i]["Weights"] = copiedWeights



                copiedTranslators = copy.deepcopy(self._localTranslators)
                if "Translators" in self._userSettings["Calibration"][i].keys():
                    if self._userSettings["Calibration"][i]["Translators"] != "Inherit":
                        for translatorKey in self._userSettings["Calibration"][i]["Translators"].keys():
                            if translatorKey in copiedTranslators.keys():
                                copiedTranslators[translatorKey] = self._userSettings["Calibration"][i]["Translators"][translatorKey]

                self._userSettings["Calibration"][i]["Translators"] = copiedTranslators


                copiedMp = copy.deepcopy(self._userSettings["MP"])
                if "MP" in self._userSettings["Calibration"][i].keys():
                    if self._userSettings["Calibration"][i]["MP"] !="Inherit":
                        for mpKey in self._userSettings["Calibration"][i]["MP"].keys():
                            valueModified = self._userSettings["Calibration"][i]["MP"][mpKey]
                            if mpKey in copiedMp["Required"].keys():
                                copiedMp["Required"][mpKey] = valueModified
                            elif mpKey in copiedMp["Optional"].keys():
                                copiedMp["Optional"][mpKey] = valueModified
                            else:
                                logging.warning ("[pyCGM2] mp key not detected")
                self._userSettings["Calibration"][i]["MP"] =copiedMp


                finalCalibration.append(self._userSettings["Calibration"][i])


        self._userSettings["Calibration"] = finalCalibration

        return self._userSettings["Calibration"]



    def getConditions(self):

        self._Calibration()

        for index in range(0,len(self._userSettings["Protocol"]["Conditions"])):

            logging.info("Traitement de la condition [%i]"%(index))

            # protocol level
            self._userSettings["Protocol"]["Conditions"][index] .update({"ResearchProtocol": self._userSettings["Protocol"]["ResearchProtocol"]})

            # Global
            self._userSettings["Protocol"]["Conditions"][index] .update({"Global": self._userSettings["Global"]})

            self._userSettings["Protocol"]["Conditions"][index] .update({"FittingTrials":list()})
            self._userSettings["Protocol"]["Conditions"][index] .update({"EmgTrials":list()})

            # EMG trials
            copiedEMGsettings = copy.deepcopy(self._emgSettings)
            if "EmgSettings"  in self._userSettings["Protocol"]["Conditions"][index].keys():
                if self._userSettings["Protocol"]["Conditions"][index] ["EmgSettings"] != "Inherit":
                    if self._userSettings["Protocol"]["Conditions"][index] ["EmgSettings"].has_key("CHANNELS"):
                        copiedEMGsettings["CHANNELS"].update(self._userSettings["Protocol"]["Conditions"][index] ["EmgSettings"]["CHANNELS"])
                    if self._userSettings["Protocol"]["Conditions"][index] ["EmgSettings"].has_key("Processing"):
                        copiedEMGsettings["Processing"].update(self._userSettings["Protocol"]["Conditions"][index] ["EmgSettings"]["Processing"])
            self._userSettings["Protocol"]["Conditions"][index] ["EmgSettings"] = copiedEMGsettings




            # ----Moment projection----

            momentProjection = self._momentProjection

            if "MomentProjection" in self._userSettings["Protocol"]["Conditions"][index].keys():
                 if self._userSettings["Protocol"]["Conditions"][index] ["MomentProjection"] != "Inherit":
                    momentProjection = self._userSettings["Protocol"]["Conditions"][index] ["MomentProjection"]
                    if not enums.isValueBelongToEnum (momentProjection,enums.MomentProjection):
                        raise Exception("[pipeCGM2] Moment projection not known" )
            self._userSettings["Protocol"]["Conditions"][index] ["MomentProjection"] = momentProjection



            #---Fitting trials---
            for fittingIt in self._userSettings["Fitting"]["Trials"]:

                #global
                fittingIt .update( self._userSettings["Global"])



                copiedTranslators = copy.deepcopy(self._localTranslators)
                if "Translators" in fittingIt.keys():
                    if fittingIt["Translators"] != "Inherit":
                        for translatorKey in fittingIt["Translators"].keys():
                            if translatorKey in copiedTranslators.keys():
                                copiedTranslators[translatorKey] = fittingIt["Translators"][translatorKey]
                fittingIt["Translators"] = copiedTranslators



                copiedWeights = copy.deepcopy(self._weights)
                if "Weights" in fittingIt.keys():
                    if fittingIt["Weights"] != "Inherit" :
                        for weightKey in fittingIt["Weights"].keys():
                            if weightKey in copiedWeights.keys():
                                copiedWeights[weightKey] = fittingIt["Weights"][weightKey]
                fittingIt["Weights"] = copiedWeights



                if fittingIt["ConditionId"] == self._userSettings["Protocol"]["Conditions"][index]["ConditionID"]:
                    # ajout de la calibration associé au fitting
                    calibrationID = fittingIt["CalibrationID"]
                    for calibrationIt in self._userSettings["Calibration"]:
                        if calibrationIt["ID"] == calibrationID:
                            fittingIt["Calibration"] = calibrationIt
                            pass
                    # ajout de fitting au protocole
                    self._userSettings["Protocol"]["Conditions"][index] ["FittingTrials"].append(fittingIt)


            #--- emg trials---
            emgFiles = list()
            for fittingIt in self._userSettings["Fitting"]["Trials"]:
                if fittingIt["ConditionId"] == self._userSettings["Protocol"]["Conditions"][index]["ConditionID"]:
                    if fittingIt["Emg"]:
                        if fittingIt["File"] not in emgFiles:
                            emgFiles.append(fittingIt["File"])

            if self._userSettings["Emg"]["Trials"] is not None:
                for emgIt in self._userSettings["Emg"]["Trials"]:
                    if emgIt["ConditionId"] == self._userSettings["Protocol"]["Conditions"][index]["ConditionID"]:
                        if emgIt["File"] not in emgFiles:
                            emgFiles.append(emgIt["File"])

            self._userSettings["Protocol"]["Conditions"][index] ["EmgTrials"] = emgFiles


        return self._userSettings["Protocol"]["Conditions"]


class CGM23PipelineBuilder(CGM22PipelineBuilder):
    """
        **Description** :
    """

    def __init__(self,userSettings):

        super(CGM23PipelineBuilder, self).__init__(userSettings)

    def _internalSettingsManipulation(self):

        if self._localInternalSettings is None:
            if os.path.isfile(pyCGM2.PYCGM2_APPDATA_PATH + "CGM2_3-pyCGM2.settings"):
                self._internSettings = files.openFile(pyCGM2.PYCGM2_APPDATA_PATH,"CGM2_3-pyCGM2.settings")
            else:
                self._internSettings = files.openFile(pyCGM2.PYCGM2_SETTINGS_FOLDER,"CGM2_3-pyCGM2.settings")
        else:
            logging.info("Local internal setting found")
            self._internSettings = self._localInternalSettings

        if self._localTranslators is None:
            self._localTranslators = self._internSettings["Translators"]
        else:
            logging.info("Local translators found")


        for key in self._userSettings["Global"].keys():
            if self._userSettings["Global"][key] is None  and self._internSettings["Global"].has_key(key):
                self._userSettings["Global"][key] = self._internSettings["Global"][key] if self._internSettings["Global"][key] !="None" else None

        self._momentProjection = self._internSettings["Fitting"]["Moment Projection"]

        self._HJC = self._internSettings["Calibration"]["HJC"]
        self._weights = self._internSettings["Fitting"]["Weight"]

class CGM24PipelineBuilder(CGM23PipelineBuilder):
    """
        **Description** :
    """

    def __init__(self,userSettings):

        super(CGM23PipelineBuilder, self).__init__(userSettings)

    def _internalSettingsManipulation(self):

        if self._localInternalSettings is None:
            if os.path.isfile(pyCGM2.PYCGM2_APPDATA_PATH + "CGM2_4-pyCGM2.settings"):
                self._internSettings = files.openFile(pyCGM2.PYCGM2_APPDATA_PATH,"CGM2_4-pyCGM2.settings")
            else:
                self._internSettings = files.openFile(pyCGM2.PYCGM2_SETTINGS_FOLDER,"CGM2_4-pyCGM2.settings")
        else:
            logging.info("Local internal setting found")
            self._internSettings = self._localInternalSettings

        if self._localTranslators is None:
            self._localTranslators = self._internSettings["Translators"]
        else:
            logging.info("Local translators found")


        for key in self._userSettings["Global"].keys():
            if self._userSettings["Global"][key] is None  and self._internSettings["Global"].has_key(key):
                self._userSettings["Global"][key] = self._internSettings["Global"][key] if self._internSettings["Global"][key] !="None" else None

        self._momentProjection = self._internSettings["Fitting"]["Moment Projection"]

        self._HJC = self._internSettings["Calibration"]["HJC"]
        self._weights = self._internSettings["Fitting"]["Weight"]


class CGM25PipelineBuilder(CGM24PipelineBuilder):
    """
        **Description** :
    """

    def __init__(self,userSettings):

        super(CGM23PipelineBuilder, self).__init__(userSettings)


    def _internalSettingsManipulation(self):

        if self._localInternalSettings is None:
            if os.path.isfile(pyCGM2.PYCGM2_APPDATA_PATH + "CGM2_4-pyCGM2.settings"):
                self._internSettings = files.openFile(pyCGM2.PYCGM2_APPDATA_PATH,"CGM2_5-pyCGM2.settings")
            else:
                self._internSettings = files.openFile(pyCGM2.PYCGM2_SETTINGS_FOLDER,"CGM2_5-pyCGM2.settings")
        else:
            logging.info("Local internal setting found")
            self._internSettings = self._localInternalSettings

        if self._localTranslators is None:
            self._localTranslators = self._internSettings["Translators"]
        else:
            logging.info("Local translators found")


        for key in self._userSettings["Global"].keys():
            if self._userSettings["Global"][key] is None  and self._internSettings["Global"].has_key(key):
                self._userSettings["Global"][key] = self._internSettings["Global"][key] if self._internSettings["Global"][key] !="None" else None

        self._momentProjection = self._internSettings["Fitting"]["Moment Projection"]

        self._HJC = self._internSettings["Calibration"]["HJC"]
        self._weights = self._internSettings["Fitting"]["Weight"]
