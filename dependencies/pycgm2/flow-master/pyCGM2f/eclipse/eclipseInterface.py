# coding: utf-8
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
import os

import logging

from pyCGM2.Utils import files
from pyCGM2.Eclipse import vskTools,eclipse
from pyCGM2 import enums
import re




def findStaticTrials(path):
    enfs = eclipse.getEnfFiles(path,enums.EclipseType.Trial)

    detected = list()
    for enf in enfs:
        enfTrial = eclipse.TrialEnfReader(path,enf)
        if enfTrial.get("TrialType") is not None:
            if enfTrial.get("TrialType") == "Static":
                detected.append(enf)

    if detected ==[] : raise Exception("No static file detected")
    return detected

def findMotionTrials(path):
    enfs = eclipse.getEnfFiles(path,enums.EclipseType.Trial)

    detected = list()
    for enf in enfs:
        enfTrial = eclipse.TrialEnfReader(path,enf)
        if enfTrial.get("TrialType") is not None:
            if enfTrial.get("TrialType") == "Motion":
                detected.append(enf)

    if detected ==[] : raise Exception("No Motion file detected")
    return detected

def findEmgTrials(path):
    enfs = eclipse.getEnfFiles(path,enums.EclipseType.Trial)

    detected = list()
    for enf in enfs:
        enfTrial = eclipse.TrialEnfReader(path,enf)
        if enfTrial.get("TrialType") is not None:
            if enfTrial.get("TrialType") == "EMG":
                detected.append(enf)

    if detected ==[] : logging.info("No emg file detected")
    return detected



def staticDetails(path,enfs):

    id_checking = None
    details = list()

    for enf in enfs:

        enfTrial = eclipse.TrialEnfReader(path,enf)
        lff = enfTrial.get("LeftFlatFoot")
        rff = enfTrial.get("RightFlatFoot")
        id = enfTrial.get("CalibID")
        if id == id_checking:
            raise Exception("Calibration ID is not unique")
        else:
            id_checking = id

        file = enf[0:enf.find(".")] +".c3d"
        staticDetails =[id,file,lff,rff,]


        details.append(staticDetails)


    return details


def FittingDetails(path,enfs):

    details = list()

    for enf in enfs:

        enfTrial = eclipse.TrialEnfReader(path,enf)
        file = enf[0:enf.find(".")] +".c3d"
        conditionID = enfTrial.get("ConditionID")
        calibID = enfTrial.get("CalibID")
        fpa = enfTrial.getForcePlateAssigment()

        fittingDetails =[file,calibID,fpa,conditionID]

        details.append(fittingDetails)


    return details


def EmgDetails(path,enfs):

    details = list()

    for enf in enfs:

        enfTrial = eclipse.TrialEnfReader(path,enf)
        file = enf[0:enf.find(".")] +".c3d"
        conditionID = enfTrial.get("ConditionID")

        emgDetails =[file,conditionID]

        details.append(emgDetails)


    return details

def getConditions(path,motionEnfs):

    conditions = list()


    for enf in motionEnfs:
        enfTrial = eclipse.TrialEnfReader(path,enf)

        if enfTrial.get("Share"):
            condition=dict()

            condition["ConditionID"] = enfTrial.get("ConditionID")
            condition["Context"] = enfTrial.get("Context")
            condition["Block"] = enfTrial.get("Block")
            condition["Task"] = enfTrial.get("Task")
            condition["Shoes"] = enfTrial.get("Shoes")
            condition["ProthesisOrthosis"] = enfTrial.get("ProthesisOrthosis")
            condition["ExternalAid"] = enfTrial.get("ExternalAid")
            condition["PersonalAid"] = enfTrial.get("PersonalAid")
            condition["EmgRepresentativeTrial"] = enf[0:enf.find(".")] +".c3d"
            condition["EmgReferenceConditionID"] = ""

            conditions.append(condition)

    sortConditions = [None]* len(conditions)

    for condition in conditions:
        pos = int(re.findall("\d+",condition["ConditionID"])[0])-1
        sortConditions[pos] = condition

        if pos != 0:
            sortConditions[pos]["EmgReferenceConditionID"] = "Condition1"


    return sortConditions


def repareEnf(DATA_PATH):

    enfs = files.getFiles(DATA_PATH,"enf")
    for enf in enfs:
        sep0 = enf.find(".")
        if enf[sep0:].find("Trial")!=-1:
            if  enf[sep0:]!=".Trial.enf":
                newName = enf[0:sep0]
                os.rename(DATA_PATH+enf,DATA_PATH+newName+".Trial.enf")

        if enf[sep0:].find("Session")!=-1:
            if enf[sep0:]!=".Session.enf":
                newName = enf[0:sep0]
                os.rename(DATA_PATH+enf,DATA_PATH+newName+".Session.enf")
