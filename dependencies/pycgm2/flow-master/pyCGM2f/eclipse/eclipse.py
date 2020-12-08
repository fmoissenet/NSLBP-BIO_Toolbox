# coding: utf-8
import warnings
import logging

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
import os
import argparse
from jinja2 import Template

from pyCGM2.Utils import files
from pyCGM2.Eclipse import vskTools,eclipse
from pyCGM2 import enums

import pyCGM2f
from pyCGM2f.eclipse import eclipseInterface



def main(CGMVersion):

    TEMPLATE_PATH = pyCGM2f.MAIN_APPS_PATH+"templates\\"
    if CGMVersion  not in ["CGM1","CGM11","CGM21","CGM22","CGM23","CGM24","CGM25"]:
        raise Exception("CGM version not known ( choice is CGM1 CGM11 CGM21  CGM22 CGM23 CGM24 CGM25)")

    template = TEMPLATE_PATH+CGMVersion+".userSettings"
    jinja2_template_string = open(template, 'rb').read()
    template = Template(jinja2_template_string.decode("utf-8"))

    DATA_PATH = os.getcwd()+"\\"
    parent = os.path.abspath(os.path.join(os.getcwd(),os.pardir))+"\\"

    try:
        vskFile = vskTools.getVskFiles(DATA_PATH)
        vsk = vskTools.Vsk(str(DATA_PATH +  vskFile))
        required_mp,optional_mp = vskTools.getFromVskSubjectMp(vsk, resetFlag=True)
    except:
        logging.warning("[pyCGM2f-eclipse] No vsk detected, mp section of the user settings not initiated")
        required_mp = None




    # DATA_PATH = "C:\\Users\\FLEBOEUF.CHU-NANTES\\Documents\\Programmation\\pyCGM2\\pyCGM2-extensions\\pyCGM2_flow\\data\\eclipse\\Session 1\\"
    # parent = "C:\\Users\\FLEBOEUF.CHU-NANTES\\Documents\\Programmation\\pyCGM2\\pyCGM2-extensions\\pyCGM2_flow\\data\\eclipse\\"

    eclipseInterface.repareEnf(DATA_PATH)

    patientEnfFile =  eclipse.getEnfFiles(parent,enums.EclipseType.Patient)
    patientInfo = eclipse.PatientEnfReader(parent,patientEnfFile)

    patient = dict()
    patient["PatientID"] = patientInfo.get("PatientID")


    sessionEnfFile =  eclipse.getEnfFiles(DATA_PATH,enums.EclipseType.Session)
    sessionInfo = eclipse.SessionEnfReader(DATA_PATH,sessionEnfFile)

    visit = dict()
    visit["Age"] = sessionInfo.get("Age")
    visit["SessionID"] = sessionInfo.get("SessionID")

    eclipseDate = sessionInfo.get("CREATIONDATEANDTIME").split(",")
    visit["Date"] = str(eclipseDate[2]) +"-"+ str(eclipseDate[1]) +"-"+ str(eclipseDate[0])

    staticTrials = eclipseInterface.findStaticTrials(DATA_PATH)

    calibs = eclipseInterface.staticDetails(DATA_PATH, staticTrials)

    motionTrials = eclipseInterface.findMotionTrials(DATA_PATH)
    emgTrials = eclipseInterface.findEmgTrials(DATA_PATH)

    fits = eclipseInterface.FittingDetails(DATA_PATH,motionTrials)
    emgs = eclipseInterface.EmgDetails(DATA_PATH,emgTrials)
    conditions = eclipseInterface.getConditions(DATA_PATH,motionTrials+emgTrials)



    data = {"Patient":patient,
            "Visit":visit,
            "Mp": required_mp,
            "Calibration":calibs,
            "Fitting":fits,
            "Emg":emgs,
            "Conditions":conditions
            }


    template.stream(data=data).dump(DATA_PATH + CGMVersion+".userSettings")
