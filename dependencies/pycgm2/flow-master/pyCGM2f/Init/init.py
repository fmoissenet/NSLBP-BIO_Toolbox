# coding: utf-8
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
import os
import argparse

import pyCGM2
import pyCGM2f

from pyCGM2.Utils import files




def main(expert,model):


    DATA_PATH = os.getcwd()+"\\"

    files.createDir(DATA_PATH+"Videos")
    files.createDir(DATA_PATH+"Exams")
    files.createDir(DATA_PATH+"Images")
    files.createDir(DATA_PATH+"Doc")

    files.copyPaste(pyCGM2f.MAIN_APPS_PATH+"reportGenerator.py",DATA_PATH+"reportGenerator.py")

    files.copyPaste(pyCGM2.PYCGM2_SETTINGS_FOLDER+"emg.settings", os.getcwd()+"\\"+"emg.settings")

    if expert:
        if model == "CGM1":
            files.copyPaste(pyCGM2.PYCGM2_SETTINGS_FOLDER+"CGM1-pyCGM2.settings", DATA_PATH+"CGM1-pyCGM2.settings")
        if model == "CGM11":
            files.copyPaste(pyCGM2.PYCGM2_SETTINGS_FOLDER+"CGM1_1-pyCGM2.settings", DATA_PATH+"CGM1_1-pyCGM2.settings")
        if model == "CGM21":
            files.copyPaste(pyCGM2.PYCGM2_SETTINGS_FOLDER+"CGM2_1-pyCGM2.settings", DATA_PATH+"CGM2_1-pyCGM2.settings")
        if model == "CGM22":
            files.copyPaste(pyCGM2.PYCGM2_SETTINGS_FOLDER+"CGM2_2-pyCGM2.settings", DATA_PATH+"CGM2_2-pyCGM2.settings")
        if model == "CGM23":
            files.copyPaste(pyCGM2.PYCGM2_SETTINGS_FOLDER+"CGM2_3-pyCGM2.settings", DATA_PATH+"CGM2_3-pyCGM2.settings")
        if model == "CGM24":
            files.copyPaste(pyCGM2.PYCGM2_SETTINGS_FOLDER+"CGM2_4-pyCGM2.settings", DATA_PATH+"CGM2_4-pyCGM2.settings")
        if model == "CGM25":
            files.copyPaste(pyCGM2.PYCGM2_SETTINGS_FOLDER+"CGM2_5-pyCGM2.settings", DATA_PATH+"CGM2_5-pyCGM2.settings")
