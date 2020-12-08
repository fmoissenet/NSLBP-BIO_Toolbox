# coding: utf-8
import argparse
import logging
import os

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

import pyCGM2
from pyCGM2 import log; log.setLogger(level=logging.INFO)
from pyCGM2.Utils import files

import pyCGM2f
from pyCGM2f.models import cgm1, cgm11,cgm21,cgm22,cgm23,cgm24,cgm25
from pyCGM2f.eclipse import eclipse
from pyCGM2f.Init import init
from pyCGM2f.check import check
from pyCGM2f.pushDB import pushDB


MODELS =["CGM1","CGM11", "CGM21",  "CGM22", "CGM23", "CGM24", "CGM25"]

def initiatingCommand():
    """ Initiation of the data folder


    :param -m, --model [str]:  CGM model (choice is CGM1 CGM11 CGM21  CGM22 CGM23 CGM24 CGM25)
    :param -es, --expert [bool]:  copy paste the file CGMi-pycgm2.settings into the data folder

    Examples:

        >>> pycgm2f-init.exe
        >>> pycgm2f-init.exe --expert -m CGM21
    """

    parser = argparse.ArgumentParser(description='pyCGM2-flow-init')
    parser.add_argument('-m','--model', type=str, required = False, help="choice is CGM1 CGM11 CGM21  CGM22 CGM23 CGM24 CGM25")
    parser.add_argument('-es','--expert', action='store_true',  help="fetch expert settings")
    args = parser.parse_args()

    if args.expert and args.model is None:
        raise Exception("[pycgm2f] - you need to select a model as input argument ex: -m CGM1")

    init.main( args.expert, args.model)



def editingCommand():

    """ Edition of the user settings

    :param -m, --model [str] - REQUIRED -:  CGM model (choice is CGM1 CGM11 CGM21  CGM22 CGM23 CGM24 CGM25)
    :param -g, --guru [bool]:  copy paste the guru settings instead of basic settings
    :param -e, --eclipse [bool]:  (NEXUS user only) - automation of the userSettings from eclipse scheme

    Examples:

        >>> pycgm2f-edit.exe -m CGM21
        equivalent to
        >>> pycgm2f-edit.exe --model CGM21
    """

    parser = argparse.ArgumentParser(description='pyCGM2-EclipseInterface')
    parser.add_argument('-e','--eclipse', action='store_true',  help="edit from Vicon Eclipse")
    parser.add_argument('-m','--model', type=str, required = True, help="choice is CGM1 CGM11 CGM21  CGM22 CGM23 CGM24 CGM25")
    parser.add_argument('-g','--guru', action='store_true',  help="guru user-settings")

    args = parser.parse_args()

    if args.model not in MODELS:
        raise Exception ("[pyCGM2f] Model not known. Choice is CGM1, CGM11, CGM21.... CGM25")
    else:
        if args.eclipse:
            eclipse.main(args.model)
        else:
            modelSettingsTemplate =  args.model + ".userSettings"
            if args.guru:
                files.copyPaste(pyCGM2f.MAIN_APPS_PATH+"userSettings\\guru\\"+modelSettingsTemplate,os.getcwd()+"\\"+modelSettingsTemplate)
            else:
                files.copyPaste(pyCGM2f.MAIN_APPS_PATH+"userSettings\\basic\\"+modelSettingsTemplate,os.getcwd()+"\\"+modelSettingsTemplate)


def processingCommand():

    """ Process your data according the userSettings file


    :param -m, --model [str] - REQUIRED -:  CGM model (choice is CGM1 CGM11 CGM21  CGM22 CGM23 CGM24 CGM25)
    :param  --userFile [str]:  force the name of the userSettings file
    :param  --expertFile [str]:  force the name of the expert Settings file
    :param  --vskFile [str]:  force the name of the vsk file (Vicon User only)
    :param -po, --ProcessingOnly [bool]:  No modelling processing. Just handle model ouputs written in your c3d files
    :param -p, --Plot [bool]:  enable plot panels
    :param -now, --NoOverWritten [bool]:  disable the overwritting mode. Data are written into a new c3d
    :param -r, --reprocess [bool]:  reprocess the data from the AQM-exam.info file generated from a previous processing

    Examples:

        >>> pycgm2f-process.exe -m CGM21
        >>> pycgm2f-process.exe --model CGM21
        >>> pycgm2f-process.exe --model CGM21 --plot
        >>> pycgm2f-process.exe --model CGM21 --plot --usersettings "CGM1-altered.Usersettings"

    """

    parser = argparse.ArgumentParser(description='CGM1-pipeline')
    parser.add_argument('-m','--model', type=str, required = True, help="choice is CGM1 CGM11 CGM21  CGM22 CGM23 CGM24 CGM25")
    parser.add_argument('--userFile', type=str, help='userSettings')
    parser.add_argument('--expertFile', type=str, help='Local expert settings')
    parser.add_argument('--emgFile', type=str, help='Local emg settings')
    parser.add_argument('--vskFile', type=str, help='Local vsk file')
    parser.add_argument('-po','--ProcessingOnly', action='store_true', help='processing only' )
    parser.add_argument('-p','--Plot', action='store_true', help='plot data' )
    parser.add_argument('-now','--NoOverWritten', action='store_true', help='overwite c3d' )
    parser.add_argument('-r','--reprocess', action='store_true', help='reprocess data from CGA.info' )


    args = parser.parse_args()

    overwrittenModeFlag = True if not args.NoOverWritten else False

    if args.model == "CGM1":
        if args.userFile is None: args.userFile = "CGM1.userSettings"
        if args.expertFile is None: args.expertFile = "CGM1-pyCGM2.settings"
        if args.emgFile is None: args.emgFile = "emg.settings"
        cgm1.main(args.reprocess,userFile = args.userFile, expertFile=args.expertFile,emgFile=args.emgFile, vskFile=args.vskFile,
            processingOnly =args.ProcessingOnly,
            plotFlag= args.Plot, overwrittenMode = overwrittenModeFlag)

    elif args.model == "CGM11":
        if args.userFile is None: args.userFile = "CGM11.userSettings"
        if args.expertFile is None: args.expertFile = "CGM1_1-pyCGM2.settings"
        if args.emgFile is None: args.emgFile = "emg.settings"
        cgm11.main(args.reprocess,userFile = args.userFile, expertFile=args.expertFile,emgFile=args.emgFile, vskFile=args.vskFile,
            processingOnly =args.ProcessingOnly,
            plotFlag= args.Plot, overwrittenMode = overwrittenModeFlag)

    elif args.model == "CGM21":
        if args.userFile is None: args.userFile = "CGM21.userSettings"
        if args.expertFile is None: args.expertFile = "CGM2_1-pyCGM2.settings"
        if args.emgFile is None: args.emgFile = "emg.settings"
        cgm21.main(args.reprocess,userFile = args.userFile,  expertFile=args.expertFile,emgFile=args.emgFile, vskFile=args.vskFile,
        processingOnly =args.ProcessingOnly,
        plotFlag= args.Plot, overwrittenMode = overwrittenModeFlag)

    elif args.model == "CGM22":
        if args.userFile is None: args.userFile = "CGM22.userSettings"
        if args.expertFile is None: args.expertFile = "CGM2_2-pyCGM2.settings"
        if args.emgFile is None: args.emgFile = "emg.settings"
        cgm22.main(args.reprocess,userFile = args.userFile,  expertFile=args.expertFile,emgFile=args.emgFile, vskFile=args.vskFile,
            processingOnly =args.ProcessingOnly,
            plotFlag= args.Plot, overwrittenMode = overwrittenModeFlag)

    elif args.model == "CGM23":
        if args.userFile is None: args.userFile = "CGM23.userSettings"
        if args.expertFile is None: args.expertFile = "CGM2_3-pyCGM2.settings"
        if args.emgFile is None: args.emgFile = "emg.settings"
        cgm23.main(args.reprocess,userFile = args.userFile, expertFile=args.expertFile,emgFile=args.emgFile, vskFile=args.vskFile,
            processingOnly =args.ProcessingOnly,
            plotFlag= args.Plot, overwrittenMode = overwrittenModeFlag)

    elif args.model == "CGM24":
        if args.userFile is None: args.userFile = "CGM24.userSettings"
        if args.expertFile is None: args.expertFile = "CGM2_4-pyCGM2.settings"
        if args.emgFile is None: args.emgFile = "emg.settings"
        cgm24.main(args.reprocess,userFile = args.userFile,  expertFile=args.expertFile,emgFile=args.emgFile, vskFile=args.vskFile,
            processingOnly =args.ProcessingOnly,
            plotFlag= args.Plot, overwrittenMode = overwrittenModeFlag)

    elif args.model == "CGM25":
        if args.userFile is None: args.userFile = "CGM25.userSettings"
        if args.expertFile is None: args.expertFile = "CGM2_5-pyCGM2.settings"
        if args.emgFile is None: args.emgFile = "emg.settings"
        cgm25.main(args.reprocess,userFile = args.userFile, expertFile=args.expertFile,emgFile=args.emgFile, vskFile=args.vskFile,
            processingOnly =args.ProcessingOnly,
            plotFlag= args.Plot, overwrittenMode = overwrittenModeFlag)

    else:
        raise Exception ("[pyCGM2f] Model not known. Choice is CGM1, CGM11, CGM21.... CGM25")

def reportingCommand():

    parser = argparse.ArgumentParser(description='Reporting')
    args = parser.parse_args()

    os.system("ipython reportGenerator.py")


def checkingCommand():

    parser = argparse.ArgumentParser(description='Checking')
    args = parser.parse_args()

    check.main()



def pushingCommand():

    parser = argparse.ArgumentParser(description='push to database')
    parser.add_argument('--mongoDB', action='store_true',  help="edit from Vicon Eclipse")
    parser.add_argument('-m','--model', type=str, required = True, help="choice is CGM1 CGM11 CGM21  CGM22 CGM23 CGM24 CGM25")
    args = parser.parse_args()

    if args.model not in MODELS:
        raise Exception ("[pyCGM2f] Model not known. Choice is CGM1, CGM11, CGM21.... CGM25")
    else:
        pushDB.main(args.model, mongoDBinsert = args.mongoDB)
