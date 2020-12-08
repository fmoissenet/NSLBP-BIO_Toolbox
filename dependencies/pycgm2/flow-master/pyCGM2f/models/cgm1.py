# coding: utf-8

import os
import argparse
import logging
from  collections import OrderedDict


import pyCGM2
from pyCGM2 import enums
from pyCGM2.Utils import files
from pyCGM2.Tools import btkTools
from pyCGM2.Lib.CGM import cgm1
from pyCGM2.Lib import analysis

import pyCGM2f
from pyCGM2f import pipelineProcedure
from pyCGM2f import pipelineFilter
from pyCGM2f.models import common
from pyCGM2f.processing import processing, gaitProcedure


def main(reprocess,userFile = "CGM1.userSettings" , processingOnly = False, expertFile="CGM1-pyCGM2.settings",emgFile="emg.settings",vskFile=None,plotFlag=False, overwrittenMode=True):


    modelInfo={"Version":"cgm1"}

    DATA_PATH = os.getcwd()+"\\"
    DATA_PATH_OUT = DATA_PATH+"Processing-CGM1\\"
    files.createDir(DATA_PATH_OUT)

    if not reprocess:
        userSettings,internalSettings,emgSettings,vsk = common.settingsManager(DATA_PATH,userFile,vskFile,expertFile,emgFile)

        builder = pipelineProcedure.CGM1PipelineBuilder(userSettings)
        builder.setLocalInternalSettings(internalSettings)
        builder.setVsk(vsk)
        if emgSettings is not None: builder.setLocalEmgSettings(emgSettings)

        pmf = pipelineFilter.PipelineFilter()
        pmf.setBuilder(builder)
        pmf.build()
        manager = pmf.pipeline

    else:
        previousManager = files.openFile(DATA_PATH_OUT,"AQM-exam.info")
        manager = pipelineFilter.Pipeline(document=previousManager)


    #---- Modelling------
    models = OrderedDict()
    jsonAnalysisInstances=OrderedDict()

    logging.info("=============Calibration=============")
    if not processingOnly:
        for calibrationIt in manager.getNestedCalibrations():

            model,finalAcqStatic = cgm1.calibrate(DATA_PATH,
                calibrationIt["StaticTrial"],
                calibrationIt["Translators"],
                calibrationIt["MP"]["Required"],
                calibrationIt["MP"]["Optional"],
                calibrationIt["LeftFlatFoot"],
                calibrationIt["RightFlatFoot"],
                calibrationIt["HeadFlat"],
                calibrationIt["MarkerDiameter"],
                calibrationIt["PointSuffix"],
                displayCoordinateSystem=True,
                noKinematicsCalculation=True)

            models[str(calibrationIt["ID"])] = model


    for conditionIt in manager.Document["Conditions"]:

        conditionID = conditionIt["ConditionID"]

        if not processingOnly:
            logging.info("=============Fitting=============")
            if conditionIt["FittingTrials"] != []:
                for fittingTrialIt in conditionIt["FittingTrials"]:
                    acqGait = cgm1.fitting(
                        models[str(fittingTrialIt["CalibrationID"])],
                        DATA_PATH,
                        fittingTrialIt["File"],
                        fittingTrialIt["Translators"],
                        fittingTrialIt["MarkerDiameter"],
                        fittingTrialIt["PointSuffix"],
                        fittingTrialIt["Mfpa"],
                        enums.enumFromtext(conditionIt["MomentProjection"],enums.MomentProjection),
                        displayCoordinateSystem=True)


                    if overwrittenMode:
                        btkTools.smartWriter(acqGait, str(DATA_PATH+fittingTrialIt["File"]))
                    else:
                        btkTools.smartWriter(acqGait, str(DATA_PATH_OUT+fittingTrialIt["File"]))

        logging.info("=============EMG Processing=============")
        if conditionIt["EmgTrials"] != []:
            outDataPath =  None if overwrittenMode else DATA_PATH_OUT
            analysis.processEMG(DATA_PATH,
                                manager.getEmgFiles(conditionID),
                                manager.getEmgConfiguration(conditionID,outputType = "dict")["Labels"],
                                highPassFrequencies = conditionIt["EmgSettings"]["Processing"]["BandpassFrequencies"],
                                envelopFrequency= conditionIt["EmgSettings"]["Processing"]["EnvelopLowpassFrequency"],
                                fileSuffix=None,
                                outDataPath=outDataPath)
            logging.info("[pyCGM2f] - Emg trial for the condition [%s]----> Done"%(conditionID))



        logging.info("=============Processing=============")
        if not overwrittenMode: DATA_PATH = DATA_PATH_OUT

        procedure = gaitProcedure.GaitProcessingProcedure(DATA_PATH,DATA_PATH_OUT, modelInfo)
        processFilter = processing.ProcessingFilter(procedure,manager,plotFlag,jsonAnalysisInstances)
        processFilter.process(conditionIt)

        jsonAnalysisInstances = processFilter.getJsonAnalyses()


    logging.info("=============AQM Info Export=============")
    if not reprocess:
        pmf.exportAQMinfos(DATA_PATH_OUT,"CGM 1.0",jsonAnalyses = jsonAnalysisInstances)
    else:
        for index in range(0,len(manager.Document["Conditions"])):
            conditionID = manager.Document["Conditions"][index]["ConditionID"]
            if conditionID in jsonAnalysisInstances.keys():
                manager.Document["Conditions"][index]["DATA"] = jsonAnalysisInstances[conditionID]


        files.saveJson(DATA_PATH_OUT, "AQM-exam.info", manager.Document)
        logging.info("[pyCGM2f] AQM-exam.info overwritten")


    return manager

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='CGM1-pipeline')
    parser.add_argument('-r','--reprocess', action='store_true', help='reprocess data from CGA.info' )
    parser.add_argument('--userFile', type=str, help='userSettings', default="CGM1.userSettings")
    parser.add_argument('--expertFile', type=str, help='Local expert settings',default="CGM1-pyCGM2.settings")
    parser.add_argument('--emgFile', type=str, help='Local emg settings',default="emg.settings")
    parser.add_argument('--vskFile', type=str, help='Local vsk file')
    parser.add_argument('-po','--Processing', action='store_true', help='processing only' )
    parser.add_argument('-p','--Plot', action='store_true', help='plot data' )

    args = parser.parse_args()

    main(args.reprocess,userFile = args.userFile, processingOnly =args.Processing, expertFile=args.expertFile,emgFile=args.emgFile, vskFile=args.vskFile,plotFlag= args.Plot)
