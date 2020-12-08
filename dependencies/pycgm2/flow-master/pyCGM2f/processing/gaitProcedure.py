# coding: utf-8
import logging
import copy
from  collections import OrderedDict
import numpy as np

from pyCGM2.Utils import files
from pyCGM2.Lib import analysis
from pyCGM2.Processing import exporter
from pyCGM2.Report import normativeDatasets
from pyCGM2.Inspect import inspectFilters, inspectProcedures
from pyCGM2.Tools import btkTools
from pyCGM2.Processing import discretePoints
from pyCGM2.EMG import discreteEMGanalysis
from pyCGM2.EMG import normalActivation
from pyCGM2.Processing import cycle
from pyCGM2 import btk


from pyCGM2.Lib import plot

class GaitProcessingProcedure(object):

    def __init__(self,DATA_PATH,DATA_PATH_OUT, modelInfo):
        self.m_DATA_PATH = DATA_PATH
        self.m_DATA_PATH_OUT = DATA_PATH_OUT
        self.m_modelInfo = modelInfo

    def process(self,manager,condition,plotFlag):


        conditionID = condition["ConditionID"]

        EMG_LABELS,EMG_MUSCLES,EMG_CONTEXT,NORMAL_ACTIVITIES = manager.getEmgConfiguration(conditionID)


        logging.info("---Events checking---")
        checkGaitEvents(self.m_DATA_PATH,manager,conditionID)
        logging.info("[pyCGM2f] -Gait Event checking for the condition [%s]----> Done"%(conditionID))


        analysisInstance = analysis.makeCGMGaitAnalysis(self.m_DATA_PATH,
                        manager.getFittingFiles(conditionID),
                        manager.getEmgFiles(conditionID),
                        EMG_LABELS,
                        subjectInfo=manager.getSubjectFirstLevelInfos(),
                        experimentalInfo=manager.getConditionFirstLevelInfos(conditionID),
                        modelInfo=self.m_modelInfo,
                        pointLabelSuffix=condition["Global"]["PointSuffix"])


        logging.info("=============EMG Normalization Processing=============")

        if manager.getEmgFiles(conditionID) != []:
            if condition["EmgReferenceConditionID"] is not None:
                analysisDenominator = files.loadAnalysis(self.m_DATA_PATH_OUT, condition["EmgReferenceConditionID"])
                analysis.normalizedEMG(analysisInstance,EMG_LABELS,EMG_CONTEXT,method="MeanMax", fromOtherAnalysis=analysisDenominator)
            else:
                analysis.normalizedEMG(analysisInstance,EMG_LABELS,EMG_CONTEXT,method="MeanMax", fromOtherAnalysis=None)

        files.saveAnalysis(analysisInstance,self.m_DATA_PATH_OUT,conditionID)

        # report temporal emg features ( values + normal emg activities) in a dict
        temporalEmgDict = None
        if manager.getRepresentativeEMGTrial(conditionID) is not None:
            temporalEmgDict = getTemporalEmgAsDict(self.m_DATA_PATH, manager, conditionID, EMG_LABELS, EMG_CONTEXT,NORMAL_ACTIVITIES)


        #----- Benedetti Processing + export du fichier
        if manager.getFittingFiles(conditionID)!=[]:
            extractBenedetti(analysisInstance,manager,conditionID,self.m_modelInfo,self.m_DATA_PATH_OUT,conditionID)
            logging.info("[pyCGM2f] -Benedetti parameters extracted for the condition [%s]----> Done"%(conditionID))
        else:
            logging.warning("[pyCGM2-flow] - No benedetti computation for the condition [%s]. It is not a  Gait condition "%(conditionID))

        #----- extract EMG + export du fichier
        #----- EMG
        if manager.getEmgFiles(conditionID) != []:
            extractEmgAmplitudes(analysisInstance,manager,conditionID,EMG_LABELS,EMG_MUSCLES,EMG_CONTEXT, self.m_DATA_PATH_OUT,conditionID)
            logging.info("[pyCGM2f] -EMG amplitudes  extracted for the condition [%s]----> Done"%(conditionID))
        else:
            logging.warning("[pyCGM2-flow] - No discrete EMG computation for the condition [%s]. It is not a  Gait condition "%(conditionID))

        #----- Analysis export
        # export xls
        copiedAnalysis = copy.deepcopy(analysisInstance)
        exporter.renameEmgInAnalysis(copiedAnalysis,EMG_LABELS,EMG_MUSCLES,EMG_CONTEXT)
        exportXlsFilter = exporter.XlsAnalysisExportFilter()
        exportXlsFilter.setAnalysisInstance(copiedAnalysis)
        exportXlsFilter.export(conditionID, path=self.m_DATA_PATH_OUT,excelFormat = "xls",mode="Advanced")


        # json file
        exportFilter = exporter.AnalysisExportFilter()
        exportFilter.setAnalysisInstance(analysisInstance)
        jsonContent = exportFilter.export(conditionID, path=self.m_DATA_PATH_OUT)

        if temporalEmgDict is not None:
            jsonContent["RepresentativeTemporalEmg"] = temporalEmgDict


        # plots
        if plotFlag:
            nds = normativeDatasets.Schwartz2008("Free") # normative dataset
            viewPlots(self.m_DATA_PATH, condition,
                        analysisInstance, nds,
                        EMG_LABELS, EMG_MUSCLES, EMG_CONTEXT, NORMAL_ACTIVITIES)

        return jsonContent

def checkGaitEvents(DATA_PATH,manager,conditionID):

    if manager.getFittingFiles(conditionID) !=[]:
        for trial in manager.getFittingFiles(conditionID):
            acq = btkTools.smartReader(DATA_PATH+trial)
            inspectprocedureEvents = inspectProcedures.GaitEventQualityProcedure(acq)
            inspector = inspectFilters.QualityFilter(inspectprocedureEvents,verbose=False)
            inspector.run()
            if not inspector.getState():
                raise Exception("[pyCGM2-flow] - gait events not correctly detected in fitting trial  [%s]"%(trial))

    if manager.getEmgFiles(conditionID) !=[]:
        for trial in manager.getEmgFiles(conditionID):
            acq = btkTools.smartReader(DATA_PATH+trial)
            inspectprocedureEvents = inspectProcedures.GaitEventQualityProcedure(acq)
            inspector = inspectFilters.QualityFilter(inspectprocedureEvents,verbose=False)
            inspector.run()
            if not inspector.getState():
                raise Exception("[pyCGM2-flow] - gait events not correctly detected in emg trial  [%s]"%(trial))


def extractBenedetti(analysisInstance,manager,conditionID,modelInfo,DATA_PATH_OUT,filename):
    #----- Benedetti Processing + export du fichier
    dpProcedure = discretePoints.BenedettiProcedure()
    dpf = discretePoints.DiscretePointsFilter(dpProcedure, analysisInstance,
                modelInfo=modelInfo,
                subjInfo=manager.getSubjectFirstLevelInfos(),
                condExpInfo=manager.getConditionFirstLevelInfos(conditionID))
    dataFrame = dpf.getOutput()

    xlsExport = exporter.XlsExportDataFrameFilter()
    xlsExport.setDataFrames(dataFrame)
    xlsExport.export(filename+"-Benedetti", path=DATA_PATH_OUT)


def extractEmgAmplitudes(analysisInstance,manager,conditionID,EMG_LABELS,EMG_MUSCLES,EMG_CONTEXT, DATA_PATH_OUT,filename):
    emgProcedure = discreteEMGanalysis.AmplitudesProcedure()
    filter = discreteEMGanalysis.DiscreteEMGFilter(emgProcedure,analysisInstance,EMG_LABELS,EMG_MUSCLES,EMG_CONTEXT,
                        subjInfo=manager.getSubjectFirstLevelInfos(),
                        condExpInfo=manager.getConditionFirstLevelInfos(conditionID))
    emgDataframe = filter.getOutput()


    xlsExport = exporter.XlsExportDataFrameFilter()
    xlsExport.setDataFrames(emgDataframe)
    xlsExport.export(filename+"-EMG", path=DATA_PATH_OUT)


def getTemporalEmgAsDict(DATA_PATH,manager,conditionID,EMG_LABELS,EMG_CONTEXT,NORMAL_ACTIVITIES):

    out = OrderedDict()

    representativeEmgTrial = manager.getRepresentativeEMGTrial(conditionID)
    acq =btkTools.smartReader(DATA_PATH+representativeEmgTrial)
    events = acq.GetEvents()

    index = 0
    for itLabel in EMG_LABELS:

        values= acq.GetAnalog(itLabel).GetValues()
        normalValues = np.zeros((values.shape[0]))
        footStrikes = list()
        footOffs = list()

        for it in btk.Iterate(events):
            if it.GetLabel() == "Foot Strike" and it.GetContext() == EMG_CONTEXT[index]:
                footStrikes.append((it.GetFrame()-acq.GetFirstFrame())*acq.GetNumberAnalogSamplePerFrame())
            if it.GetLabel() == "Foot Off" and it.GetContext() == EMG_CONTEXT[index]:
                footOffs.append((it.GetFrame()-acq.GetFirstFrame())*acq.GetNumberAnalogSamplePerFrame())

        try:
            gaitCycles = cycle.construcGaitCycle(acq)
            for cycleIt  in gaitCycles:
                if cycleIt.context == EMG_CONTEXT[index]:
                    pos,burstDuration=normalActivation.getNormalBurstActivity_fromCycles(NORMAL_ACTIVITIES[index],cycleIt.firstFrame,cycleIt.begin, cycleIt.m_contraFO, cycleIt.end, cycleIt.appf)
                    for i in range(0,len(pos)):
                        normalValues[int(pos[i]):int(pos[i]+burstDuration[i])] = 1
        except:
            logging.warning("[pyCGM2f] - Normal emg from gait events not calculated")

        out[itLabel] = {
            "Values": values.tolist(),
            "NormalValues":  normalValues.tolist(),
            "FootStrikes" : footStrikes ,
            "footOffs" : footOffs }

        index+=1

    return out

def viewPlots(DATA_PATH,condition,analysisInstance,nds, EMG_LABELS,EMG_MUSCLES,EMG_CONTEXT, NORMAL_ACTIVITIES):
    if condition["EmgTrials"] != []:
        plot.plotConsistencyEnvelopEMGpanel(DATA_PATH,analysisInstance, EMG_LABELS,EMG_MUSCLES,EMG_CONTEXT, NORMAL_ACTIVITIES,
                                title=condition["ConditionID"],
                                normalized=True)

    plot.plot_DescriptiveKinematic(DATA_PATH,analysisInstance,"LowerLimb",
            nds,
            pointLabelSuffix=None,
            exportPdf=False,outputName=None,show=True,title=condition["ConditionID"])

    plot.plot_DescriptiveKinetic(DATA_PATH,analysisInstance,"LowerLimb",
            nds,
            pointLabelSuffix=None,
            exportPdf=False,outputName=None,show=True,title=condition["ConditionID"])
