# coding: utf-8
from __future__ import unicode_literals
import os

import matplotlib.pyplot as plt

import pyCGM2
import pyCGM2f

from pyCGM2.Utils import files

from pyCGM2.Report import normativeDatasets
from pyCGM2.Lib import configuration
from pyCGM2.Lib import plot
from pyCGM2.Lib import analysis


from pyCGM2f.report import report
from pyCGM2f import pipelineFilter


# Report instance
GAITREPORT = report.Report()

def main():

    # CHOICES
    nds = normativeDatasets.Schwartz2008("Free")
    processingFolder = "Processing-CGM23"


    # PATHS
    DATA_PATH =  os.getcwd()+"\\" #"C:\\Users\\fleboeuf\\Documents\\Programmation\\pyCGM2\\pyCGM2-extensions\\pyCGM2_flow\\data\\cgm1-nerveBlock\\Session 2\\"

    PROCESSING_PATH = DATA_PATH+processingFolder+"\\"
    DATA_PATH_OUT = PROCESSING_PATH+"ReportFiles\\"
    files.createDir(DATA_PATH_OUT)

    # Settings
    aqmInfos = files.openFile(PROCESSING_PATH,"AQM-exam.info")
    manager = pipelineFilter.Pipeline(document=aqmInfos)

    intersettings =  files.openFile(DATA_PATH,"emg.settings")
    EMG_LABELS,EMG_MUSCLES,EMG_CONTEXT,NORMAL_ACTIVITIES = configuration.getEmgConfiguration(None, intersettings)

    #___________________________________________________________________________
    # ---- ANALYSES------
    #___________________________________________________________________________

    analyses = []
    conditionIDs = []

    # analysis 1     ------
    #----------------------
    conditionID = "Condition1"
    index  = 1


    GAITREPORT.add_analysisSection("Analyse Marche Spontanée") # add a new section

    analysis1 = files.loadAnalysis(PROCESSING_PATH, conditionID)# analysis object loaded

    # --- pst ---
    fig_pst = plot.plot_spatioTemporal(PROCESSING_PATH, analysis1,  title=conditionID, exportPdf=False,outputName=conditionID,show=False) # plot
    fig_pst.savefig(DATA_PATH_OUT+conditionID+"-stp.png") # save plot
    GAITREPORT.add_pageToAnalysis(index-1, conditionID+"-stp.png" , "Parametres spatio-temporels") #add page to the created section

    # --- kinematic ---
    fig_descpKinematic = plot.plot_DescriptiveKinematic(PROCESSING_PATH,analysis1,"LowerLimb",nds,title=conditionID,show=False)
    fig_descpKinematic.savefig(DATA_PATH_OUT+conditionID+"-descKinematic.png")
    GAITREPORT.add_pageToAnalysis(index-1, conditionID+"-descKinematic.png" , "Cinematique descriptive")


    # --- kinetic ---
    fig_descpKinetic = plot.plot_DescriptiveKinetic(PROCESSING_PATH,analysis1,"LowerLimb",nds,title=conditionID,show=False)
    fig_descpKinetic.savefig(DATA_PATH_OUT+conditionID+"--descKinetic.png")
    GAITREPORT.add_pageToAnalysis(index-1, conditionID+"-descKinetic.png" , "Dynamique descriptive")

    # --- representative emg ---
    #EMG_LABELS,EMG_MUSCLES, EMG_CONTEXT, NORMAL_ACTIVITIES = manager.getEmgConfiguration(conditionID)
    representativeEmgTrial = manager.getRepresentativeEMGTrial(conditionID)
    figs_emg = plot.plotTemporalEMG(DATA_PATH,representativeEmgTrial, EMG_LABELS,EMG_MUSCLES, EMG_CONTEXT, NORMAL_ACTIVITIES,rectify=True,show=False)
    i=0
    for fig in figs_emg:
        filename = conditionID+"-representativeEMG[%i].png"%i
        fig.savefig(DATA_PATH_OUT+filename) # save plot
        GAITREPORT.add_pageToAnalysis(index-1, filename , "Emg Representatif")
        i+=1

    analyses.append(analysis1)
    conditionIDs.append(conditionID)


    # # analysis 2     ------
    # #----------------------
    conditionID = "Condition2"
    index  = 2


    GAITREPORT.add_analysisSection("Analyse Marche Spontanée-Post") # add a new section

    analysis2 = files.loadAnalysis(PROCESSING_PATH, conditionID)# analysis object loaded

    # --- pst ---
    fig_pst = plot.plot_spatioTemporal(PROCESSING_PATH, analysis2,  title=conditionID, exportPdf=False,outputName=conditionID,show=False) # plot
    fig_pst.savefig(DATA_PATH_OUT+conditionID+"-stp.png") # save plot
    GAITREPORT.add_pageToAnalysis(index-1, conditionID+"-stp.png" , "Parametres spatio-temporels") #add page to the created section

    # --- kinematic ---
    fig_descpKinematic = plot.plot_DescriptiveKinematic(PROCESSING_PATH,analysis2,"LowerLimb",nds,title=conditionID,show=False)
    fig_descpKinematic.savefig(DATA_PATH_OUT+conditionID+"-descKinematic.png")
    GAITREPORT.add_pageToAnalysis(index-1, conditionID+"-descKinematic.png" , "Cinematique descriptive")


    # --- kinetic ---
    fig_descpKinetic = plot.plot_DescriptiveKinetic(PROCESSING_PATH,analysis2,"LowerLimb",nds,title=conditionID,show=False)
    fig_descpKinetic.savefig(DATA_PATH_OUT+conditionID+"--descKinetic.png")
    GAITREPORT.add_pageToAnalysis(index-1, conditionID+"-descKinetic.png" , "Dynamique descriptive")

    # --- representative emg ---
    #EMG_LABELS,EMG_MUSCLES, EMG_CONTEXT, NORMAL_ACTIVITIES = manager.getEmgConfiguration(conditionID)
    representativeEmgTrial = manager.getRepresentativeEMGTrial(conditionID)
    figs_emg = plot.plotTemporalEMG(DATA_PATH,representativeEmgTrial, EMG_LABELS,EMG_MUSCLES, EMG_CONTEXT, NORMAL_ACTIVITIES,rectify=True,show=False)
    i=0
    for fig in figs_emg:
        filename = conditionID+"-representativeEMG[%i].png"%i
        fig.savefig(DATA_PATH_OUT+filename) # save plot
        GAITREPORT.add_pageToAnalysis(index-1, filename , "Emg Representatif")
        i+=1


    analyses.append(analysis2)
    conditionIDs.append(conditionID)


    #___________________________________________________________________________
    # ---- COMPARISONS------
    #___________________________________________________________________________


    # ---- comparaison Intra-Session------
    #-------------------------------------

    # *** Comparaison Intra 1***

    name = "Comparaison PRE vs POST"
    shortName = "preVsPost"
    legends =["C1","C2"]
    analysisToCompare = [analyses[0], analyses[1]]
    index = 1

    # --- kinematic ---
    GAITREPORT.add_comparison(name)

    fig = plot.compareKinematic(analysisToCompare,legends,"Left","LowerLimb",nds,plotType="Descriptive",type="Gait",pointSuffixes=None,show=False,title=shortName)
    figName = "ComparaisonIntra"+index+"-Kinematics_Left_"+shortName+".png"
    fig.savefig(DATA_PATH_OUT+figName)
    GAITREPORT.add_pageToComparison(index-1,figName, "cinematique-Gauche")

    fig = plot.compareKinematic(analysisToCompare,legends,"Right","LowerLimb",nds,plotType="Descriptive",type="Gait",pointSuffixes=None,show=False,title=shortName)
    fig = "ComparaisonIntra"+index+"-Kinematics_Right_"+shortName+".png"
    fig.savefig(DATA_PATH_OUT+figName)
    GAITREPORT.add_pageToComparison(index-1,figName, "cinematique-Droite")

    #--- emg ----
    fig = plot.compareEmgEnvelops(analysisToCompare,
                            legends,
                          EMG_LABELS,EMG_MUSCLES,EMG_CONTEXT,NORMAL_ACTIVITIES,
                          normalized=True,
                          plotType="Descriptive",show=False)
    figName = "ComparaisonIntra"+index+"-Emg_"+shortName+".png"
    fig.savefig(DATA_PATH_OUT+figName)

    GAITREPORT.add_pageToComparison(index-1,figName,"comparaison EMG")


    # ---- comparaison Inter-Session------
    #-------------------------------------


    # *** Comparaison Inter 1***

    processingPath_previousSession = "??\\Processing-CGM11\\"
    analysisName = "Condition1"
    previousAnalysis = files.loadAnalysis(processingPath_previousSession, analysisName)

    name = "Comparaison avec session Precedente"
    shortName = "ActualVsPrev"
    legends =["Act","prev"]
    analysisToCompare = [analyses[0], previousAnalysis]
    index = 1

    GAITREPORT.add_comparison(name)

    # --- kinematic ---
    fig = plot.compareKinematic(analysisToCompare,legends,"Left","LowerLimb",nds,plotType="Descriptive",type="Gait",pointSuffixes=None,show=False,title=shortName)
    figName = "ComparaisonInter"+index+"-Kinematics_Left_"+shortName+".png"
    fig.savefig(DATA_PATH_OUT+figName)
    GAITREPORT.add_pageToComparison(index-1,figName,"cinematique-Gauche")

    fig = plot.compareKinematic(analysisToCompare,legends,"Right","LowerLimb",nds,plotType="Descriptive",type="Gait",pointSuffixes=None,show=False,title=shortName)
    figName = "ComparaisonInter"+index+"-Kinematics_Right_"+shortName+".png"
    fig.savefig(DATA_PATH_OUT+figName)
    GAITREPORT.add_pageToComparison(index-1,figName,"cinematique-Gauche")


    #___________________________________________________________________________
    # Generation
    GAITREPORT.generate(DATA_PATH_OUT,processedInfos=aqmInfos)

    # export as Docx
    files.copyPaste(pyCGM2f.MAIN_APPS_PATH+"templates\\"+"wordTemplate.docx",DATA_PATH_OUT+"wordTemplate.docx")
    os.chdir(DATA_PATH_OUT)
    command = "pandoc  --reference-doc=wordTemplate.docx  \"%s\" -o \"%s\" "%("gaitReport.md", DATA_PATH+"Doc//gaitReport.docx")
    os.system(command)

    os.startfile( DATA_PATH+"Doc//gaitReport.docx")


if __name__ == '__main__':
    main()
