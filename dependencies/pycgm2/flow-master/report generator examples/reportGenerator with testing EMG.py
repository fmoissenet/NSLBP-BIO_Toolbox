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

def main(export=True):

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

    # analysis 1     ------
    #----------------------
    GAITREPORT.add_analysisSection("Analyse Marche Spontanée") # add a new section

    analysis1 = files.loadAnalysis(PROCESSING_PATH, "Condition1")# analysis object loaded

    # --- pst ---
    fig1_pst = plot.plot_spatioTemporal(PROCESSING_PATH, analysis1,  title="Condition1", exportPdf=False,outputName="Condition1",show=False) # plot
    fig1_pst.savefig(DATA_PATH_OUT+"Condition1-stp.png") # save plot

    GAITREPORT.add_pageToAnalysis(0, "Condition1-stp.png" , "Parametres spatio-temporels") #add page to the created section

    # --- kinematic ---
    fig1_descpKinematic = plot.plot_DescriptiveKinematic(PROCESSING_PATH,analysis1,"LowerLimb",nds,title="Condition1",show=False)
    fig1_descpKinematic.savefig(DATA_PATH_OUT+"Condition1-descKinematic.png")

    GAITREPORT.add_pageToAnalysis(0, "Condition1-descKinematic.png" , "Cinematique descriptive")

    # --- representative emg ---
    #EMG_LABELS,EMG_MUSCLES, EMG_CONTEXT, NORMAL_ACTIVITIES = manager.getEmgConfiguration("Condition1")
    representativeEmgTrial = manager.getRepresentativeEMGTrial("Condition1")
    figs1_emg = plot.plotTemporalEMG(DATA_PATH,representativeEmgTrial, EMG_LABELS,EMG_MUSCLES, EMG_CONTEXT, NORMAL_ACTIVITIES,rectify=True,show=False)
    i=0
    for fig in figs1_emg:
        filename = DATA_PATH_OUT+"Condition1-representativeEMG[%i].png"%i
        fig.savefig(filename) # save plot
        GAITREPORT.add_pageToAnalysis(0, filename , "Emg Representatif")
        i+=1

    # UNCOMMENT IF NECESSARY---------------------------------------------------
    # analysis 2     ------
    #----------------------
    GAITREPORT.add_analysisSection("Analyse Marche Spontanée - Post Bloc RFd")

    analysis2 = files.loadAnalysis(PROCESSING_PATH, "Condition2")

    # --- pst ---
    fig2_pst = plot.plot_spatioTemporal(PROCESSING_PATH, analysis2,  title="Condition1", exportPdf=False,outputName="Condition2",show=False)
    fig2_pst.savefig(DATA_PATH_OUT+"Condition2-stp.png")

    GAITREPORT.add_pageToAnalysis(1, "Condition2-stp.png" , "Spatio-temporal parameter")

    # --- kinematic ---
    fig2_descpKinematic = plot.plot_DescriptiveKinematic(PROCESSING_PATH,analysis2,"LowerLimb",nds,title="Condition2",show=False)
    fig2_descpKinematic.savefig(DATA_PATH_OUT+"Condition2-descKinematic.png")

    GAITREPORT.add_pageToAnalysis(1, "Condition2-descKinematic.png" , "Descriptive Kinematics")

    # --- representative emg ---
    #EMG_LABELS,EMG_MUSCLES, EMG_CONTEXT, NORMAL_ACTIVITIES = manager.getEmgConfiguration("Condition1")
    representativeEmgTrial = manager.getRepresentativeEMGTrial("Condition2")
    figs2_emg = plot.plotTemporalEMG(DATA_PATH,representativeEmgTrial, EMG_LABELS,EMG_MUSCLES, EMG_CONTEXT, NORMAL_ACTIVITIES,rectify=True,show=False)
    i=0
    for fig in figs2_emg:
        filename = DATA_PATH_OUT+"Condition2-representativeEMG[%i].png"%i
        fig.savefig(filename) # save plot
        GAITREPORT.add_pageToAnalysis(1, filename , "Emg Representatif")
        i+=1



    #___________________________________________________________________________
    # ---- COMPARISONS------
    #___________________________________________________________________________


    # ---- comparaison 1------
    #-------------------------

    GAITREPORT.add_comparison("Comparaison PRE vs POST")


    figA1_emg = plot.compareEmgEnvelops([analysis1, analysis2],
                            ["Cond1", "Cond2"],
                          EMG_LABELS,EMG_MUSCLES,EMG_CONTEXT,NORMAL_ACTIVITIES,
                          normalized=True,
                          plotType="Descriptive",show=False)
    figA1_emg.savefig(DATA_PATH_OUT+"ComparaisonEMG-C1vsC2.png")

    GAITREPORT.add_pageToComparison(0,"ComparaisonEMG-C1vsC2.png","comparaisonEMG")

    fig_comp_c1vsC2 = plot.compareKinematic([analysis1, analysis2],["Cond1", "Cond2"],"Left","LowerLimb",nds,plotType="Descriptive",type="Gait",pointSuffixes=None,show=False,title="C1 vs C2")
    fig_comp_c1vsC2.savefig(DATA_PATH_OUT+"Comparaison-C1vsC2.png")

    GAITREPORT.add_pageToComparison(0,"Comparaison-C1vsC2.png","cinematique-Gauche")

    fig_comp_c1vsC2 = plot.compareKinematic([analysis1, analysis2],["Cond1", "Cond2"],"Right","LowerLimb",nds,plotType="Descriptive",type="Gait",pointSuffixes=None,show=False,title="C1 vs C2")
    fig_comp_c1vsC2.savefig(DATA_PATH_OUT+"Comparaison-C1vsC2-Right.png")

    GAITREPORT.add_pageToComparison(0,"Comparaison-C1vsC2-Right.png","cinematique-Droite")


    # ---- comparaison 2 PRE-POST Hip Flexion------
    #---------------------------------------------
    GAITREPORT.add_comparison("Comparaison Testing - Flexion hanche - PRE vs POST")

    # ------ testing 1 --------
    emgTrials_testing1_1= ["20200923-PC-PRE-hip flexion 01.c3d"]

    analysis.processEMG(DATA_PATH, emgTrials_testing1_1, EMG_LABELS,
                        highPassFrequencies=[20,200], envelopFrequency=6,fileSuffix=None)

    analysisTesting1_1 = analysis.makeEmgAnalysis(DATA_PATH,emgTrials_testing1_1,EMG_LABELS,type="testing")
    analysis.normalizedEMG(analysisTesting1_1,EMG_LABELS,EMG_CONTEXT,method="MeanMax", fromOtherAnalysis=None)

    # ------ testing 2 --------
    emgTrials_testing1_2 = ["20200923-PC-POST- hip flexion 01.c3d"]

    analysis.processEMG(DATA_PATH, emgTrials_testing1_2, EMG_LABELS,
                        highPassFrequencies=[20,200], envelopFrequency=6,fileSuffix=None)

    analysisTesting1_2 = analysis.makeEmgAnalysis(DATA_PATH,emgTrials_testing1_2,EMG_LABELS,type="testing")
    analysis.normalizedEMG(analysisTesting1_2,EMG_LABELS,EMG_CONTEXT,method="MeanMax", fromOtherAnalysis=analysisTesting1_1)

    # ------ plot --------
    fig_comp_testing1 = plot.compareEmgEnvelops([analysisTesting1_1,analysisTesting1_2],["Pre", "Post"],EMG_LABELS,EMG_MUSCLES,EMG_CONTEXT,NORMAL_ACTIVITIES,type="testing",show=False)
    fig_comp_testing1.savefig(DATA_PATH_OUT+"Comparaison-testing-HipFlexion.png")
    GAITREPORT.add_pageToComparison(1,"Comparaison-testing-HipFlexion.png","EMG")

    #plot.compareSelectedEmgEvelops([analysisTesting1,analysisTesting2],[labelId1,labelId2],["Voltage.EMG1","Voltage.EMG1"],["Left","Left"],normalized=True)


    # ---- comparaison 2 PRE-POST Hip Flexion------
    #---------------------------------------------
    GAITREPORT.add_comparison("Comparaison Testing - Stretching - PRE vs POST")

    # ------ testing 1 --------
    emgTrials_testing2_1= ["20200923-PC-PRE- stretching RF.c3d"]

    analysis.processEMG(DATA_PATH, emgTrials_testing2_1, EMG_LABELS,
                        highPassFrequencies=[20,200], envelopFrequency=6,fileSuffix=None)

    analysisTesting2_1 = analysis.makeEmgAnalysis(DATA_PATH,emgTrials_testing2_1,EMG_LABELS,type="testing")
    analysis.normalizedEMG(analysisTesting2_1,EMG_LABELS,EMG_CONTEXT,method="MeanMax", fromOtherAnalysis=None)

    # ------ testing 2 --------
    emgTrials_testing2_2 = ["20200923-PC-POST- stretching.c3d"]

    analysis.processEMG(DATA_PATH, emgTrials_testing2_2, EMG_LABELS,
                        highPassFrequencies=[20,200], envelopFrequency=6,fileSuffix=None)

    analysisTesting2_2 = analysis.makeEmgAnalysis(DATA_PATH,emgTrials_testing2_2,EMG_LABELS,type="testing")
    analysis.normalizedEMG(analysisTesting2_2,EMG_LABELS,EMG_CONTEXT,method="MeanMax", fromOtherAnalysis=analysisTesting2_1)

    # ------ plot --------
    fig_comp_testing2 = plot.compareEmgEnvelops([analysisTesting2_1,analysisTesting2_2],["Pre", "Post"],EMG_LABELS,EMG_MUSCLES,EMG_CONTEXT,NORMAL_ACTIVITIES,type="testing",show=False)
    fig_comp_testing2.savefig(DATA_PATH_OUT+"Comparaison-testing-stretching.png")
    GAITREPORT.add_pageToComparison(2,"Comparaison-testing-stretching.png","EMG")






    #___________________________________________________________________________
    # Generation
    GAITREPORT.generate(DATA_PATH_OUT,processedInfos=aqmInfos)

    # export as Docx
    if export:
        files.copyPaste(pyCGM2f.MAIN_APPS_PATH+"templates\\"+"wordTemplate.docx",DATA_PATH_OUT+"wordTemplate.docx")
        os.chdir(DATA_PATH_OUT)
        command = "pandoc  --reference-doc=wordTemplate.docx  \"%s\" -o \"%s\" "%("GAITREPORT.md", DATA_PATH+"GAITREPORT.docx")
        os.system(command)

if __name__ == '__main__':
    main(export=True)
