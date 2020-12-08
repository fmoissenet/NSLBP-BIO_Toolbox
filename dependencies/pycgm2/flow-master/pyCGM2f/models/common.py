import os
import logging

from pyCGM2.Eclipse import vskTools
from pyCGM2.Utils import files



from openpyxl import load_workbook


def overloadUserSettingsFromXls(userSettings,xlsFilename,laboSettings):

    # excel
    workbook = load_workbook(filename=xlsFilename)
    sheet = workbook.get_sheet_by_name(laboSettings["Sheet"])


    if sheet[laboSettings["Ipp"]].value is not None:
        userSettings["SubjectInfo"]["Ipp"] = sheet[laboSettings["Ipp"]].value


    if sheet[laboSettings["Name"]].value is not None:
        userSettings["SubjectInfo"]["Name"] = sheet[laboSettings["Name"]].value


    if sheet[laboSettings["FirstName"]].value is not None:
        userSettings["SubjectInfo"]["FirstName"] = sheet[laboSettings["FirstName"]].value



    if sheet[laboSettings["Age"]].value is not None:
        userSettings["VisitInfo"]["Age"] = sheet[laboSettings["Age"]].value


    if sheet[laboSettings["Goal"]].value is not None:
        userSettings["ExamInfo"]["Goal"] = sheet[laboSettings["Goal"]].value

    if sheet[laboSettings["Comments"]].value is not None:
        userSettings["ExamInfo"]["Comments"] = sheet[laboSettings["Comments"]].value



def settingsManager(DATA_PATH,userFile,vskFile,expertFile,emgFile):
    # User Settings
    if os.path.isfile(DATA_PATH + userFile):
        userSettings = files.openFile(DATA_PATH,userFile)
    else:
        raise Exception ("user setting file not found")

    if vskFile is not None:
        vsk = vskTools.Vsk(str(DATA_PATH +  vskFile))
    else:
        vsk=None

    # internal (expert) Settings
    if os.path.isfile(DATA_PATH + expertFile):
        logging.info("[pyCGM2-pipe] expert settings found in the data folder")
        internalSettings = files.openFile(DATA_PATH,expertFile)
    else:
        internalSettings = None


    # emgSettings
    if os.path.isfile(DATA_PATH + emgFile):
        logging.info("[pyCGM2-pipe] local emgSettings found in the data folder")
        emgSettings = files.openFile(DATA_PATH,emgFile)
    else:
        emgSettings = None

    return userSettings,internalSettings,emgSettings,vsk
