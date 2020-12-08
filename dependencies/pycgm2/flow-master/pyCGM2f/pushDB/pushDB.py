# coding: utf-8
# from __future__ import unicode_literals
import logging
logging.basicConfig(level=logging.INFO)
import os
import shutil
import json

from pyCGM2.Utils import files

import pyCGM2f
from pyCGM2f import pipelineFilter

from pyCGM2f.mongodb import connection
# from pprintpp import pprint

DB_MEDIA_GAITDATA_PATH = pyCGM2f.DB_MEDIA_PATH+"gaitdata\\"

def fixEmgChannels_forMongo(aqmInfos):
    for index in range(0, len(aqmInfos["Conditions"])):

        for key in aqmInfos["Conditions"][index]["EmgSettings"]["CHANNELS"]:
            if "." in key:
                flag = True
                aqmInfos["Conditions"][index]["EmgSettings"]["CHANNELS"][key.replace(".","-")] =  aqmInfos["Conditions"][index]["EmgSettings"]["CHANNELS"].pop(key)

        for key in aqmInfos["Conditions"][index]["DATA"].keys():
            if "." in key:
                flag = True
                aqmInfos["Conditions"][index]["DATA"][key.replace(".","-")] =  aqmInfos["Conditions"][index]["DATA"].pop(key)
            if key =="RepresentativeTemporalEmg":
                for key2 in aqmInfos["Conditions"][index]["DATA"]["RepresentativeTemporalEmg"].keys():
                    if "." in key2:
                        aqmInfos["Conditions"][index]["DATA"]["RepresentativeTemporalEmg"][key2.replace(".","-")] =  aqmInfos["Conditions"][index]["DATA"]["RepresentativeTemporalEmg"].pop(key2)


    if flag: logging.warning("[pyCGM2f] - emg channel labels  \".\" replace by \"-\"  ")

    return aqmInfos


def main(CGMVersion,mongoDBinsert=False):

    processingFolder = "Processing-" + CGMVersion


    DATA_PATH = os.getcwd()+"\\"#"C:\\Users\\fleboeuf\\Documents\\Programmation\\pyCGM2\\pyCGM2-extensions\\pyCGM2_flow\\data\\cgm1-nerveBlock\Session 2\\" #
    PROCESSING_PATH = DATA_PATH+processingFolder+"\\"

    aqmInfos = files.openFile(PROCESSING_PATH,"AQM-exam.info")
    manager = pipelineFilter.Pipeline(document=aqmInfos)

    ipp = manager.Document["Ipp"]
    visitNumber = manager.Document["SessionNumber"]

    DB_MEDIA_GAITDATA_PATH_SUBJECT_VISIT = DB_MEDIA_GAITDATA_PATH  + ipp +"//"+str(visitNumber)+"//"

    # check presence of a subject folder in media
    if ipp in os.listdir(DB_MEDIA_GAITDATA_PATH):
        logging.info("Subject %s already detected in the media folder")
        if str(visitNumber) in os.listdir(DB_MEDIA_GAITDATA_PATH+ipp):
            raise Exception ("[pyCGM2f] STOP - visit #%i of the patient %s already exists"%(visitNumber, ipp))
        else:
            files.createDir(DB_MEDIA_GAITDATA_PATH_SUBJECT_VISIT)
            logging.info( "folder %s // %i created in media"%(ipp,visitNumber))
    else:
        files.createDir(DB_MEDIA_GAITDATA_PATH_SUBJECT_VISIT)
        logging.info( "folder %s // %i created in media"%(ipp,visitNumber))

    # copy paste folders
    shutil.copytree(DATA_PATH+"Doc", DB_MEDIA_GAITDATA_PATH_SUBJECT_VISIT+"Doc")
    shutil.copytree(DATA_PATH+"Videos", DB_MEDIA_GAITDATA_PATH_SUBJECT_VISIT+"Videos")
    shutil.copytree(DATA_PATH+"Images", DB_MEDIA_GAITDATA_PATH_SUBJECT_VISIT+"Images")
    shutil.copytree(DATA_PATH+"Exams", DB_MEDIA_GAITDATA_PATH_SUBJECT_VISIT+"Exams")
    shutil.copytree(PROCESSING_PATH, DB_MEDIA_GAITDATA_PATH_SUBJECT_VISIT+processingFolder)

    if mongoDBinsert:
        # update mongoDB - CGA
        manDB = connection.MongoConnection()
        cursor_one = manDB.get_collection('CGA').find_one({"Ipp":ipp,"SessionNumber":visitNumber})

        if cursor_one is not None:
            raise Exception("[pyCGM2f] STOP - visit #%i of the patient %s  already stored in the database"%(visitNumber, ipp))
        else:
            aqmInfos = fixEmgChannels_forMongo(aqmInfos) # remove . from fieldname
            manDB.get_collection('CGA').insert_one(aqmInfos)
            logging.info(" aqm infos pushed into the database")




if __name__ == '__main__':
    main() #main("CGM1",mongoDBinsert=True)

    # [ name for name in os.listdir(thedir) if os.path.isdir(os.path.join(thedir, name)) ]
