# coding: utf-8
from __future__ import unicode_literals
import pyCGM2f
import logging
logging.basicConfig(level=logging.INFO)
import os

def main():

    DATA_PATH = os.getcwd()+"\\" #C:\\Users\\fleboeuf\\Documents\\Programmation\\pyCGM2\\pyCGM2-extensions\\pyCGM2_flow\\data\\cgm1-nerveBlock\Session 2\\"

    print("=================CHECKING==================")

    # Videos
    print("--- Videos ----")
     #
    VIDEO_PATH = DATA_PATH+"Videos\\"
    folders = os.listdir(VIDEO_PATH)
    if len(folders) == 0:
        logging.warning(" !!! Main video folder EMPTY !! ")
    else:
        for foldIt in folders:
            fileLst = os.listdir(VIDEO_PATH+foldIt)
            if fileLst == []:
                logging.warning("Video subfolder [%s] EMPTY"%(foldIt))
            else:
                logging.info("Video subfolder [%s] - %i files detected "%(foldIt, len(fileLst)))


    print("--- Images ----")
    IMAGE_PATH = DATA_PATH+"Images\\"
    folders = os.listdir(IMAGE_PATH)
    if len(folders) == 0:
        logging.warning(" !!! Main image  folder EMPTY !! ")
    else:
        for foldIt in folders:
            fileLst = os.listdir(VIDEO_PATH+foldIt)
            if fileLst == []:
                logging.warning("Video subfolder [%s] EMPTY"%(foldIt))
            else:
                logging.info("Video subfolder [%s] - %i files detected "%(foldIt, len(fileLst)))


    print("--- Exams ----")
    # exam
    EXAM_PATH = DATA_PATH+"Exams\\"
    if "CGA.xlsx" not in os.listdir(EXAM_PATH):
        logging.error("The Clinical Exam file (CGA.xls) has not been imported in the Exam folder")


    print("--- gaitReport ----")
    DOC_PATH = DATA_PATH+"Doc\\"
    if "gaitReport.pdf" not in os.listdir(DOC_PATH):
        logging.error("there is no gaitReport.pdf in the doc folder")



if __name__ == '__main__':
    main()
