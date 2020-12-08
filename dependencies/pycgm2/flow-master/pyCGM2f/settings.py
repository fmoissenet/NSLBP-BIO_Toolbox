import os
from openpyxl import load_workbook



# excel
if os.path.isfile(DATA_PATH + "CGA.xlsx"):
    workbook = load_workbook(filename=DATA_PATH + "CGA.xlsx")
    coverSheet = workbook.get_sheet_by_name("Cover")

    if coverSheet["B3"].value is not None:
        manager.Subject["Ipp"] = coverSheet["B3"].value
        manager.Document["Ipp"] = coverSheet["B3"].value

    if coverSheet["B4"].value is not None:
        manager.Subject["Ipp"] = coverSheet["B3"].value
        manager.Document["Name"] = coverSheet["B4"].value

    if coverSheet["B5"].value is not None:
        manager.Subject["Ipp"] = coverSheet["B3"].value
        manager.Document["FirstName"] = coverSheet["B5"].value

    if coverSheet["B7"].value is not None:
        manager.Visit["Age"] = coverSheet["B7"].value
        manager.Document["AQM"]["Age"] = coverSheet["B7"].value

    if coverSheet["B22"].value is not None:
        manager.Session["Goal"] = coverSheet["B3"].value
        manager.Document["AQM"]["Goal"] = coverSheet["B22"].value

    if coverSheet["B23"].value is not None:
        manager.Session["Comments"] = coverSheet["B3"].value
        manager.Document["AQM"]["Comments"] = coverSheet["B23"].value
