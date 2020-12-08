# coding: utf-8
# pytest -s --disable-pytest-warnings  test_CGM23Settings.py::Test_CGM23::test_GuruSettings
# coding: utf-8
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
import os

from pyCGM2.Utils import files

from pyCGM2f import pipelineProcedure
from pyCGM2f import pipelineFilter

import pytest


class Test_CGM23:
    def test_basicSettings(self):

        DATA_PATH = os.getcwd()+"\\settings\\cgm23\\"
        # User Settings
        userSettings = files.openFile(DATA_PATH,"CGM23-basic.userSettings")

        internalSettings = None
        emgSettings = None
        vsk=None

        builder = pipelineProcedure.CGM23PipelineBuilder(userSettings)
        builder.setLocalInternalSettings(internalSettings)
        if emgSettings is not None: builder.setLocalEmgSettings(emgSettings)
        builder.setVsk(vsk)

        pmf = pipelineFilter.PipelineFilter()
        pmf.setBuilder(builder)
        pmf.build()
        manager = pmf.pipeline

        files.saveJson(DATA_PATH,"manager-basic.info",manager.Document)

    def test_AdvancedSettings(self):

        DATA_PATH = os.getcwd()+"\\settings\\cgm23\\"
        # User Settings
        userSettings = files.openFile(DATA_PATH,"CGM23-basic.userSettings")
        internalSettings = files.openFile(DATA_PATH,"CGM2_3-pyCGM2.settings")
        emgSettings = files.openFile(DATA_PATH,"emg.settings")
        vsk=None

        builder = pipelineProcedure.CGM23PipelineBuilder(userSettings)
        builder.setLocalInternalSettings(internalSettings)
        if emgSettings is not None: builder.setLocalEmgSettings(emgSettings)
        builder.setVsk(vsk)

        pmf = pipelineFilter.PipelineFilter()
        pmf.setBuilder(builder)
        pmf.build()
        manager = pmf.pipeline

        files.saveJson(DATA_PATH,"manager-advanced.info",manager.Document)

    def test_GuruInheritedSettings(self):

        DATA_PATH = os.getcwd()+"\\settings\\cgm23\\"
        # User Settings
        userSettings = files.openFile(DATA_PATH,"CGM23-guru-inherit.userSettings")
        internalSettings = files.openFile(DATA_PATH,"CGM2_3-pyCGM2.settings")
        emgSettings = files.openFile(DATA_PATH,"emg.settings")
        vsk=None

        builder = pipelineProcedure.CGM23PipelineBuilder(userSettings)
        builder.setLocalInternalSettings(internalSettings)
        if emgSettings is not None: builder.setLocalEmgSettings(emgSettings)
        builder.setVsk(vsk)

        pmf = pipelineFilter.PipelineFilter()
        pmf.setBuilder(builder)
        pmf.build()
        manager = pmf.pipeline

        files.saveJson(DATA_PATH,"manager-guruInherited.info",manager.Document)

    def test_GuruSettings(self):

        DATA_PATH = os.getcwd()+"\\settings\\cgm23\\"
        # User Settings
        userSettings = files.openFile(DATA_PATH,"CGM23-guru-modified.userSettings")
        internalSettings = files.openFile(DATA_PATH,"CGM2_3-pyCGM2.settings")
        emgSettings = files.openFile(DATA_PATH,"emg.settings")
        vsk=None

        builder = pipelineProcedure.CGM23PipelineBuilder(userSettings)
        builder.setLocalInternalSettings(internalSettings)
        if emgSettings is not None: builder.setLocalEmgSettings(emgSettings)
        builder.setVsk(vsk)

        pmf = pipelineFilter.PipelineFilter()
        pmf.setBuilder(builder)
        pmf.build()
        manager = pmf.pipeline

        files.saveJson(DATA_PATH,"manager-guru.info",manager.Document)
