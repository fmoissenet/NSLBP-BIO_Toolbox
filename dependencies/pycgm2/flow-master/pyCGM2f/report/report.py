# coding: utf-8
from __future__ import unicode_literals

import pyCGM2f
from jinja2 import Environment, PackageLoader, select_autoescape,FileSystemLoader







TEMPLATE_PATH = pyCGM2f.MAIN_APPS_PATH+"templates\\"
ENV = Environment(
    loader=FileSystemLoader(TEMPLATE_PATH),
    autoescape=select_autoescape()
)



class Report(object):
    def __init__(self, name = "gaitReport"):
        self.name = name
        self.Analyses = []
        self.Comparisons = []


    def add_analysisSection(self, sectionTitle):

        analysis = dict()
        analysis["section"] = sectionTitle
        analysis["pages"] = list()
        analysis["pageTitles"] = list()

        self.Analyses.append(analysis)

    def add_pageToAnalysis(self, pos,page,pageTiltle):

        self.Analyses[pos]["pages"].append(page)
        self.Analyses[pos]["pageTitles"].append(pageTiltle)

    def add_comparison(self, sectionTitle):

        comparison = dict()
        comparison["section"] = sectionTitle
        comparison["pages"] = list()
        comparison["pageTitles"] = list()

        self.Comparisons.append(comparison)

    def add_pageToComparison(self, pos,page,pageTiltle):

        self.Comparisons[pos]["pages"].append(page)
        self.Comparisons[pos]["pageTitles"].append(pageTiltle)

    def generate(self, DATA_PATH_OUT,processedInfos=None):

        template = ENV.get_template("gaitReport.md")

        output = template.render(data={"Report": self, "infos":processedInfos})

        with open(DATA_PATH_OUT + self.name+".md", 'w') as f:
            f.write(output.encode("utf-8"))
