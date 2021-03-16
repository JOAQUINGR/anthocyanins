from django.core.management.base import BaseCommand
from gestion_pedidos.models import (
    Blastp,
    FimoClasification,
    FimoTotals,
    MotifOrder,
    MotifGeneralInfo,
    HomePage,
)
from django.db.models import Q
from csv import DictReader
import re
import requests
import pandas as pd
import numpy as np
import statsmodels.stats.multitest as fdr
from scipy.stats import hypergeom
import json
import math


class Command(BaseCommand):
    def add_arguments(self, parser):

        parser.add_argument(
            "species_names",
            type=str,
            help="The name of all the pathways that you are going to analyce",
        )
        parser.add_argument(
            "genesfile",
            type=str,
            help="file with genemanes",
        )

    def getBlastpResultForSpecies(self, geneNamesQuery, speciesName):
        tabresult = {}
        orthologsList = set()

        for gene in geneNamesQuery:
            # print(gene)
            genes = Blastp.objects.filter(
                Q(specie_name__scientific_name__icontains=speciesName)
                & Q(arabidopsis_gene__gene_name__icontains=gene)
            ).values("specie_gene__gene_name")
            tabresult[gene] = []
            for result in genes:
                tabresult[gene].append(dict(result)["specie_gene__gene_name"])
                orthologsList.add(dict(result)["specie_gene__gene_name"])
        return orthologsList

    def loadParameters(self, jsonWithParameters):
        with open(jsonWithParameters) as j:
            jsonParameters = json.load(j)
            return jsonParameters

    def getSortedMotifNames(self):

        databaseQueryForThisTable = MotifOrder.objects.filter(id=1).values(
            "orederedmotifnames"
        )
        return databaseQueryForThisTable[0]["orederedmotifnames"]

    def removeEmptyRows(self, pandasDataFrame):
        return pandasDataFrame.loc[:, (pandasDataFrame != 0).any(axis=0)]

    def generateDataFrameForGenesValues(self, geneValues):
        listOfNameLabels = self.getSortedMotifNames()
        df = pd.DataFrame.from_dict(
            geneValues, columns=listOfNameLabels, dtype=int, orient="index"
        )
        return self.removeEmptyRows(df)

    def getTfbsMatchsForEachGeneInDataFrame(self, preferences, geneSet):
        databaseQueryForThisTable = (
            f"all_motifs__{preferences['upsize']}__{ preferences['score']}"
        )
        genesMatchValuesDic = {}
        for gene in geneSet:
            # añado el i contains para que funcione a pesar de que los nombres esten incompletos
            geneMatchsArray = FimoClasification.objects.filter(
                speciegene__gene_name__icontains=gene
            ).values(databaseQueryForThisTable)
            try:
                if len(dict(geneMatchsArray[0])[databaseQueryForThisTable]) > 0:
                    genesMatchValuesDic[gene] = dict(geneMatchsArray[0])[
                        databaseQueryForThisTable
                    ]
                else:
                    print(
                        "There is no info for the gene "
                        + gene
                        + "with this parameters in Fimo"
                    )
            except:
                print("There is no info for the gene " + gene + "in Fimo")

        return self.generateDataFrameForGenesValues(genesMatchValuesDic)

    def getTotalGeneNumberForThisGenome(self, speciesName):

        databaseQueryForThisTable = HomePage.objects.filter(
            specie_name__scientific_name=speciesName
        ).values("total_genes")
        return dict(databaseQueryForThisTable[0])["total_genes"]

    def labelValuesWithNames(self, motifsValues, motifsMames):
        motifTotalsDict = {}
        for motif, totalMatchesNumber in zip(motifsMames, motifsValues):
            motifTotalsDict[motif] = totalMatchesNumber
        return motifTotalsDict

    def getTotalMatchesPerMotifInThisGenome(self, preferences, speciesName):

        databaseQueryForThisTable = (
            f"totalvalues__{preferences['upsize']}__{preferences['score']}"
        )
        totalMatchesForParameters = FimoTotals.objects.filter(
            specie_name__scientific_name=speciesName
        ).values(databaseQueryForThisTable)
        listOfTotalMotifValues = dict(totalMatchesForParameters[0])[
            databaseQueryForThisTable
        ]
        listOfNameLabels = self.getSortedMotifNames()

        return self.labelValuesWithNames(listOfTotalMotifValues, listOfNameLabels)

    def getMotifsExtraInfo(self):
        with open("motifinfo.json", "r") as fp:
            data = json.load(fp)
            return pd.DataFrame.from_dict(data, orient="index")

    def totalGenesInSample(self, genesDataFrame):
        return len(genesDataFrame)

    def getDataFrameOutput(self, genesdf, totalMatchesMotif, totalGenes):
        motifNameList = list(genesdf.columns)
        firstResultsList = []

        total_sample = self.totalGenesInSample(genesdf)
        total_in_genome = totalGenes
        for queryMotif in motifNameList:
            total_maches_insample = genesdf[queryMotif].sum()  # x
            total_motif_genome = totalMatchesMotif[queryMotif]  # n
            ###########Calculate the hipergeometric
            pvalues = hypergeom.pmf(  # x,M,n,N
                total_maches_insample,
                total_in_genome,
                total_motif_genome,
                total_sample,
            )
            ###########Calculate Fold enrichment as log2((x/N)/(n/M))
            FoldEn = math.log2(
                (total_maches_insample / total_sample)
                / (total_motif_genome / total_in_genome)
            )
            firstResultsList.append(  # add a new result to a list of each motif result
                [
                    queryMotif,
                    total_maches_insample,
                    total_sample,
                    total_motif_genome,
                    total_in_genome,
                    pvalues,
                    FoldEn,
                ]
            )

        finaldf = pd.DataFrame(
            firstResultsList,
            columns=[
                "Motif_ID",
                "x",
                "N",
                "n",
                "M",
                "pvalues",
                "FoldEn",
            ],
        )
        """Calculate the adjusted pval"""
        fdrinput = finaldf["pvalues"].to_numpy()
        fdrinput = fdrinput.transpose()
        calculo_fdr = fdr.fdrcorrection(
            fdrinput,
            alpha=float(1),
            method="indep",
            is_sorted=False,
        )
        finaldf["FDR"] = list(calculo_fdr[1])
        finaldf["clasification"] = list(calculo_fdr[0])
        del calculo_fdr
        del fdrinput
        finaldf = pd.merge(finaldf, self.getMotifsExtraInfo(), on="Motif_ID")

        return finaldf.sort_values(by=["FDR"], ascending=True, ignore_index=True)

    def getMYBsOnly(self, totalDataFrame):
        return totalDataFrame[totalDataFrame["tffamily"] == "Myb/SANT"]

    def getOnlyTfNameAndvalue(self, valueColumn, mybs):
        return mybs[[valueColumn, "tf"]].sort_values(by=[valueColumn], ascending=False)

    def processingQueryInDatabaseAndOutput(
        self,
        speciesName,
        geneSet,
        parametersDic,
    ):

        """Perfomm a search using fimo for the parameters 1000 bp and 3000 bq
        and redirect to the generation of the files for the MYB factors which are
        the one that matters"""

        for parametersCase in parametersDic["parameters"]:
            # print(parametersDic)
            # exit()
            print("Computing search with this parameters: ", parametersCase)
            genesMatchValuesDataFrame = self.getTfbsMatchsForEachGeneInDataFrame(
                parametersCase, geneSet
            )
            totalGenesNumber = self.getTotalGeneNumberForThisGenome(speciesName)
            totalMatchesMotifGenome = self.getTotalMatchesPerMotifInThisGenome(
                parametersCase, speciesName
            )
            totalValuesAndCalculations = self.getDataFrameOutput(
                genesMatchValuesDataFrame, totalMatchesMotifGenome, totalGenesNumber
            )
            print(totalValuesAndCalculations.head(20))
            onlyMYBs = self.getMYBsOnly(totalValuesAndCalculations)

            relevantRows = self.getOnlyTfNameAndvalue("FoldEn", onlyMYBs)
            print(relevantRows.head())

            totalValuesAndCalculations.to_csv(
                "{}_{}_{}_{}.csv".format(
                    speciesName,
                    "fimo",
                    parametersCase["score"],
                    parametersCase["upsize"],
                )
            )

    def handle(self, *args, **options):

        print("program start")

        genesfileInfo = options["genesfile"]
        speciesFileInfo = options["species_names"]
        speciesName = []
        with open(speciesFileInfo) as species:
            for selected in species:
                speciesName.append(selected.strip())
        geneInfo = []
        with open(genesfileInfo) as AtOrthologs:
            for gene in AtOrthologs:
                geneInfo.append(gene.strip())

        parametersDic = self.loadParameters("parameterAutomaticBlastp.json")

        for targetSpeciesName in speciesName:
            print(".....................")
            print("Computing Blastp for " + targetSpeciesName)

            orthologsSet = self.getBlastpResultForSpecies(
                geneInfo,
                targetSpeciesName,
            )
            print(
                "A total of ",
                len(orthologsSet),
                " orthologs found for ",
                len(geneInfo),
                " Arabidopsis genes",
            )
            self.processingQueryInDatabaseAndOutput(
                targetSpeciesName, orthologsSet, parametersDic
            )

            # TODO import parameters of search, launch query
            # print(orthologsSet)
            print(".....................")
