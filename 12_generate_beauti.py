

def generate_taxon(acc, date, direction, units, **kwargs):
    taxa_string = ""
    taxa_string += '\t\t<taxon id="%s">\n' % acc
    taxa_string += '\t\t\t<date value="%s" direction="%s" units="%s"/>\n' % (date, direction, units)
    for k, v in kwargs.items():
        taxa_string += '\t\t\t<attr name="%s">\n\t\t\t\t%s\n\t\t\t</attr>\n' % (k, v)
    taxa_string += '\t\t</taxon>\n'
    return taxa_string


def generate_alignment(acc, seq):
    align_string = ""
    align_string += '\t\t<sequence>\n'
    align_string += '\t\t\t<taxon idref="%s"/>\n' % acc
    align_string += '\t\t\t%s\n' % seq
    align_string += '\t\t</sequence>\n'
    return align_string


s1 = generate_taxon('EPI_ISL_984690', '2020.8469', 'forwards', 'years', Subtype='A / H5N8', host2='Other', Subregion='North_Europe', Latitude=55.33, Longitude=9.09, loc='55.33 9.09')
s2 = generate_taxon('EPI_ISL_444455', '2020.8469', 'forwards', 'years', Subtype='A / H5N8', host2='Other', Subregion='North_Europe', Latitude=55.33, Longitude=9.09, loc='55.33 9.09')


des = """<?xml version="1.0" standalone="yes"?>

<beast version="1.10.4">
"""

print(des)
print('\t<taxa id="taxa">')

print(s1)
print(s2)

print('\t</taxa>')
print('\t<alignment id="alignment" dataType="nucleotide">')


a1 = generate_alignment('EPI_ISL_984690', 'ATGGAAGACTTTGTGCGACAATGCTTCAATCCAATGATTGTCGAGCTTGCGGAAAAGGCAATGAAAGAATATGGGGAAGATCCTAAGATCGAAACAAACAAGTTTGCCGCAATATGTACACACTTAGAAGTCTGCTTCATGTATTCGGATTTCCATTTCATTGATGAACGAGGCGAATCAATAATCGTAGAATCTGGCGACCCGAATGCATTATTGAAGCACCGATTTGAGATAATTGAAGGGAGAGACCGAACAATGGCCTGGACAGTGGTGAATAGTATATGCAACACTACAGGAGTCGAAAAGCCCAAGTTCCTTCCTGATTTGTATGACTACAGAGAAAACCGATTCATTGAAATTGGAGTAACGCGAAGGGAAGTTCACATATACTATCTTGAAAAAGCCAACAAGATAAAATCAGAGAAAACACACATTCACATATTCTCATTCACTGGGGAGGAAATGGCCACCAAGGCAGATTACACTCTTGATGAAGAGAGCAGAGCAAGAATAAAAACCAGGCTATTCACTATAAGACAAGAGATGGCCAATAGGGGTCTATGGGATTCCTTTCGTCAGTCCGAGAGAGGCGAAGAGACAATTGAAGAAAGATTTGAAATCACAGGAACCATGCGCAGGCTTGCTGACCAAAGTCTCCCACCGAACTTCGCCAGCCTTGAAAACTTTAGAGCCTATGTGGATGGATTTGAACCGAACGGCTGCATTGAGGGCAAGCTTTCTCAAATGTCAAAAGAAGTGAATGCCAGAATTGAGCCATTTCTGAAGACAACACCACGCCCTCTCAGATTACCTGATGGACCTCCCTGTTCCCAGCGGTCGAAATTCTTGCTAATGGATGCCCTTAAATTGAGCATTGAAGACCCGAGCCATGAGGGGGAGGGTATACCGCTGTACGATGCAACCAAATGCATGAAGACATTTTTTGGCTGGAAGGAGCCCAACATCGTGAAACCACATGAAAAGGGCATAAACCCTAATTACCTCCTGGCTTGGAAGCAGGTGCTATCAGAACTCCAAGATATTGAAAACGAGGA')
a2 = generate_alignment('EPI_ISL_444455', 'ATGGAAGACTTTGTGCGACAATGCTTCAATCCAATGATTGTCGAGCTTGCGGAAAAGGCAATGAAAGAATATGGGGAAGATCCTAAGATCGAAACAAACAAGTTTGCCGCAATATGTACACACTTAGAAGTCTGCTTCATGTATTCGGATTTCCATTTCATTGATGAACGAGGCGAATCAATAATCGTAGAATCTGGCGACCCGAATGCATTATTGAAGCACCGATTTGAGATAATTGAAGGGAGAGACCGAACAATGGCCTGGACAGTGGTGAATAGTATATGCAACACTACAGGAGTCGAAAAGCCCAAGTTCCTTCCTGATTTGTATGACTACAGAGAAAACCGATTCATTGAAATTGGAGTAACGCGAAGGGAAGTTCACATATACTATCTTGAAAAAGCCAACAAGATAAAATCAGAGAAAACACACATTCACATATTCTCATTCACTGGGGAGGAAATGGCCACCAAGGCAGATTACACTCTTGATGAAGAGAGCAGAGCAAGAATAAAAACCAGGCTATTCACTATAAGACAAGAGATGGCCAATAGGGGTCTATGGGATTCCTTTCGTCAGTCCGAGAGAGGCGAAGAGACAATTGAAGAAAGATTTGAAATCACAGGAACCATGCGCAGGCTTGCTGACCAAAGTCTCCCACCGAACTTCGCCAGCCTTGAAAACTTTAGAGCCTATGTGGATGGATTTGAACCGAACGGCTGCATTGAGGGCAAGCTTTCTCAAATGTCAAAAGAAGTGAATGCCAGAATTGAGCCATTTCTGAAGACAACACCACGCCCTCTCAGATTACCTGATGGACCTCCCTGTTCCCAGCGGTCGAAATTCTTGCTAATGGATGCCCTTAAATTGAGCATTGAAGACCCGAGCCATGAGGGGGAGGGTATACCGCTGTACGATGCAACCAAATGCATGAAGACATTTTTTGGCTGGAAGGAGCCCAACATCGTGAAACCACATGAAAAGGGCATAAACCCTAATTACCTCCTGGCTTGGAAGCAGGTGCTATCAGAACTCCAAGATATTGAAAACGAGGA')

print(a1)
print(a2)

print('\t</alignment>')




class Seq:
    def __init__(self, acc=None, seq=None, date=None, host=None, region=None, latitude=None, longitude=None, **kwargs):
        self.acc = acc
        self.seq = seq
        self.date = date
        self.host = host
        self.region = region
        self.latitude = latitude
        self.longitude = longitude
        for name, value in kwargs.items():
            self.name = value

    def __repr__(self):
        return ">%s\n%s\n" % (self.acc, self.seq)


class partitionos:
    def __init__(self, partion=None):
        if partion not in [None, 2, 3]:
            raise ValueError("Partion must be 2 or 3")
        self.partion = partion

    def __str__(self):
        if self.partion == 2:
            return self._patterns2()
        elif self.partion == 2:
            pass
        elif not self.partion:
            return self._patterns1()

    def _patterns2(self):
        partition_string = ''
        partition_string += '\t<mergePatterns id="CP1+2.patterns">\n'
        partition_string += '\t\t<patterns from="1" every="3" strip="false">\n'
        partition_string += '\t\t\t<alignment idref="alignment"/>\n'
        partition_string += '\t\t</patterns>\n'
        partition_string += '\t\t<patterns from="2" every="3" strip="false">\n'
        partition_string += '\t\t\t<alignment idref="alignment"/>\n'
        partition_string += '\t\t</partterns>\n'
        partition_string += '\t</mergePatterns>\n\n'
        partition_string += '\t<partterns id="CP3.patterns" from="3" every="3" strip="false">\n'
        partition_string += '\t\t<alignment idref="alignment"/>\n'
        partition_string += '\t</patterns>\n\n'
        return partition_string

    def _patterns1(self):
        text = '''
               \t<patterns id="patterns" from="1" strip="false">
               \t\t<alignment idref="alignment"/>
               \t</patterns>
               '''
        return text


# par = partitionos()
# print(par)

par = partitionos(partion=2)
print(par)


class taxaSet:
    def __init__(self, name, values):
        self.name = name
        self.values = values

    def __str__(self):
        taxa_string = ''
        taxa_string += '\t<taxa id="%s">\n' % self.name
        for value in self.values:
            taxa_string += '\t\t<taxon idref="%s"/>\n' % value
        taxa_string += '\t</taxa>\n'
        return taxa_string

taxaset = taxaSet('in_group', ['EPI1', 'EPI2', 'EPI3'])
# print(taxaset)



class discreteTraitModel:
    def __init__(self, name=None, values=None):
        self.name = name
        self.values = values

    def generalDataType(self):
        trait_string = ''
        trait_string += '\t<generalDataType id="%s.dataType">\n' % self.name
        for value in self.values:
            trait_string += '\t\t<state code="%s"/>\n' % value
        trait_string += '\t</generalDataType>\n'
        return trait_string

    def attributePatterns(self):
        attr_string = ''
        attr_string += '\t<attributePatterns id="%s.pattern" attribute="%s">\n' % (self.name, self.name)
        attr_string += '\t\t<taxa idref="taxa"/>\n'
        attr_string += '\t\t<generalDataType idref="%s.dataType"/>\n' % self.name
        attr_string += '\t</attributePatterns>\n'
        return attr_string


dt1 = discreteTraitModel('host2', ['Waterfowl', 'Wild-gal', 'Dom-ans', 'Dom-gal', 'Seabird', 'Human', 'Other'])
print(dt1.generalDataType())
print(dt1.attributePatterns())



class treeModel:
    def __init__(self, tree_prior=None, starting_tree=None):
        self.tree_prior = tree_prior
        self.starting_tree = starting_tree

    def startingTree(self):
        starting_string = ''
        starting_string += '\t<constantSize id="constant" units="years">\n'
        starting_string += '\t\t<populationSize>\n'
        starting_string += '\t\t\t<parameter id="constant.popSize" value="1.0" lower="0.0"/>\n'
        starting_string += '\t\t</populationSize>\n'
        starting_string += '\t</constantSize>\n\n'
        # starting_string += '\t<coalescentSimulator id="startingTree">\n'
        # starting_string += '\t\t<taxa idref="taxa"/>\n'
        # starting_string += '\t\t<constantSize idref="constant"/>\n'
        # starting_string += '\t</coalescentSimulator>\n\n'
        return starting_string

    def treePrior(self, fileName):
        empiricalDis_string = ''
        empiricalDis_string += '\t<empiricalTreeDistributionModel id="treeModel" fileName="%s">\n' % fileName
        empiricalDis_string += '\t\t<taxa idref="taxa"/>\n'
        empiricalDis_string += '\t</empiricalTreeDistributionModel>\n\n'
        return empiricalDis_string

    def treeStatistic(self):
        statistic_string = ''
        statistic_string += '\t<statistic id="treeModel.currentTree" name="Current Tree">\n'
        statistic_string += '\t\t<empiricalTreeDistributionModel idref="treeModel"/>\n'
        statistic_string += '\t</statistic>\n\n'

        statistic_string += '\t<treeLengthStatistic id="treeLength">\n'
        statistic_string += '\t\t<treeModel idref="treeModel"/>\n'
        statistic_string += '\t</treeLengthStatistic>\n\n'
        return statistic_string

    def treePriorLikelihood(self):
        treePriorLH_string = ''
        treePriorLH_string += '\t<coalescentLikelihood id="coalescent">\n'
        treePriorLH_string += '\t\t<model>\n'
        treePriorLH_string += '\t\t\t<constantSize idref="constant"/>\n'
        treePriorLH_string += '\t\t</model>\n'
        treePriorLH_string += '\t\t<populationTree>\n'
        treePriorLH_string += '\t\t\t<treeModel idref="treeModel"/>\n'
        treePriorLH_string += '\t\t</populationTree>\n'
        treePriorLH_string += '\t</coalescentLikelihood>\n\n'
        return treePriorLH_string


tm = treeModel()
# print(tm.startingTree())
# print(tm.treePrior('resample3_ext_PA.align.trees'))
# print(tm.treeStatistic())
# print(tm.treePriorLikelihood())


class strictClockModel:
    def __init__(self, name):
        self.name = name

    def strictClockBranchRates(self):
        branchRates_string = ''
        branchRates_string += '\t<strictClockBranchRates id="%s.branchRates">\n' % self.name
        branchRates_string += '\t\t<rate>\n'
        branchRates_string += '\t\t\t<parameter id="%s.clock.rate" value="1.0" lower="0.0"/>\n' % self.name
        branchRates_string += '\t\t</rate>\n'
        branchRates_string += '\t</strictClockBranchRates>\n\n'
        return branchRates_string

    def rateStatistic(self):
        rateStatistic_string = ''
        rateStatistic_string += '\t<rateStatistic id="%s.meanRate" name="%s.meanRate" mode="mean" internal="true" external="true">\n' % (self.name, self.name)
        rateStatistic_string += '\t\t<treeModel idref="treeModel"/>\n'
        rateStatistic_string += '\t\t<strictClockBranchRates idref="%s.branchRates"/>\n' % self.name
        rateStatistic_string += '\t</rateStatistic>\n\n'
        return rateStatistic_string


# sc = strictClockModel('Subtype')
# print(sc.strictClockBranchRates())
# print(sc.rateStatistic())

# sc = strictClockModel('host2')
# print(sc.strictClockBranchRates())
# print(sc.rateStatistic())


class multivariateDiffusionModel():
    def __init__(self, name=None):
        self.name = name

    def multivariateDiffusionModel(self):
        multivariateDM_stirng = ''
        multivariateDM_stirng += '\t<multivariateDiffusionModel id="%s.diffusionModel">\n' % self.name
        multivariateDM_stirng += '\t\t<precisionMatrix>\n'
        multivariateDM_stirng += '\t\t\t<matrixParameter id="%s.precision">\n' % self.name
        multivariateDM_stirng += '\t\t\t\t<parameter id="%s.precision.col1" value="0.05 0.002"/>\n' % self.name
        multivariateDM_stirng += '\t\t\t\t<parameter id="%s.precision.col2" value="0.002 0.05"/>\n' % self.name
        multivariateDM_stirng += '\t\t\t</matrixParameter>\n'
        multivariateDM_stirng += '\t\t</precisionMatrix>\n'
        multivariateDM_stirng += '\t</multivariateDiffusionModel>\n\n'
        return multivariateDM_stirng

    def wishartPrior(self):
        prior_string = ''
        prior_string += '\t<multivariateWishartPrior id="%s.precisionPrior" df="2">\n' % self.name
        prior_string += '\t\t<scaleMatrix>\n'
        prior_string += '\t\t\t<matrixParameter>\n'
        prior_string += '\t\t\t\t<parameter value="1.0 0.0"/>\n'
        prior_string += '\t\t\t\t<parameter value="0.0 1.0"/>\n'
        prior_string += '\t\t\t</matrixParameter>\n'
        prior_string += '\t\t</scaleMatrix>\n'
        prior_string += '\t\t<data>\n'
        prior_string += '\t\t\t<parameter idref="%s.precision"/>\n' % self.name
        prior_string += '\t\t</data>\n'
        prior_string += '\t</multivariateWishartPrior>\n\n'
        return prior_string

    def branchRates(self):
        branchRates_string = ''
        branchRates_string += '\t<arbitraryBranchRates id="%s.diffusion.branchRates">\n' % self.name
        branchRates_string += '\t\t<treeModel idref="treeModel"/>\n'
        branchRates_string += '\t\t<rates>\n'
        branchRates_string += '\t\t\t<parameter id="%s.diffusion.rates" lower="0.0"/>\n' % self.name
        branchRates_string += '\t\t</rates>\n'
        branchRates_string += '\t</arbitraryBranchRates>\n'

        branchRates_string += '\t<distributionLikelihood id="%s.diffusion.prior">\n' % self.name
        branchRates_string += '\t\t<data>\n'
        branchRates_string += '\t\t\t<parameter idref="%s.diffusion.rates"/>\n' % self.name
        branchRates_string += '\t\t</data>\n'
        branchRates_string += '\t\t<distribution>\n'
        branchRates_string += '\t\t\t<onePGammaDistributionModel>\n'
        branchRates_string += '\t\t\t\t<shape>\n'
        branchRates_string += '\t\t\t\t\t<parameter value="0.5"/>\n'
        branchRates_string += '\t\t\t\t</shape>\n'
        branchRates_string += '\t\t\t</onePGammaDistributionModel>\n'
        branchRates_string += '\t\t</distribution>\n'
        branchRates_string += '\t</distributionLikelihood>\n\n'

        branchRates_string += '\t<multivariateTraitLikelihood id="%s.traitLikelihood" traitName="%s" useTreeLength="true" scaleByTime="true" reportAsMultivariate="true" reciprocalRates="true" integrateInternalTraits="true">\n'  % (self.name, self.name)
        branchRates_string += '\t\t<multivariateDiffusionModel idref="%s.diffusionModel"/>\n'  % self.name
        branchRates_string += '\t\t<treeModel idref="treeModel"/>\n'
        branchRates_string += '\t\t<traitParameter>\n'
        branchRates_string += '\t\t\t<parameter id="leaf.%s"/>\n' % self.name
        branchRates_string += '\t\t</traitParameter>\n'
        branchRates_string += '\t\t<jitter window="0.01 0.01" duplicatesOnly="true">\n'
        branchRates_string += '\t\t\t<parameter id="leaf.%s"/>\n' % self.name
        branchRates_string += '\t\t</jitter>\n'
        branchRates_string += '\t\t<conjugateRootPrior>\n'
        branchRates_string += '\t\t\t<meanParameter>\n'
        branchRates_string += '\t\t\t\t<parameter value="0.0 0.0">\n'
        branchRates_string += '\t\t\t</meanParameter>\n'
        branchRates_string += '\t\t\t<priorSampleSize>\n'
        branchRates_string += '\t\t\t\t<parameter value="0.000001"/>\n'
        branchRates_string += '\t\t\t</priorSampleSize>\n'
        branchRates_string += '\t\t</conjugateRootPrior>\n'
        branchRates_string += '\t\t<arbitraryBranchRates idref="%s.diffusion.branchRates"/>\n' % self.name
        branchRates_string += '\t</multivariateTraitLikehood>\n\n'

        branchRates_string += '\t<correlation id="%s.correlation" dimension1="1" dimension2="2">\n' % self.name
        branchRates_string += '\t\t<matrixParameter idref="%s.precision"/>\n' % self.name
        branchRates_string += '\t</correlation>\n'

        branchRates_string += '\t<matrixInverse id="%s.varCovar">\n' % self.name
        branchRates_string += '\t\t<matrixParameter idref="%s.precision"/>\n' % self.name
        branchRates_string += '\t</matrixInverse>\n'

        branchRates_string += '\t<continuousDiffusionStatistic id="%s.diffusionRate" greatCircleDistance="true">\n' % self.name
        branchRates_string += '\t\t<multivariateTraitLikelihood idref="%s.traitLikelihood"/>\n' % self.name
        branchRates_string += '\t</continuousDiffusionStatistic>\n\n'
        return branchRates_string


# mll = multivariateDiffusionModel('loc')
# print(mll.multivariateDiffusionModel())
# print(mll.wishartPrior())
# print(mll.branchRates())



class ctmcModel:
    def __init__(self, name):
        self.name = name

    def generalSubstitutionModel(self):
        gtm_string = ''
        gtm_string += '\t<generalSubstitutionModel id="%s.model">\n' % self.name
        gtm_string += '\t\t<generalDataType idref="%s.dataType"/>\n' % self.name
        gtm_string += '\t\t<frequencies>\n'
        gtm_string += '\t\t\t<frequencyModel id="%s.frequencyModel" normalize="true">\n' % self.name
        gtm_string += '\t\t\t\t<generalDataType idref="Subtype.dataType"/>\n'
        gtm_string += '\t\t\t\t<frequencies>\n'
        gtm_string += '\t\t\t\t\t<parameter id="%s.frequencies" dimension="10"/>\n' % self.name
        gtm_string += '\t\t\t\t</frequencies>\n'
        gtm_string += '\t\t\t</frequencyModel>\n'
        gtm_string += '\t\t</frequencies>\n'
        gtm_string += '\t\t<rates>\n'
        gtm_string += '\t\t\t<parameter id="%s.rates" dimension="45" value="1.0" lower="0.0"/>\n' % self.name
        gtm_string += '\t\t</rates>\n'
        gtm_string += '\t\t<rateIndicator>\n'
        gtm_string += '\t\t\t<parameter id="%s.indicators" dimension="45" value="1.0"/>\n' % self.name
        gtm_string += '\t\t</rateIndicator>\n'
        gtm_string += '\t</generalSubstitutionModel>\n'

        gtm_string += '\t<sumStatistic id="%s.nonZeroRates" elementwise="true">\n' % self.name
        gtm_string += '\t\t<parameter idref="%s.indicators"/>\n' % self.name
        gtm_string += '\t</sumStatistic>\n'
        gtm_string += '\t<productStatistic id="%s.actualRates" elementwise="false">\n' % self.name
        gtm_string += '\t\t<parameter idref="%s.indicators"/>\n' % self.name
        gtm_string += '\t\t<parameter idref="%s.rates"/>\n' % self.name
        gtm_string += '\t</productStatistic>\n'
        gtm_string += '\t<siteModel id="%s.siteModel">\n' % self.name
        gtm_string += '\t\t<substitutionModel>\n'
        gtm_string += '\t\t\t<generalSubstitutionModel idref="%s.model"/>\n' % self.name
        gtm_string += '\t\t</substitutionModel>\n'
        gtm_string += '\t</siteModel>\n\n'
        return gtm_string

    def ancestralTreeLikelihood(self):
        ancestralllh_string = ''
        ancestralllh_string += '\t<ancestralTreeLikehood id="%s.treeLikehood" stateTagName="%s.states"  useUniformization="true" saveCompleteHistory="false" logCompleteHistory="false">\n' % (self.name, self.name)
        ancestralllh_string += '\t\t<attributePatterns idref="%s.pattern"/>\n' % self.name
        ancestralllh_string += '\t\t<treeModel idref="treeModel"/>\n'
        ancestralllh_string += '\t\t<siteModel idref="%s.siteModel"/>\n' % self.name
        ancestralllh_string += '\t\t<generalSubstitutionModel idref="%s.model"/>\n' % self.name
        ancestralllh_string += '\t\t<generalSunstitutionModel idref="%s.model"/>\n' % self.name
        ancestralllh_string += '\t\t<strictClockBranchRates idref="%s.branchRates"/>\n' % self.name
        ancestralllh_string += '\t</ancestralTreeLikelihood>\n\n'
        return ancestralllh_string


dfM1 = ctmcModel('Subtype')
dfM2 = ctmcModel('host2')
dfM3 = ctmcModel('Subregion')

# print(dfM1.generalSubstitutionModel())
# print(dfM2.generalSubstitutionModel())
# print(dfM3.generalSubstitutionModel())

# print(dfM1.ancestralTreeLikelihood())
# print(dfM2.ancestralTreeLikelihood())
# print(dfM3.ancestralTreeLikelihood())


class operators:
    def __init__(self):
        pass

    def _scaleOperator(self, name, sF=0.75, weight=3, scaleAllIndependently=None):
        sO = ''
        if scaleAllIndependently:
            sO += '\t\t<scaleOperator scaleFactor="%s" weight="%s" scaleAllIndependently="true">\n' % (sF, weight)
        else:
            sO += '\t\t<scaleOperator scaleFactor="%s" weight="%s">\n' % (sF, weight)
        sO += '\t\t\t<parameter idref="%s"/>\n' % name
        sO += '\t\t</scaleOperator>\n\n'
        return sO

    def _bitFlipOperator(self, name, weight):
        bO = ''
        bO += '\t\t<bitFlipOperator weight="%s">\n' % weight
        bO += '\t\t\t<parameter idref="%s"/>\n' % name
        bO += '\t\t</bitFilpOperator>\n\n'
        return bO

    def _empiricalTreeDistributionOperator(self, weight):
        eO = ''
        eO += '\t\t<empiricalTreeDistributionOperator weight="%s">\n' % weight
        eO += '\t\t\t<empiricalTreeDistributionModel idref="treeModel"/>\n'
        eO += '\t\t</empiricalTreeDistributionOperator>\n\n'
        return eO

    def _precisionGibbsOper(self, name, weight):
        gO = ''
        gO += '\t\t<precisionGibbsOper weight="%s">\n' % weight
        gO += '\t\t\t<multivariateTraitLikehood idref="%s.traitLikelihood"/>\n' % name
        gO += '\t\t\t<miltivariateWishartPrior idref="%s.precisionPrior"/>\n' % name
        gO += '\t\t</precisionGibbsOper>\n\n'
        return gO


    def main(self):
        string = ''
        string += '\t<operators id="operators" optimizationSchedule="log">\n\n'
        string += self._empiricalTreeDistributionOperator(weight=3)
        string += self._scaleOperator(name='Subtype.clock.rate', sF=0.75, weight=3)
        string += self._scaleOperator(name='host2.clock.rate', sF=0.75, weight=3)
        string += self._scaleOperator(name='Subregion.clock.rate', sF=0.75, weight=3)
        string += self._scaleOperator(name='loc.diffusion.rates', sF=0.75, weight=30)
        string += self._scaleOperator(name='Subtype.rates', sF=0.75, weight=15, scaleAllIndependently="true")
        string += self._bitFlipOperator(name='Subtype.indicators', weight=7)
        string += self._scaleOperator(name='host2.rates', sF=0.75, weight=15, scaleAllIndependently="true")
        string += self._bitFlipOperator(name='host2.indicators', weight=7)
        string += self._scaleOperator(name='Subregion.rates', sF=0.75, weight=15, scaleAllIndependently="true")
        string += self._bitFlipOperator(name='Subregion.indicators', weight=7)

        string += self._precisionGibbsOper(name='loc', weight=2)

        string += '\t</operators>'
        return string


opp = operators()
# print(opp.main())


class mcmc():
    def __init__(self, chainLength=None, logFileF=None, logScreenF=None):
        self.chainLength = chainLength
        self.logFileF = logFileF
        self.logScreenF = logScreenF


    def ctmcScalePrior(self, name):
        string = ''
        string += '\t\t\t\t<ctmcScalePrior>\n'
        string += '\t\t\t\t\t<ctmcScale>\n'
        string += '\t\t\t\t\t\t<parameter idref="%s.clock.rate"/>\n' % name
        string += '\t\t\t\t\t</ctmcScale>\n'
        string += '\t\t\t\t\t<treeModel idref="treeModel"/>\n'
        string += '\t\t\t\t</ctmcScalePrior>\n'
        return string

    def poissonPrior(self, name, mean, offset):
        string = ''
        string += f'\t\t\t\t<poissonPrior mean="{mean}" offset="{offset}">\n'
        string += '\t\t\t\t\t<statistic idref="{name}.nonZeroRates/>"\n'
        string += '\t\t\t\t</poissonPrior>\n'
        return string

    def uniformPrior(self, name, lower, upper):
        string = ''
        string += f'\t\t\t\t<uniformPrior lower="{lower}" upper="{upper}">\n'
        string += f'\t\t\t\t\t<parameter idref="{name}.frequencies"/>\n'
        string += '\t\t\t\t</uniformPrior>\n'
        return string

    def cachedPrior(self, name, shape, scale, offset):
        string = ''
        string += '\t\t\t\t<cachedPrior>\n'
        string += f'\t\t\t\t\t<gammaPrior shape="{shape}" scale="{scale}" offset="{offset}">\n'
        string += f'\t\t\t\t\t\t<parameter idref="{name}.rates"/>\n'
        string += '\t\t\t\t\t</gammaPrior>\n'
        string += f'\t\t\t\t\t<parameter idref="{name}.rates"/>\n'
        string += '\t\t\t\t</cachedPrior>\n'
        return string



    def main(self):
        string = ''
        string += '\t<mcmc id="mcmc" chainLength="%s" autoOptimize="true">\n' % self.chainLength
        string += '\t\t<joint id="joint">\n'
        string += '\t\t\t<prior id="prior">\n'
        string += self.ctmcScalePrior(name='Subtype')
        string += self.ctmcScalePrior(name='host2')
        string += self.ctmcScalePrior(name='Subregion')
        string += self.poissonPrior(name='Subtype', mean=0.6931471805599453, offset=9.0)
        string += self.uniformPrior(name='Subtype', lower=0.0, upper=1.0)
        string += self.cachedPrior(name='Subtype', shape=1.0, scale=1.0, offset=0.0)
        string += self.poissonPrior(name='host2', mean=0.6931471805599453, offset=4.0)
        string += self.uniformPrior(name='host2', lower=0.0, upper=1.0)
        string += self.cachedPrior(name='host2', shape=1.0, scale=1.0, offset=0.0)
        string += self.poissonPrior(name='Subregion', mean=0.6931471805599453, offset=11.0)
        string += self.uniformPrior(name='Subregion', lower=0.0, upper=1.0)
        string += self.cachedPrior(name='Subregion', shape=1.0, scale=1.0, offset=0.0)

        string += '\t\t\t\t<coalescentLikehood idref="coalescent"/>\n'
        string += f'\t\t\t\t<strictClockBranchRates idref="{"Subtype"}.branchRates"/>\n'
        string += f'\t\t\t\t<strictClockBranchRates idref="{"host2"}.branchRates"/>\n'
        string += f'\t\t\t\t<strictClockBranchRates idref="{"Subregion"}.branchRates"/>\n'

        string += f'\t\t\t\t<distributionLikehood idref="{"loc"}.diffusion.prior"/>\n'
        string += f'\t\t\t\t<multivariateWishartPrior idref="{"loc"}.precisionPrior"/>\n'

        string += f'\t\t\t\t<generalSubstitutionModel idref="{"Subtype"}.model"/>\n'
        string += f'\t\t\t\t<generalSubstitutionModel idref="{"host2"}.model"/>\n'
        string += f'\t\t\t\t<generalSubstitutionModel idref="{"Subregion"}.model"/>\n'
        string += '\t\t\t</prior>\n'

        string += '\t\t\t<likelihood id="likelihood">\n'
        string += '\t\t\t\t<multivariateTraitLikelihood idref="loc.traitLikelihood"/>\n'
        string += '\t\t\t\t<ancestralTreeLikelihood idref="Subtype.treeLikelihood"/>\n'
        string += '\t\t\t\t<ancestralTreeLikelihood idref="host2.treeLikelihood"/>\n'
        string += '\t\t\t\t<ancestralTreeLikelihood idref="Subregion.treeLikelihood"/>\n'
        string += '\t\t\t</likelihood>\n'
        string += '\t\t</joint>\n'
        string += '\t\t<operators idref="operators"/>\n'

        string += '\t</mcmc>\n'
        return string




mcmc = mcmc()
# print(mcmc.main())
