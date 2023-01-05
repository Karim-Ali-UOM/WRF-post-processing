from src.utilities import filename, message, newline, prepareloc, turbCoordCode, splitline, turbTypeCode,\
    extractDynamicString, wfnamesCode
import glob

class globe:
    outdir = str
    ntLabel = str
    dirs = []
    labels = []
    bases = []
    domainIndex= []
    processedTurbLoc = []
    turbLatDict = {}
    turbLonDict = {}
    turbWindFarmIndexDict = {}
    turbTypesDict = {}
    windFarmTurbTypeDict = {}
    turbTypeHHDict = {}
    turbTypeDiameterDict = {}
    turbTypePower = {}
    windfarmsNamesDict = {}
    windfarmsRatedPowerDict = {}
    windFarmNumTurbinesDict = {}
    numWindFarms = 0
    windfarmIndices = []

    def __init__(self,dirs,labels,bases, ntLabel,outdir,domainIndex, processedTurbLoc):
        self.dirs = dirs.copy()
        self.labels = labels.copy()
        self.ntLabel = ntLabel
        self.bases = bases.copy()
        self.outdir = outdir
        self.domainIndex = domainIndex.copy()
        self.processedTurbLoc = processedTurbLoc.copy()

        msglvl = 0
        newline()
        self.initOutFolder(msglvl)
        self.readTurbineCoords(msglvl)
        self.readTurbineTypes(msglvl)
        self.readFarmsNames(msglvl)
        self.initiateRatedPowers(msglvl)
        self.numWindFarms = len(self.windfarmsNamesDict)
        self.windfarmIndices = list(self.windfarmsNamesDict.keys())
        message("Finished creating the class: globe", msglvl)

    def initOutFolder(self, msglvl):
        message("Creating output folder", msglvl)
        prepareloc(self.outdir)

    def readTurbineCoords(self, msglvl):
        file = self.processedTurbLoc[self.domainIndex.index(turbCoordCode())]
        message("Reading turbine coordinates file (" + filename(file) + ")", msglvl)
        datafile = open(file,"r")
        counter = 0
        while True:
            try: line = splitline(datafile.readline()," ")
            except: break
            if len(line) == 0 or line[0] in ["", "\n"]: break
            counter += 1
            self.turbLatDict[counter] = float(line[0])
            self.turbLonDict[counter] = float(line[1])
            self.turbTypesDict[counter] = int(line[2])
            self.turbWindFarmIndexDict[counter] = int(line[3])
            if not int(line[3]) in self.windFarmTurbTypeDict.keys():
                self.windFarmTurbTypeDict[int(line[3])] = int(line[2])
            if int(line[3]) in self.windFarmNumTurbinesDict.keys():
                self.windFarmNumTurbinesDict[int(line[3])] += 1
            else:
                self.windFarmNumTurbinesDict[int(line[3])] = 1
        datafile.close()

    def readTurbineTypes(self, msglvl):
        dynamicFilename = self.processedTurbLoc[self.domainIndex.index(turbTypeCode())]
        files =  glob.glob(dynamicFilename)
        for file in files:
            message("Reading turbine type file (" + filename(file) + ")", msglvl)
            turbType = int(extractDynamicString(file,dynamicFilename))
            data = open(file,"r")
            data.readline() # number of table entries
            splits = splitline(data.readline()," ")
            self.turbTypePower[turbType] = float(splits[3]) * 1000000
            self.turbTypeDiameterDict[turbType] = float(splits[1])
            self.turbTypeHHDict[turbType] = float(splits[0])
            data.close()

    def readFarmsNames(self, msglvl):
        file = self.processedTurbLoc[self.domainIndex.index(wfnamesCode())]
        message("Reading wind farms names from " + filename(file), msglvl)
        datafile = open(file ,"r")
        counter = 0
        while True:
            try: line = splitline(datafile.readline(),"\n")
            except: break
            if len(line) == 0 or line[0] in ["", "\n"]: break
            counter += 1
            self.windfarmsNamesDict[counter] = line[0]
        datafile.close()

    def initiateRatedPowers(self, msglvl):
        message("Initiating wind farms rated power dictionary", msglvl)
        self.windfarmsRatedPowerDict = {"Gode Wind":582, "Veja Mate":402, "Bard Offshore 1":400, "Riffgat": 108,
            "Borkum Riffgrund 1":312, "Meerwind":288, "Nordsee Ost":295.2, "Amrumbank West":288, "Alpha Ventus":60,
            "Trianel Windpark Borkum I":200, "Global Tech I":400, "Nordsee One":332, "Gemini":600}

    def windFarmTurbineType(self, windFarmIndex):
        return self.windFarmTurbTypeDict[windFarmIndex]

    def turbineTypePower(self, turbType):
        return self.turbTypePower[turbType]

    def turbineIndexToWindFarmIndex(self,turb):
        return self.turbWindFarmIndexDict[turb]

    def farmRatedPower(self, wf):
        return self.windfarmsRatedPowerDict[wf]

    def windFarmName(self, index):
        return self.windfarmsNamesDict[index]

    def numTurbsPerfarm(self, wf):
        return self.windFarmNumTurbinesDict[wf]

    def turbineTypeDiameter(self, turbType):
        return self.turbTypeDiameterDict[turbType]
    
    def turbTypeHH(self, turbType):
        return self.turbTypeHHDict[turbType]