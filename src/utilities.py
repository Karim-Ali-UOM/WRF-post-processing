from os.path import isdir, isfile, join, dirname
from os import mkdir, listdir
from netCDF4 import Dataset	
from wrf import ll_to_xy, getvar, latlon_coords, tk, to_np, get_cartopy, destagger
import csv
import numpy as np
from sys import exit
import matplotlib.pyplot as plt
import matplotlib
from math import pi

def fataError(message):
    exit(message)

def splitline(l, char):
    return [x.strip() for x in l.split(char) if x != ""]

def essentialKeywords():
    return [dirsKeyword(), outkeyword(), turbLocKeyword(), ordersKeyword()]

def validOrderIdentifiers():
    return [activeIdentifier(), typeIdentifier(), wfIndexIdentifier(), intervalsIdentifier()]

def subpacketIdentifiers():
    return [assistFilesIdentifier(), intervalsIdentifier(), pairsIdentifier(), durIdentifier(), segmentsIdentifier()]

def mainPacketKeyword():
    return "WRF-postProcessing"

def dirsKeyword():
    return "Directories"

def refCaseKeyword():
    return "Reference"

def outkeyword():
    return "Output"

def ordersKeyword():
    return "Orders"

def durIdentifier():
    return "Duration"

def segmentsIdentifier():
    return "Segments"

def locIdentifier():
    return "Location"

def turbLocKeyword():
    return "TurbineData"

def startIdentifier():
    return "Start"

def endIdentifier():
    return "End"

def activeIdentifier():
    return "Active"

def typeIdentifier():
    return "Type"

def turbineTypesTablesIdentifier():
    return "TurbineTypes"

def turbCoordIdentifier():
    return "Coordinates"

def windFarmsnamesIdentifier():
    return "WindFarmNames"

def intervalsIdentifier():
    return "Intervals"

def pairsIdentifier():
    return "Pairs"

def levelsIdentifier():
    return "Levels"

def wakeContourName():
    return "WakeContourPlot"

def verticalProfilesName():
    return "VerticalProfilesPlot"

def mtendName():
    return "MTEND"

def utendname():
    return "UTEND"

def vtendname():
    return "VTEND"

def heightIdentifier():
    return "Height"

def turbCoordCode():
    return -1

def turbTypeCode():
    return -2

def wfnamesCode():
    return -3

def wfIndexIdentifier():
    return "WindFarmsIndices"

def wfcapacityName():
    return "windFarmsCapacityFactors"

def qbudgetName():
    return "qVerticalAverage"

def wfstatsName():
    return "windFarmStats"

def inputFileName():
    return "wrfPostProcessingDict"

def calcTurbLocIdentifier():
    return "calculate"

def readTurbLocIdentifier():
    return "readLocs"

def turbinePowerFieldname():
    return "turbinepower".upper()

def powerfieldName():
    return "power".upper()

def q2BudgetFields():
    return [qvTransportName(), qsName(), qbName(), qdName(), qtName()]

def qvTransportName():
    return "qwt".upper()

def qsName():
    return "qshear".upper()

def qbName():
    return "qbuoy".upper()

def qdName():
    return "qdiss".upper()

def qtName():
    return "qtend".upper()

def listFiles(wd):
    return [wd + "/" + f for f in listdir(wd) if isfile(join(wd, f))]

def prepareloc(folder):
    if (not isdir(folder)): mkdir(folder)

def select(files, base):
    selected = []
    for file in files:
        splitname = splitline(splitline(file,"/")[-1],"_")
        if splitname[0] + "_" + splitname[1] == base: selected.append(file)
    return selected

def checkField(file, field):
    return field.upper() in Dataset(file).variables.keys()

def extractDomainIndex(base):
    return int(splitline(base,"_d")[-1])

def loadGrid(file):
    z = getvar(Dataset(file), "z")
    lats, lons= latlon_coords(z)
    return to_np(lats), to_np(lons), to_np(z), get_cartopy(z)

def fillInIndexMap(map, i, j, farmIndex, nv):
    assigned = False
    for k in range(nv):
        if map[j,i,2*k] == farmIndex:
            map[j,i,2*k+1] += 1
            assigned = True
    if not assigned:
        for k in range(nv):
            if map[j,i,2*k] == 0:
                map[j,i,2*k] = farmIndex
                map[j,i,2*k+1] += 1
    return map

def createIndexMap(wrfout, caseindex, method, globe, n1, n2, nv):
    turbs = np.zeros([n1, n2, nv])
    nv2 = int(nv / 2)

    if method == readTurbLocIdentifier:
        turbdata = csv.reader(open(globe.processedTurbLoc[globe.domainIndex.index(extractDomainIndex(globe.bases[caseindex]))]))

        for turbine in turbdata:
            if turbine[1] == " Lat": continue
            i = int(turbine[3])-1
            j = int(turbine[4])-1
            fillInIndexMap(turbs, i, j, int(turbine[8]), nv2)

    elif method == calcTurbLocIdentifier():
        locsFile = open(globe.processedTurbLoc[globe.domainIndex.index(turbCoordCode())],"r")
        while True:
            try: 
                line = splitline(locsFile.readline()," ")
            except:
                break
            if len(line) == 1: break
            i, j = ll_to_xy(Dataset(wrfout),float(line[0]),float(line[1]))
            fillInIndexMap(turbs, i, j, int(line[2]), nv2)
        locsFile.close()
    return turbs

def first_last_nz(x):							
    nonzeros = [i for i, e in enumerate(x) if e != 0]		
    try:								
        return int(min(nonzeros)), int(max(nonzeros))		
    except:								
    	return 0, 1000

def smooth(x,y) :
    from scipy.interpolate import make_interp_spline

    fz, lz = first_last_nz(y)
    # Smooth profile
    xsub = x[fz:lz+1]
    ysub = y[fz:lz+1]
    xnew = np.linspace(min(xsub), max(xsub), 300) 					
    spl = make_interp_spline(np.array(xsub), np.array(ysub), k=2)	
    ynew = spl(xnew)												
    return np.concatenate((np.concatenate((x[0:fz], xnew)), x[lz+1:])), \
           np.concatenate((np.concatenate((y[0:fz], ynew)), y[lz+1:]))

def uname():
    return "ua"

def vname():
    return "va"

def speed10Name():
    return "SPEED10"

def speedName():
    return "Speed"

def tkename():
    return "TKE"

def qkeName():
    return "QKE"

def directionName():
    return "DIR"

def tempName():
    return "Temp"

def TperturbName():
    return "T"

def pressureName():
    return "P"

def basePressureName():
    return "PB"

def maskName():
    return "XLAND"

def sarPlotName():
    return "SARplot"

def assistFilesIdentifier():
    return "AssistingFiles"

def timeInstantIdentifier():
    return "Time"

def capsIdentifier():
    return "Limits"

def rowCol(i, ncols):
    row = i // ncols
    col = i - row * ncols
    return row, col

def maskField(sol, mask, val):
    for j in range(sol.shape[0]):
        for i in range(sol.shape[1]):
            if mask[j,i] != val: sol[j,i] = float("Nan")
    return sol

def nearest(array,val):
    return (np.abs(np.array(array)-val)).argmin()

def field_diff_percent(field, ref):
    diff = field.copy() * 0
    for j in range(field.shape[0]):
        for i in range(field.shape[1]):
            if ref[j,i] == 0 and field[j,i] == 0:
                diff[j,i] == 0
            elif ref[j,i] == 0:
                diff[j,i] == float("Nan")
            else: 
                diff[j,i] = (-field[j,i] + ref[j,i]) / ref[j,i] * 100.0
    return diff

def field_diff_percent2(field, ref):
    return 100*(1 - np.divide(np.array(field), np.array(ref)))

def fieldExtract(file, field, wrfInputFile = None):
    if field in zStagFields():
        return to_np(destagger(getvar(Dataset(file), field), 0))
    elif field == speed10Name():
        return np.sqrt(fieldExtract(file, "U10")**2 + fieldExtract(file, "V10")**2)
    elif field == speedName():
        return np.sqrt(fieldExtract(file, uname())**2 + fieldExtract(file, vname())**2)
    elif field == mtendName():
        return np.sqrt(fieldExtract(file, utendname())**2 +
                       fieldExtract(file, vtendname())**2)
    elif field == tempName():
        tp = fieldExtract(file, TperturbName())
        pp = fieldExtract(file, pressureName())
        pb = fieldExtract(file, basePressureName())
        return to_np(tk(pp+pb,tp+300.0))-273.15
    elif field == tkename():
        return fieldExtract(file, qkeName()) / 2.0
    elif field == directionName():
        u = fieldExtract(file, uname())
        v = fieldExtract(file, vname())
        cosalpha = Dataset(wrfInputFile).variables["COSALPHA"][0]
        sinalpha = Dataset(wrfInputFile).variables["SINALPHA"][0]
        ue = u.copy()
        ve = v.copy()
        for k in range(u.shape[0]):
            ue[k,:,:] = u[k,:,:] * cosalpha + v[k,:,:] * sinalpha
            ve[k,:,:] = v[k,:,:] * cosalpha - u[k,:,:] * sinalpha
        return np.arctan2(ve, ue) * 180.0/np.pi + 180
    else:
        return to_np(getvar(Dataset(file), field))

def standingAverage(data, window):
    # first entry in each row is assumed a time instant
    # time instants are assumed to be arranged
    c1 = 0
    c2 = -1
    timer = 0
    OverallAvgData = []
    for j in range(len(data)-1):
        c2 += 1
        timer += data[j+1][0] - data[j][0]
        if timer>= window:
            avgData = []
            timer = 0
            avgData.append(data[int(c1 + (c2-c1+1)/2)][0])
            for k in range(len(data[0])-1):
                avgData.append(np.mean(np.array([data[l][k+1] for l in np.arange(c1,c2+1,1)])))
            c1 = c2 + 1
            OverallAvgData.append(avgData)
    return OverallAvgData

def extractDirectory(file):
    return dirname(file)

def extractDynamicString(name, dynName):
    for char in dynName:
        if char == "*": return name[dynName.index(char)]

def markWindFarmCells(file):
    power = fieldExtract(file, powerfieldName())
    for j in range(power.shape[0]):
        for i in range(power.shape[1]):
            if power[j,i] > 0: power[j,i] = 1
    lats, lons, _, _ = loadGrid(file)
    return lats, lons, power

def createPanel(ncols, nrows, nplots, l1, l2):
    matplotlib.use('Agg')
    fig, axgs = plt.subplots(nrows=nrows,ncols=ncols, figsize=(l1,l2))
    for i in np.arange(nplots, nrows*ncols, 1): fig.delaxes(axgs[i])
    return fig, axgs

def createCrsPanel(ncols, nrows, nplots, l1, l2, crss):
    matplotlib.use('Agg')
    fig, axgs = plt.subplots(nrows=nrows,ncols=ncols, subplot_kw={'projection': crss}, figsize=(l1,l2))
    for i in np.arange(nplots, nrows*ncols, 1): fig.delaxes(axgs[i-1])
    return fig, axgs

def createFigure(l1, l2):
    matplotlib.use('Agg')
    fig = plt.figure(figsize=(l1, l2))
    return fig

def dxname():
    return "DX"

def extractDirectVariable(file, field):
    if field == dxname():
        return Dataset(file).DX
    else:
        return -1

def dronPathIdentifier():
    return "DroneFiles"

def windFarmNameLocation(name):
    if name == "Gode Wind": return [7.117, 54.02]
    elif name == "Veja Mate": return [5.762, 54.18]
    elif name == "Bard Offshore 1": return [6.058, 54.38]
    elif name == "Riffgat": return [6.562, 53.69]
    elif name == "Borkum Riffgrund 1": return [6.3, 53.779]
    elif name == "Meerwind": return [7.545, 54.285]
    elif name == "Nordsee Ost": return [7.18, 54.42]
    elif name == "Amrumbank West": return [7.26, 54.57]
    elif name == "Alpha Ventus": return [6.618, 54.04]
    elif name ==  "Trianel Windpark Borkum I": return [6.2, 54.11]
    elif name == "Global Tech I": return [6.475, 54.5]
    elif name == "Nordsee One": return [6.755, 53.8]
    elif name == "Gemini": return [5.85, 53.93]

def printedWindFarmsNames(name):
    if name == "Gode Wind": return "Gode Wind 1 & 2"
    elif name == "Veja Mate": return "Veja Mate"
    elif name == "Bard Offshore 1": return "Bard"
    elif name == "Riffgat": return "Riffgat"
    elif name == "Borkum Riffgrund 1": return "Borkum\nRiffgrund 1"
    elif name == "Meerwind": return "Meerwind"
    elif name == "Nordsee Ost": return "Nordsee Ost"
    elif name == "Amrumbank West": return "Amrumbank West"
    elif name == "Alpha Ventus": return "Alpha\nVentus"
    elif name ==  "Trianel Windpark Borkum I": return "Trianel\nWindpark\nBorkum I"
    elif name == "Global Tech I": return "Global Tech I"
    elif name == "Nordsee One": return "Nordsee\nOne"
    elif name == "Gemini": return "Gemini"
   
def zoomIdentifier():
    return "Zoom"

def kmToDeg(val):
    return val * 360 / (2 * pi * earthRadiusKm())

def degToKm(val):
    return val / 360 * (2 * pi * earthRadiusKm())

def earthRadiusKm():
    return 6371

def boundsIdentifier():
    return "Boundaries"

def refFileName():
    return "RefFile"

def domainPlotName():
    return "domainPlot"

def loadDefaultColors():
    return ['black', 'red', 'blue', 'green', 'maroon', 'slategrey', 'yellow', 'lime', 'cyan', 'navy', 'purple', 'lightcoral',
    'sienna','peru','tan','orange','darkolivegreen','deeppink','crimson']

def loadIntervalsColors():
    return ["darkred", "tab:blue"]

def locStrToFloat(location):
    splitdata = splitline(location, ",")
    for j in range(2):
        data = splitline(splitdata[j], " ")
        if data[0] in ["N", "E"]:
            s = 1
        elif data[0] in ["S", "W"]:
            s = -1
        else:
            fataError("Unkown direction " + data[0] + ". Supported direction are 'N', 'S', 'E', and 'W'.")
        deg = float(splitline(data[1],"d")[0])
        min = float(splitline(data[2],"m")[0])
        sec = float(splitline(data[3],"s")[0])
        value = (deg + (min + sec / 60.0) / 60.0) * s
        if j == 0: lat = value
        elif j == 1: lon = value 

    return lat, lon

def refWindFarmName():
    return "RefWindFarm"

def closestCell(locLat, locLon, lats, lons):
    minDist = 1000000
    for j in range(lats.shape[0]):
        for i in range(lats.shape[1]):
            ds2 = (lats[j,i]-locLat)**2 + (lons[j,i]-locLon)**2
            if ds2 < minDist:
                minDist = ds2
                isol = i
                jsol = j
    return isol,jsol

def wakeStatsNames():
    return "windFarmWakeStats"

def transectPlotName():
    return "TransectFlightPlot"

def Interpolate1D(array, val):
    kn = nearest(array, val)
    if array[kn] > val:
        ibefore = kn-1
        iafter = kn
    elif array[kn] < val:
        ibefore = kn
        iafter = kn+1
    else:
        return kn, 1

    return ibefore, (array[iafter] - val) / (array[iafter] - array[ibefore])

def horizontalInterpolation(data, lat, lon, lats, lons, ic, jc):
    if lon >= lons[jc,ic]:
        islope = (data[:,jc,ic+1] - data[:,jc,ic]) / ((lons[jc,ic+1] - lons[jc,ic]))
    else:
        islope = (data[:,jc,ic] - data[:,jc,ic-1]) / ((lons[jc,ic] - lons[jc,ic-1]))
    
    if lat >= lats[jc,ic]:
        jslope = (data[:,jc+1,ic] - data[:,jc,ic]) / ((lats[jc+1,ic] - lats[jc,ic]))
    else:
        jslope = (data[:,jc,ic] - data[:,jc-1,ic]) / ((lats[jc,ic] - lats[jc-1,ic]))
    dphi = islope * ((lon - lons[jc,ic])) + jslope * ((lat - lats[jc,ic]))
    return data[:,jc,ic] + dphi[:]

def supportedSegmentKeywords():
    return [lowerTipKeyword(), upperTipKeyword(), hubKeyword(), aboveHubKeyword()]

def lowerTipKeyword():
    return "lowerTip"

def upperTipKeyword():
    return "upperTip"

def hubKeyword():
    return "hub"

def aboveHubKeyword():
    return "aboveHub"

def zMassName():
    return "z"

def zStagName():
    return "zstag"

def temporalPlotName():
    return "TemporalPlot"

def removeChars(string, char):
    return string.replace(char, "")

def closestBelow(vec, val):
    # the vector is assumed to be monotically increasing
    for i in range(len(vec)):
        if vec[i] < val: continue
        return max(0, i - 1)

def metersToIndex(lat, lon, file, ic, jc, limMeters, includeLimits):
    limitIndices = limMeters.copy() * 0
    lats, lons, _, _ = loadGrid(file)
    hStag = horizontalInterpolation(fieldExtract(file, zStagName()), lat, lon, lats, lons, ic, jc)
    for isegment in range(limMeters.shape[0]):
        h1 = limMeters[isegment, 0]
        h2 = limMeters[isegment, 1]
        i1 = includeLimits[isegment, 0]
        i2 = includeLimits[isegment, 1]
        limitIndices[isegment, 0] = closestBelow(hStag, h1) + int(not i1)
        limitIndices[isegment, 1] = closestBelow(hStag, h2) + i2
    return limitIndices

def segmentToMeters(segment, globe):
    b1 = segment[0]
    b2 = segment[-1]
    segment = removeChars(removeChars(segment, b1), b2)
    splitdata = splitline(segment, ",")
    heights = [0, 0]
    for j in range(2):
        # check if it is numeric
        try: heights[j] = float(splitdata[j])
        except:
            # not numeric
            for keyword in supportedSegmentKeywords():
                if not keyword in splitdata[j]: continue
                val = removeChars(removeChars(removeChars(splitdata[j], keyword),"("),")")
                if keyword == aboveHubKeyword():
                    s1 = splitline(val, "/")
                    turbineType = int(s1[0])
                    heights[j] = globe.turbTypeHH(turbineType) + float(s1[1]) * globe.turbineTypeDiameter(turbineType)
                elif keyword == lowerTipKeyword():
                    heights[j] = globe.turbTypeHH(int(val)) - globe.turbineTypeDiameter(int(val)) / 2.0
                elif keyword == upperTipKeyword():
                    heights[j] = globe.turbTypeHH(int(val)) + globe.turbineTypeDiameter(int(val)) / 2.0
                break
    return heights[0], heights[1], int(b1=="["), int(b2=="]")


def q2BudgetColor(field):
    if field == qvTransportName(): return "darkgoldenrod"
    elif field == qsName(): return "darkolivegreen"
    elif field == qbName(): return "darkred"
    elif field == qdName(): return "darksalmon"
    elif field == qtName(): return "darkseagreen"
    else: return "black"

def q2BudgetTag(field):
    if field == qvTransportName(): return "$\mathcal{Q}_w$"
    elif field == qsName(): return "$\mathcal{Q}_s$"
    elif field == qbName(): return "$\mathcal{Q}_b$"
    elif field == qdName(): return "$\mathcal{Q}_d$"
    elif field == qtName(): return "$\mathcal{Q}$"
    else: return "None"

def toLetters(i):
    import string
    return list(string.ascii_lowercase)[i]

def zStagFields():
    return [qvTransportName(), qsName(), qbName(), qdName()]

def vgradPlotName():
    return "speedVerticalGradient"

def message(text, level):
    mes = ""
    for j in range(max(0, level)): mes += tab()
    print(mes + text)

def tab():
    return "    "

def filename(file):
    return splitline(file, "/")[-1]

def newline():
    message("  ", 0)
