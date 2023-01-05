from src.timeMod import selectTimeFrames
from src.utilities import createPanel, fataError, fieldExtract, filename, intervalsIdentifier, listFiles, loadIntervalsColors, message, newline, pairsIdentifier, powerfieldName, prepareloc, select, checkField, calcTurbLocIdentifier, readTurbLocIdentifier,\
    createIndexMap, loadGrid, splitline, turbinePowerFieldname, wfIndexIdentifier
import numpy as np
from src.orders import order
from math import ceil

def powerFolder():
    return "power"

def resolveIntersection(minHolder, maxHolder):
    #timeAveragedStats[:,farmIndex,:,0]
    whichStop = np.zeros([minHolder.shape[1], minHolder.shape[0], 2]) - 1 # nIntervals x nModels x 2 (for min and max)
    for km in range(minHolder.shape[0]):
        for i in range(minHolder.shape[1]):
            for j in [k+i+1 for k in range(minHolder.shape[1]-i-1)]:
                if maxHolder[km,i] > maxHolder[km,j] and minHolder[km,i] > minHolder[km,j] and minHolder[km,i] < maxHolder[km,j]:
                    # min i and max j overlap
                    whichStop[i, km, 0] = j
                    whichStop[j, km, 1] = i

    return whichStop

def timeAveragedPowerStats(files, globe, msglvl):
    minperwf = []
    maxperwf = []
    avgperwf = []
    stdperwf = []
    totperwf = []
    MinPs = [[] for k in files]
    MaxPs = [[] for k in files]
    AvgPs = [[] for k in files]
    stdPs = [[] for k in files]
    totPs = [[] for k in files]
    stored_i = []
    stored_j = []
    for file in files:
        minP, maxP, avgP, stdP, totP, stored_i, stored_j = singleFramePowerStats(file, stored_i, stored_j, globe, msglvl)
        MinPs[files.index(file)] = minP
        MaxPs[files.index(file)] = maxP
        AvgPs[files.index(file)] = avgP
        stdPs[files.index(file)] = stdP
        totPs[files.index(file)] = totP

    for k in range(globe.numWindFarms):
        try:
            minperwf.append(np.mean(np.array(MinPs)[:,k]))
            maxperwf.append(np.mean(np.array(MaxPs)[:,k]))
            avgperwf.append(np.mean(np.array(AvgPs)[:,k]))
            stdperwf.append(np.mean(np.array(stdPs)[:,k]))
            totperwf.append(np.mean(np.array(totPs)[:,k]))
        except: continue
    
    return minperwf, maxperwf, avgperwf, stdperwf, totperwf

def singleFramePowerStats(file, stored_i, stored_j, globe, msglvl):

    message("Extracting power stats from " + filename(file), msglvl)
    msglvl += 1

    windFarmPowerCells = [[] for k in range(globe.numWindFarms)]
    mins = []
    maxs = []
    avgs = []
    stds = []
    tots = []

    if checkField(file,turbinePowerFieldname()):
        # this is more accurate option by counting power generation from all
        # turbines belonging to this windfarm
        singlePowerField = fieldExtract(file, turbinePowerFieldname())
        
        if len(stored_j) == 0:
            message("Storing the indices of grid-cells with turbines", msglvl)
            for i in range(singlePowerField.shape[0]):
                for j in range(singlePowerField.shape[1]):
                    if singlePowerField[i,j] == 0: continue
                    stored_i.append(i)
                    stored_j.append(j)

        for ij in range(len(stored_j)):
            i = stored_i[ij]
            j = stored_j[ij]
            if singlePowerField[i,j] < 0: wfIndex = globe.turbineIndexToWindFarmIndex(abs(singlePowerField[i,j])) - 1
            elif singlePowerField[i,j] > 0: windFarmPowerCells[wfIndex].append(singlePowerField[i,j])
        
        for wf in globe.windfarmIndices:
            k = globe.windfarmIndices.index(wf)
            if len(windFarmPowerCells[k]) == 0: continue
        #    if len(windFarmPowerCells[k]) != globe.windFarmNumTurbinesDict[wf]:
        #        message("Mis-match in number of turbines for wind farm " + str(wf) + 
        #        " (" + str(globe.windfarmsNamesDict[wf]) + ").", msglvl)
            mins.append(np.min(windFarmPowerCells[k]))
            maxs.append(np.max(windFarmPowerCells[k]))
            avgs.append(np.mean(windFarmPowerCells[k]))
            stds.append(np.std(windFarmPowerCells[k]))
            tots.append(np.sum(windFarmPowerCells[k]))

        msglvl -= 1
        return mins, maxs, avgs, stds, tots, stored_i, stored_j
        
    else:
        fataError("Cannot find " + turbinePowerFieldname() + " in " + file)
        
def statsNames():
    return ["min", "max", "mean", "std", "total"]

class powerStat(order):
    
    def __init__(self, name, wfs, intervals):
        self.name = name
        self.windfarms = wfs.copy()
        self.intervals = intervals
        self.folder = powerFolder()

    def writecsv(self, dataHolder, globe):
        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)
        statnames = statsNames()

        header = "Case,"
        for interval in self.intervals: header += interval.name + ","

        for wf in self.windfarms:
            outfile = open(outdir + "/" + globe.windFarmName(wf) + "-" + self.name + ".csv", "w")
            outfile.write(header + "\n")

            for label in globe.labels:
                for k in range(5):
                    line = label  + "(" + statnames[k] + "),"
                    for j in range(dataHolder.shape[2]): line += str(dataHolder[globe.labels.index(label),globe.windfarmIndices.index(wf),j,k]) + ","
                    outfile.write(line + "\n")
            outfile.close()

    def execute(self, globe, msglvl):
        dataHolder = np.zeros([len(globe.dirs), globe.numWindFarms, len(self.intervals), 5])

        for dir in globe.dirs:
            caseindex = globe.dirs.index(dir)
            message("Calculating power stats for " + globe.labels[caseindex], msglvl)
            msglvl += 1
            for interval in self.intervals:
                intervalIndex = self.intervals.index(interval)
                message("Calculating power stats for interval " + str(intervalIndex+1), msglvl)
                msglvl += 1
                files = selectTimeFrames(select(listFiles(dir), globe.bases[caseindex]), interval)            
                minP, maxP, avgP, stdP, totP = timeAveragedPowerStats(files, globe, msglvl)
                for wfIndex in range(globe.numWindFarms):
                    try:
                        dataHolder[caseindex, wfIndex, intervalIndex, 0] = minP[wfIndex]
                        dataHolder[caseindex, wfIndex, intervalIndex, 1] = maxP[wfIndex]
                        dataHolder[caseindex, wfIndex, intervalIndex, 2] = avgP[wfIndex]
                        dataHolder[caseindex, wfIndex, intervalIndex, 3] = stdP[wfIndex]
                        dataHolder[caseindex, wfIndex, intervalIndex, 4] = totP[wfIndex]
                    except: continue
                msglvl -= 1
            msglvl -= 1

        message("Exporting power stats to csv", msglvl)
        self.writecsv(dataHolder, globe)
        return dataHolder

class capacityFactorPlot(order):
    stats = powerStat
    pairs = []

    def __init__(self, name, identValues, identNames):
        self.name = name
        self.folder = powerFolder()
        for identifier in identNames:
            k = identNames.index(identifier)
            if identifier == intervalsIdentifier():
                self.intervals = identValues[k].copy()
            elif identifier == wfIndexIdentifier():
                self.windfarms = identValues[k].copy()
            elif identifier == pairsIdentifier():
                self.pairs = identValues[k]
        self.stats = powerStat(name, self.windfarms, self.intervals)
        
    def calcPairs(self, timeAveragedStats, globe, msglvl):
        if len(self.pairs) == 0: return
        message("Calculating power stats pairs", msglvl)
        msglvl += 1

        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)
        outfile = open(outdir + "/Pairs-" + self.name + ".csv", "w")
        dataholder = np.zeros([len(self.intervals), len(self.pairs)])
        header = "Interval,"
        wfIndices = [wf-1 for wf in self.windfarms]
        for pair in self.pairs:
            message(pair, msglvl)
            splitdata = splitline(pair, ",")
            caseLabel = splitdata[0]
            refLabel = splitdata[1]
            header += caseLabel + "/" + refLabel + ","
            caseindex = globe.labels.index(caseLabel) 
            refindex = globe.labels.index(refLabel) 
            for kinterval in range(len(self.intervals)):
                caseValue = np.mean(timeAveragedStats[caseindex,wfIndices,kinterval,4])
                refValue = np.mean(timeAveragedStats[refindex,wfIndices,kinterval,4])
                dataholder[kinterval,self.pairs.index(pair)] = (caseValue - refValue) / refValue * 100

        outfile.write(header+ "\n")
        for kinterval in range(len(self.intervals)):
            line = self.intervals[kinterval].name + ","
            for ipair in range(len(self.pairs)): line += str(dataholder[kinterval,ipair]) + ","
            outfile.write(line + "\n")
        outfile.close()
        msglvl -= 1

    def execute(self, globe):
        msglvl = 0
        newline()
        message("Executing " + self.name, msglvl)
        msglvl += 1
        
        timeAveragedStats = self.stats.execute(globe, msglvl)
        self.calcPairs(timeAveragedStats, globe, msglvl)

        message("Creating the figure", msglvl)
        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)
        labels = globe.labels
        ntIndex = globe.labels.index(globe.ntLabel)
        labels.pop(ntIndex)

        # Style
        ytickslabels = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        yticks = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
        limits = [0.1,1.05]
        ncols = 3
        labelsize = 6
        legendsize = 7
        circleSize = 14
        barwidth = 0.25
        intervalColors = loadIntervalsColors()
        legendHandles = []

        nrows = ceil(len(self.windfarms)/ncols)
        fig, axgs = createPanel(ncols, nrows, len(self.windfarms), 6, 7)
        axgs = axgs.flatten()

        namesindex = [i+1 for i in range(len(labels))]

        for wf in self.windfarms:
            plotIndex = self.windfarms.index(wf)
            farmIndex = globe.windfarmIndices.index(wf)
            wfRatedPower = globe.farmRatedPower(globe.windFarmName(wf)) * 1000000
            tubineRatedPower = globe.turbineTypePower(globe.windFarmTurbineType(wf))
            ax = axgs[plotIndex]
            row = plotIndex // ncols
            col = plotIndex - row * ncols

            whichStop = resolveIntersection(np.array(timeAveragedStats[:,farmIndex,:,0]), np.array(timeAveragedStats[:,farmIndex,:,1]))

            for kinterval in range(timeAveragedStats.shape[2]):
                mins = np.array(timeAveragedStats[:,farmIndex,kinterval,0]) / tubineRatedPower
                maxs = np.array(timeAveragedStats[:,farmIndex,kinterval,1]) / tubineRatedPower
                avgs = np.array(timeAveragedStats[:,farmIndex,kinterval,4]) / wfRatedPower

                mins = np.delete(mins, ntIndex)
                maxs = np.delete(maxs, ntIndex)
                avgs = np.delete(avgs, ntIndex)
                nt = np.array(timeAveragedStats[ntIndex,farmIndex,kinterval,4]) / wfRatedPower

                for km in range(len(labels)):
                    # check if min intersection
                    if whichStop[kinterval, km, 0] != -1:
                        minVal = timeAveragedStats[km,farmIndex,int(whichStop[kinterval, km, 0]), 1] / tubineRatedPower
                    else:
                        minVal = mins[km]

                    # check if max intersection
                    if whichStop[kinterval, km, 1] != -1:
                        maxVal = timeAveragedStats[km,farmIndex,int(whichStop[kinterval, km, 1]), 0] / tubineRatedPower
                    else:
                        maxVal = maxs[km]

                    handle = ax.vlines(x = namesindex[km], ymin=minVal, ymax=maxVal, zorder = 50, linewidth=1.5, 
                        color=intervalColors[kinterval], label = self.intervals[kinterval].name)
                    if plotIndex == 0 and km == 0: legendHandles.append(handle)

                ax.scatter(namesindex, avgs, facecolor = "white", edgecolor = intervalColors[kinterval], s = circleSize, zorder = 100)
                ax.axhline(y = nt, color = intervalColors[kinterval], linestyle = "--", linewidth=1)

            # After main plots, loop to plot intersection connectors
            for kinterval in range(timeAveragedStats.shape[2]):
                for km in range(len(labels)):
                    if whichStop[kinterval, km, 0] != -1:
                        top = timeAveragedStats[km,farmIndex,int(whichStop[kinterval, km, 0]), 1] / tubineRatedPower
                        bot = timeAveragedStats[km,farmIndex,kinterval, 0] / tubineRatedPower
                        ax.vlines(x=namesindex[km], ymin=bot, ymax=top, linewidth=1, color="k", linestyle = "dotted")
                        ax.hlines(y = top, xmin = namesindex[km]-barwidth/1,  xmax = namesindex[km]+barwidth/1, color = intervalColors[int(whichStop[kinterval, km, 0])], linewidth=1.5, zorder = 60)
                        ax.hlines(y = bot, xmin = namesindex[km]-barwidth/1,  xmax = namesindex[km]+barwidth/1, color = intervalColors[kinterval], linewidth=1.5, zorder = 60)
                        whichStop[int(whichStop[kinterval, km, 0]), km, 1] = -1

            ax.set_xticks(namesindex)
            if row == nrows-1: ax.set_xticklabels(labels, fontsize = labelsize, rotation = 75) 
            else: ax.xaxis.set_ticklabels([])
            
            ax.set_ylim(limits)
            ax.set_yticks(yticks)
            ax.set_yticklabels(ytickslabels, fontsize = labelsize) 
            if col != 0: ax.yaxis.set_ticklabels([])

            ax.set_title(globe.windFarmName(wf), fontsize = labelsize)  

        fig.legend(handles = legendHandles, loc = [0.33, 0.96], ncol=max(4,len(self.intervals)), prop={'size': legendsize})
        fig.subplots_adjust(bottom=0.075, top=0.93, left=0.05, right=0.975, wspace=0.05, hspace=0.2)
        fig.savefig(outdir + "/" + self.name + ".pdf")
        fig.clf()