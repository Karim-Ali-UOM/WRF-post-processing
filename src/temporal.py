from wrf import interplevel
from src.orders import order
from src.timeMod import timeInstant
from src.utilities import assistFilesIdentifier, createPanel, directionName, durIdentifier, filename, locIdentifier, message, nearest, newline,\
    prepareloc, refFileName, removeChars, rowCol, speedName, splitline, standingAverage, tempName, tkename, toLetters,\
        Interpolate1D
from src.verticalMod import vertical
from math import ceil
from numpy import array

def temporalFolder():
    return "temporal"

class temporalPlot(order):
    files = []
    verticalProfile = vertical 
    reFile = str

    def __init__(self, name, identValues, identNames):
        self.name = name
        self.folder = temporalFolder()
        for identifier in identNames:
            k = identNames.index(identifier)
            if identifier == locIdentifier():
                location = identValues[k]
            elif identifier == durIdentifier():
                duration = identValues[k]
            if identifier == assistFilesIdentifier():
                self.files = identValues[k]
            elif identifier == refFileName():
                self.reFile = identValues[k]
        self.verticalProfile = vertical(name, location, duration)

    def execute(self, globe):
        from matplotlib.ticker import FormatStrFormatter

        msglvl = 0
        newline()
        message("Executing " + self.name, msglvl)
        msglvl += 1

        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)

        #Style
        ncols = 2
        nplots = 8
        linewidth = 0.8
        labelsize = 6
        legendsize = 7
        ticksize = 6
        hubwidth = 0.5
        colordict = {"A100":"tab:red", "A90":"tab:red", "A80":"tab:red", "F100":"tab:blue", "F25":"tab:blue", "F100-NE":"tab:blue", 
                "V100":"tab:green", "V120":"tab:green", "V80":"tab:green", "PAN":"tan", "R100":"purple", "R25":"purple",
                "NT":"black", "EXP":"black"}
        styledict = {"EXP":"--", "A100":"-", "A90":"--", "A80":":", "F100":"-", "F25":"--", "F100-NE":":", 
                "V100":"-", "V120":"--", "V80":":", "PAN":"-", "P25":"--", "R100":"-", "R25":"--", "NT":"-"}

        nrows = ceil(nplots / ncols)
        fig, axgs = createPanel(ncols, nrows, nplots, 6.5, 6)
        axgs = axgs.flatten()
        legendHandles = []
        levelsObtained = False

        for file in self.files:
            # Assumed to be FINO files
            message("Working on " + filename(file), msglvl)
            msglvl += 1
            ifile = self.files.index(file)
            finoFile = open(file, "r")
            allLines = finoFile.readlines()
            lines = []
            header = []
            headerFound = False
            for line in allLines:
                if line[0] == "#": continue
                if splitline(line, "\t")[0] == "Time" and not headerFound: 
                    header = splitline(line, "\t")
                    headerFound = True
                else:
                    lines.append(line)
            data = []
            for line in lines:
                split = splitline(line,"\t")
                timesplit = splitline(split[0]," ")
                s1 = splitline(timesplit[0],"-")
                day = int(s1[2])
                month = int(s1[1])
                year = int(s1[0])
                time = timeInstant(str(day) + "-" + str(month) + "-" + str(year) + ", " + timesplit[1])
                if self.verticalProfile.intervals[0].within(time):
                    data.append([time] + [float(s) for s in split[1:]])
            finoFile.close()
            avgData = standingAverage(data, 10/60.0)
            measTime = [avgData[k][0].hour + avgData[k][0].min / 60.0 + avgData[k][0].sec / 3600 for k in range(len(avgData))]
            
            parts = splitline(splitline(file, "/")[-1],"_")
            if parts[1] == "winddirection":
                fields = [directionName()]
            elif parts[1] == "windspeed":
                fields = [speedName(), tkename()]
            elif parts[1] == "airtemperature":
                fields = [tempName()]
            
            # find height of measurment
            for part in parts:
                try:
                    height = float(removeChars(part, "m"))
                    break
                except: continue

            for field in fields:
                message(field, msglvl)
                msglvl += 1
                ifield = fields.index(field)
                meas = []; minMeas = []; maxMeas = []; devMeas = []; noonloc = -1
                
                if field == speedName():
                    iax = 0
                    ax = axgs[iax]
                    meas = [avgData[k][header.index("Value")] for k in range(len(avgData))]
                    minMeas = [avgData[k][header.index("Minimum")] for k in range(len(avgData))]
                    maxMeas = [avgData[k][header.index("Maximum")] for k in range(len(avgData))]
                    devMeas = [avgData[k][header.index("Deviation")] for k in range(len(avgData))]
                    ax.fill_between(measTime, array(minMeas), array(meas)-array(devMeas), color="slategrey", 
                            alpha=0.4, facecolor="none", linewidth = 0) 
                    mins = ax.fill_between(measTime, array(meas)+array(devMeas), array(maxMeas), color="slategrey", 
                                        alpha=0.4, facecolor="none", linewidth = 0, label="min-max") 
                    devs = ax.fill_between(measTime, array(meas)-array(devMeas), array(meas)+array(devMeas), 
                                    color="navy", alpha=0.5, facecolor="none", linewidth = 0, label="standard deviation") 
                    ax.legend(handles = [devs, mins], loc = [0.05, 0.85], ncol=2, prop={'size': legendsize-2})
                elif field == directionName():
                    iax = 2
                    ax = axgs[iax]
                    meas = [avgData[k][header.index("Value")] for k in range(len(avgData))]
                    devMeas = [avgData[k][header.index("Deviation")] for k in range(len(avgData))]
                    devs = ax.fill_between(measTime, array(meas)-array(devMeas), array(meas)+array(devMeas), 
                            color="navy", alpha=0.5, facecolor="none", linewidth = 0, label="standard deviation")
                    ax.legend(handles = [devs], loc = [0.05, 0.85], ncol=2, prop={'size': legendsize-2})
                elif field == tempName():
                    iax = 4
                    meas = [avgData[k][header.index("Value")] for k in range(len(avgData))]
                elif field == tkename():
                    iax = 6
                    meas = [float("Nan") for k in range(len(avgData))]
                
                ax = axgs[iax]
                
                # scatter measurement
                message("Plotting measurement", msglvl)
                curve = ax.scatter(measTime, meas, color="k", label="FINO-1", s = 2)
                if ifield == 0 and ifile == 0: legendHandles.append(curve)  
                
                # Get the refernece field
                message("Obtaining reference field", msglvl)
                msglvl += 1
                zref, dataref, reftime = self.verticalProfile.timeSeries(globe, field, 
                    globe.dirs[globe.labels.index(globe.ntLabel)], self.reFile, msglvl)
                msglvl -= 1
                if not levelsObtained:
                    levelsObtained = True
                    ib, w = Interpolate1D(zref, height)
                refField = [w * dataref[k][ib] + (1-w) * dataref[k][ib+1] for k in range(len(dataref))]
                curveRef, = ax.plot(reftime, refField, color=colordict[globe.ntLabel], 
                        linestyle = styledict[globe.ntLabel], label=globe.ntLabel, linewidth = linewidth)
                
                for dir in globe.dirs:
                    idir  = globe.dirs.index(dir)
                    label = globe.labels[idir]
                    if label == globe.ntLabel: continue
                    message("Working on " + label, msglvl)
                    msglvl += 1

                    # left column
                    ax = axgs[iax]
                    _, data, time = self.verticalProfile.timeSeries(globe, field, dir, self.reFile, msglvl)
                    msglvl -= 1

                    simField = [w * data[k][ib] + (1-w) * data[k][ib+1] for k in range(len(data))]
                    curve, = ax.plot(time, simField, color=colordict[label], 
                        linestyle = styledict[label], label=label, linewidth = linewidth)  
                    if ifield == 0 and ifile == 0: legendHandles.append(curve) 

                    # right column
                    ax = axgs[iax + 1]
                    noonloc = nearest(time, 12)
                    time2 = array(time[noonloc:])
                    if field == speedName():
                        curve, = ax.plot(time2, (array(refField)-array(simField))[noonloc:],
                            color=colordict[label], linestyle = styledict[label], label=label, 
                            linewidth = linewidth)  
                    else:
                        curve, = ax.plot(time2, (array(simField)-array(refField))[noonloc:],
                            color=colordict[label], linestyle = styledict[label], label=label, 
                            linewidth = linewidth) 

                if ifield == 0 and ifile == 0: legendHandles.append(curveRef) 

                # styling
                # left column      
                message("Formatting figure", msglvl) 
                ax = axgs[iax]
                row, _ = rowCol(iax, ncols)
                if field == speedName():
                    ax.set_ylabel("U [$\mathrm{m\,s^{-1}}$]", fontsize = labelsize)
                    ax.set_ylim(5,20)
                elif field == directionName():
                    ax.set_ylabel("$\mathrm{\phi}$ [$\mathrm{^\circ}$]", fontsize = labelsize)
                    ax.set_ylim(205,245)
                    ax.set_yticks([205,213,221,229,237,245])
                elif field == tempName():
                    ax.set_ylabel("T [$\mathrm{^\circ C}$]", fontsize = labelsize)
                    ax.set_ylim(14.5,16)
                    ax.set_yticks([14.5,14.75,15,15.25,15.5,15.75,16])
                elif field == tkename():
                    ax.set_ylabel("TKE [$\mathrm{m^2\,s^{-2}}$]", fontsize = labelsize)
                    ax.set_ylim(0,1.3)
                    ax.set_yticks([0,0.2,0.4,0.6,0.8,1.0,1.2])
                ax.set_xlim(-0.5,24.5)
                ax.set_xticks([0,2,4,6,8,10,12,14,16,18,20,22,24])
                if row == nrows-1:
                    ax.set_xlabel("Time [hr]", fontsize = labelsize)
                else:
                    ax.xaxis.set_ticklabels([])
                ax.tick_params(axis='both', labelsize=ticksize)
                ax.text(0.93, 0.89, "(" + toLetters(iax) + ")", fontsize = labelsize, transform=ax.transAxes)

                # right column: relative to NT
                ax = axgs[iax+1]
                ax.yaxis.set_label_position("right")
                ax.yaxis.tick_right()
                ax.yaxis.set_major_formatter(FormatStrFormatter("%g"))
                if field == speedName():
                    ax.set_ylabel("$\mathrm{-\Delta U}$ [$\mathrm{m\,s^{-1}}$]", fontsize = labelsize)
                    ax.set_ylim(0.5,2.5)
                elif field == directionName():
                    ax.set_ylabel("$\mathrm{\Delta \phi}$ [$\mathrm{^\circ}$]", fontsize = labelsize)
                    ax.set_ylim(-5,3.2)
                    ax.set_yticks([-5,-4,-3,-2,-1,0,1,2,3])
                    ax.axhline(y=0, color='grey', linestyle='--', linewidth = hubwidth)
                elif field == tempName():
                    ax.set_ylabel("$\mathrm{\Delta T}$ [$\mathrm{^\circ C}$]", fontsize = labelsize)
                    ax.set_ylim(-0.25,0.25)
                    ax.set_yticks([-0.25, -0.125, 0, 0.125, 0.25])
                    ax.set_yticklabels(["$\mathrm{-1/4}$", "$\mathrm{-1/8}$", "$\mathrm{0}$",
                                    "$\mathrm{1/8}$", "$\mathrm{1/4}$"])
                    ax.axhline(y=0, color='grey', linestyle='--', linewidth = hubwidth)
                elif field == tkename():
                    ax.set_ylabel("$\mathrm{\Delta}$TKE [$\mathrm{m^2\,s^{-2}}$]", fontsize = labelsize)
                    ax.set_ylim(-0.1,0.8)
                    ax.set_yticks([0,0.2,0.4,0.6,0.8])
                ax.set_xlim(11.5,24.5)
                ax.set_xticks([12,13,14,15,16,17,18,19,20,21,22,23,24])
                if row == nrows-1:
                    ax.set_xlabel("Time [hr]", fontsize = labelsize)
                else:
                    ax.xaxis.set_ticklabels([])
                ax.tick_params(axis='both', labelsize=ticksize)
                ax.text(0.93, 0.89, "(" + toLetters(iax+1) + ")", fontsize = labelsize, transform=ax.transAxes)
            msglvl -= 1
        msglvl -= 1
        message("Saving figure", msglvl)
        fig.legend(handles = legendHandles, loc = [0.115, 0.93], ncol=7, prop={'size': legendsize})
        fig.subplots_adjust(bottom=0.075, top=0.92, left=0.1, right=0.9, wspace=0.05, hspace=0.1)
        fig.savefig(outdir + "/" + self.name + ".pdf")
        fig.clf()
        
