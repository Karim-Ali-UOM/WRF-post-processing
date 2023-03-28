from src.orders import order
from src.timeMod import fileToHour, selectTimeFrames
from src.utilities import capsIdentifier, closestCell, createPanel, directionName, durIdentifier, fieldExtract, filename,\
    horizontalInterpolation, listFiles, loadGrid, locIdentifier, locStrToFloat, message, mtendName,\
    nearest, newline, prepareloc, qtName, refFileName, refWindFarmName, rowCol, select, smooth, speedName,\
    tempName, tkename, toLetters, zMassName
from numpy import array
from math import ceil

def verticalFolder():
    return "vertical"

def singleFrameVerticalProfile(file, field, locLat, locLon, lats, lons, ic, jc, reFile):
    return array(horizontalInterpolation(fieldExtract(file, field, reFile), locLat, locLon, lats, lons, ic, jc))

class vertical(order):
    locLat = 0
    locLon = 0
    ic = []
    jc = []

    def __init__(self, name, location, duration):
        self.locLat, self.locLon = locStrToFloat(location)
        self.name = name
        self.intervals.append(duration)
        self.folder = verticalFolder()

    def closestGridCell(self, globe):
        for dir in globe.dirs:
            files = selectTimeFrames(select(listFiles(dir), globe.bases[globe.dirs.index(dir)]), self.intervals[0])
            lats, lons, _, _ = loadGrid(files[0])
            ic, jc  = closestCell(self.locLat, self.locLon, lats, lons)
            self.ic.append(ic)
            self.jc.append(jc)

    def execute(self, globe, field, interval, dir, reFile):
        idir = globe.dirs.index(dir)
        files = selectTimeFrames(select(listFiles(dir), globe.bases[idir]), interval)
        if len(self.ic) == 0: self.closestGridCell(globe)
        lats, lons, _, _ = loadGrid(files[0])
            
        vpAvg = 0
        zholder = singleFrameVerticalProfile(files[0], zMassName(), self.locLat, self.locLon, lats, lons, self.ic[idir],
            self.jc[idir], reFile)
        for file in files:
            vpAvg += singleFrameVerticalProfile(file, field, self.locLat, self.locLon, lats, lons, self.ic[idir], 
            self.jc[idir], reFile)
        dataholder = vpAvg / float(len(files))
        self.writecsv(zholder, dataholder, globe, field, interval.name, dir)
        return zholder, dataholder

    def timeSeries(self, globe, field, dir, reFile, msglvl):
        idir = globe.dirs.index(dir)
        files = selectTimeFrames(select(listFiles(dir), globe.bases[idir]), self.intervals[0])
        if len(self.ic) == 0:
            message("Calculating the closest grid-cell for all directories", msglvl)
            self.closestGridCell(globe)
        lats, lons, _, _ = loadGrid(files[0])
            
        vpAvg = 0
        zholder = singleFrameVerticalProfile(files[0], zMassName(), self.locLat, self.locLon, lats, lons, self.ic[idir],
            self.jc[idir], reFile)
        vpAvg = []
        time = []
        for file in files:
            message("Extracting data from " + filename(file), msglvl)
            vpAvg.append(singleFrameVerticalProfile(file, field, self.locLat, self.locLon, lats, lons, self.ic[idir], 
            self.jc[idir], reFile))
            time.append(fileToHour(file, self.intervals[0].start))

        # sort time array
        message("Sorting time series", msglvl)
        indices = array(time).argsort()
        timeSorted = [time[ind] for ind in indices]
        vpAvgSorted = [vpAvg[ind] for ind in indices]

        return zholder, vpAvgSorted, timeSorted

    def writecsv(self, zholder, dataholder, globe, field, iname, dir):
        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)
        # prepare header
        label = globe.labels[globe.dirs.index(dir)]
        header = "z," + field + "(" + label + "),"

        outfielname = outdir + "/" + self.name + "(" + field + "-" + label + ")-I(" + iname + ")-" + ".csv"
        outfile = open(outfielname, "w")
        outfile.write(header + "\n")
        for iz in range(len(dataholder)):
            line = str(zholder[iz]) + "," + str(dataholder[iz]) + ","
            outfile.write(line + "\n")
        outfile.close

class verticalProfilesPlot(order):
    verticalProfile = vertical
    wfindex = 0
    reFile = str

    def __init__(self, name, identValues, identNames):
        self.name = name
        self.folder = verticalFolder()
        for identifier in identNames:
            k = identNames.index(identifier) 
            if identifier == locIdentifier():
                location = identValues[k]
            elif identifier == durIdentifier():
                duration = identValues[k]
            elif identifier == refWindFarmName():
                self.wfindex = int(identValues[k])
            elif identifier == refFileName():
                self.reFile = identValues[k]
            elif identifier == capsIdentifier():
                self.limits = identValues[k]
        self.verticalProfile = vertical(name, location, duration)

    def execute(self, globe):
        from matplotlib.ticker import FormatStrFormatter

        msglvl = 0
        newline()
        message("Executing " + self.name, msglvl)
        msglvl += 1

        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)
        capsIdentified = False

        # Style
        ticksize = 6
        labelsize = 7
        legendsize = 7
        linewidth = 0.9
        hubwidth = 0.5
        arrowyloc = {speedName():300, directionName():300, tempName():355}
        arrowsize = 65
        ncols = 2

        # Ref turbine hub and diameter
        turbType = globe.windFarmTurbineType(self.wfindex)
        hubheight = globe.turbTypeHH(turbType)
        diameter = globe.turbineTypeDiameter(turbType)

        fields = [mtendName(), qtName(), speedName(), tkename(), directionName(), tempName()]
        symbols = ["$f_U$", "$\mathcal{Q}$", "$-\mathrm{\Delta U}$", "TKE", "$\mathrm{\Delta \phi}$", "$\mathrm{\Delta}$T"]
        units = ["[$\mathrm{cm\,s^{-2}}$]", "[$\mathrm{m^2\,s^{-2}}$]", "[$\mathrm{m\,s^{-1}}$]", "[$\mathrm{m^2\,s^{-2}}$]",
                "[$\mathrm{^\circ}$]", "[$\mathrm{^\circ\!C}$]"]
        colordict = {"A100":"tab:red", "A90":"tab:red", "A80":"tab:red", "F100":"tab:blue", "F25":"tab:blue", "F100-NE":"tab:blue", 
                    "V100":"tab:green", "V120":"tab:green", "V80":"tab:green", "PAN":"tan", "R100":"purple", "R25":"purple",
                    "NT":"black"}
        styledict = {"A100":"-", "A90":"--", "A80":":", "F100":"-", "F25":"--", "F100-NE":":", 
                    "V100":"-", "V120":"--", "V80":":", "PAN":"-", "P25":"--", "R100":"-", "R25":"--", "NT":"-"}
        
        
        nrows = ceil(len(fields) / ncols)
        fig, axgs = createPanel(ncols, nrows, len(fields), 6, 7)
        axgs = axgs.flatten()
        legendHandles = []

        for field in fields:
            message("Field " + field, msglvl)
            msglvl += 1
            ifield = fields.index(field)
            ax = axgs[ifield]
            row, col = rowCol(ifield, ncols)
            if field in [directionName(), speedName(), tempName()]: ax2 = ax.twiny()

            if field == mtendName():
                gain = 100
            else:
                gain = 1

            # reference field
            message("Reference field", msglvl)
            if not field in [qtName(), mtendName()]:
                zref, dataref = self.verticalProfile.execute(globe, field, 
                    self.verticalProfile.intervals[0], globe.dirs[globe.labels.index(globe.ntLabel)],
                    self.reFile)

            for dir in globe.dirs:
                idir = globe.dirs.index(dir)
                label = globe.labels[idir]
                if label == globe.ntLabel and field in [qtName(), mtendName()]: continue
                if field == qtName() and label in ["V80", "V100", "V120"]: continue
                message(label, msglvl)

                if label == globe.ntLabel:
                    if not capsIdentified:
                        capsIdentified = True
                        message("Identify nearest levels", msglvl+1)
                        l1 = nearest(zref, self.limits[0])
                        l2 = nearest(zref, self.limits[1])
                    z2, toplot2 = smooth(array(zref[l1:l2+1]), dataref[l1:l2+1] * gain)
                else:
                    z, data = self.verticalProfile.execute(globe, field, self.verticalProfile.intervals[0], dir, self.reFile)
                    if not capsIdentified:
                        capsIdentified = True
                        message("Identify nearest levels", msglvl+1)
                        l1 = nearest(z, self.limits[0])
                        l2 = nearest(z, self.limits[1])
                    if  field in [directionName(), tempName()]:
                        toplot = array(data) - array(dataref)
                    elif field == speedName():
                        toplot = -array(data) + array(dataref)
                    else:
                        toplot = array(data)

                    
                    z2, toplot2 = smooth(array(z[l1:l2+1]), toplot[l1:l2+1] * gain)
                
                if field == tkename():
                    curve, = ax.plot(toplot2,z2, color=colordict[label], 
                        linestyle = styledict[label], label=label, linewidth = linewidth)
                    legendHandles.append(curve)
                elif field in [directionName(), speedName(), tempName()] and label == globe.ntLabel:
                    ax2.plot(toplot2,z2, color=colordict[label], linestyle = styledict[label], label=label,
                        linewidth = linewidth)
                    arrowxloc = toplot2[nearest(z2,arrowyloc[field])]
                    ax2.annotate("", xy=(arrowxloc, arrowyloc[field]+arrowsize), xycoords='data', 
                        xytext=(arrowxloc, arrowyloc[field]+5), textcoords='data', 
                        arrowprops=dict(arrowstyle="->", connectionstyle="arc3"), 
                        color=colordict[label])
                else:
                    ax.plot(toplot2,z2, color=colordict[label], linestyle = styledict[label], label=label, linewidth = linewidth)
            msglvl -= 1

            message("Formatting figure", msglvl)
            ax.set_xlabel(symbols[ifield] + " " + units[ifield], fontsize = labelsize)
            ax.text(0.91, 0.93, "(" + toLetters(ifield) + ")", fontsize = labelsize, transform=ax.transAxes)
            if field == directionName():
                ax2.set_xlabel("$\mathrm{\phi}$ " + units[ifield], fontsize = labelsize)
            elif field == speedName():
                ax2.set_xlabel("U " + units[ifield], fontsize = labelsize)
            elif field == tempName():
                ax2.set_xlabel("T " + units[ifield], fontsize = labelsize)

            if col == 0: 
                ax.set_ylabel("z [m]", fontsize = labelsize)
            else:
                ax.yaxis.set_ticklabels([])

            ax.set_ylim(0, 420)
            if field == mtendName(): ax.set_xlim(-0.01,0.3)
            elif field == qtName(): ax.set_xlim(-0.002,0.05)
            elif field == speedName():
                ax2.set_xlim(6,14) 
                ax.set_xlim(-0.5,2.1)
                ax2.tick_params(axis='both', labelsize=ticksize)                
            elif field == tkename():
                ax.set_xlim(0,0.8)
            elif field == directionName():
                ax2.set_xlim(200,225)
                ax2.tick_params(axis='both', labelsize=ticksize)
                ax.set_xlim(-3,3)
            elif field == tempName():
                ax2.set_xlim(14,15.7)
                ax2.tick_params(axis='both', labelsize=ticksize)
                ax.set_xlim(-0.4,0.15)
                ax.xaxis.set_major_formatter(FormatStrFormatter("%g"))
                ax2.xaxis.set_major_formatter(FormatStrFormatter("%g"))
            ax.tick_params(axis='both', labelsize=ticksize)

        iax = -1
        for ax in axgs:
            iax += 1
            ax.axhline(y=hubheight, color='grey', linestyle='--', linewidth = hubwidth)
            ax.axhline(y=hubheight+diameter/2.0, color='grey', linestyle='--', linewidth = hubwidth)
            ax.axhline(y=hubheight-diameter/2.0, color='grey', linestyle='--', linewidth = hubwidth)  
            if fields[iax] in [mtendName(), qtName(), tkename()]: continue
            ax.axvline(x=0, color='k', linestyle='--', linewidth = hubwidth)    
        
        message("Saving figure", msglvl)
        fig.subplots_adjust(bottom=0.05, top=0.93, left=0.1, right=0.95, wspace=0.08, hspace=0.45)
        fig.legend(handles = legendHandles, loc = [0.11, 0.94], ncol=7, prop={'size': legendsize})
        fig.savefig(outdir + "/" + self.name + ".pdf", dpi=500)
