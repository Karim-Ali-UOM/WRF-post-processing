from src.orders import order
from src.timeMod import clockFormat, selectNearestTimeFrames, timeInstant
from src.utilities import Interpolate1D, assistFilesIdentifier, closestCell, createPanel, filename,\
    listFiles, loadGrid, message, newline, prepareloc, removeChars, rowCol, select, speedName, splitline,\
    tkename, toLetters, zMassName
from glob import glob
from numpy import sqrt
from src.verticalMod import singleFrameVerticalProfile

def transectFolder():
    return "transect"

class TransectFlightsPlot(order):
    files = []

    def __init__(self, name, identValues, identNames):
        self.name = name
        self.folder = transectFolder()
        for identifier in identNames:
            k = identNames.index(identifier)
            if identifier == assistFilesIdentifier():
                dynFiles = identValues[k]
        try:
            for dynFile in dynFiles:
                files = glob(dynFile)
                for file in files: self.files.append(file)
        except:
            self.files = []

    def writecsv(self, flightIndex, droneClock, droneLat, droneLon, droneAlt, droneSpeed, droneTKE, simSpeeds, simTKEs, globe):
        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)

        outfile = open(outdir + "/Speed-transect" + str(flightIndex) + ".csv", "w")
        line = "Time, Lat, Lon, Alt[m], OBS,"
        for label in globe.labels: line += label + ","
        outfile.write(line + "\n")
        for k in range(len(droneClock)):
            line = str(droneClock[k]) + "," + str(droneLat[k]) + "," + str(droneLon[k]) + "," + str(droneAlt[k]) + "," + str(droneSpeed[k]) + ","
            for j in range(len(globe.labels)): line += str(simSpeeds[j][k]) + ","
            outfile.write(line + "\n")
        outfile.close()

        outfile = open(outdir + "/TKE-transect" + str(flightIndex) + ".csv", "w")
        line = "Time, Lat, Lon, Alt[m], OBS,"
        for label in globe.labels: line += label + ","
        outfile.write(line + "\n")
        for k in range(len(droneClock)):
            line = str(droneClock[k]) + "," + str(droneLat[k]) + "," + str(droneLon[k]) + "," + str(droneAlt[k]) + "," + str(droneTKE[k]) + ","
            for j in range(len(globe.labels)): line += str(simTKEs[j][k]) + ","
            outfile.write(line + "\n")
        outfile.close()
                
    def execute(self, globe):
        msglvl = 0
        newline()
        message("Executing " + self.name, msglvl)
        msglvl += 1
        
        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)

        # style
        ncols = 2
        nrows = len(self.files)
        fig, axgs = createPanel(ncols, nrows, len(self.files)*2, 6, 7)
        axgs = axgs.flatten()
        legendHandles = []
        linewidth = 0.9
        labelsize = 6
        legendsize = 7
        ticksize = 6
        colordict = {"A100":"tab:red", "A90":"tab:red", "A80":"tab:red", "F100":"tab:blue", "F25":"tab:blue", "F100-NE":"tab:blue", 
                "V100":"tab:green", "V120":"tab:green", "V80":"tab:green", "PAN":"tan", "R100":"purple", "R25":"purple",
                "NT":"black", "OBS":"black"}
        styledict = {"OBS":"--", "A100":"-", "A90":"--", "A80":":", "F100":"-", "F25":"--", "F100-NE":":", 
                "V100":"-", "V120":"--", "V80":":", "PAN":"-", "R100":"-", "R25":"--", "NT":"-"}

        # Assume same grid for all dirs to save time. modify this when different grids as used
        wrfFiles = select(listFiles(globe.dirs[0]), globe.bases[0])
        lats, lons, _, _ = loadGrid(wrfFiles[0])
        
        # read drone files
        for file in self.files:
            message("Working on " + filename(file), msglvl)
            msglvl += 1
            ifile = self.files.index(file)

            cutName = splitline(splitline(file,"/")[-1],"-")
            for part in cutName:
                try:
                    flightIndex = int(removeChars(part, "transect"))
                    break
                except: continue
            timeStamp = splitline(cutName[-1], "_")[0]
            day = timeStamp[6:8] + "-" + timeStamp[4:6] + "-" +  timeStamp[0:4]
            droneClock = []; droneLat = []; droneLon = []; droneAlt = []; droneSpeed = []; droneTKE = []
            filehandle = open(file, "r")
            alldata = filehandle.readlines()
            for line in alldata:
                splitdata = splitline(line, ",")
                try: vals = [float(splitdata[k]) for k in range(len(splitdata))]
                except:
                    header = splitdata
                    continue
                droneClock.append(vals[header.index("Time")])
                droneLat.append(vals[header.index("Lat[deg]")])
                droneLon.append(vals[header.index("Lon[deg]")])
                droneAlt.append(vals[header.index("Alt[m]")])
                droneSpeed.append(sqrt(vals[header.index("u[m/s]")]**2 + vals[header.index("v[m/s]")]**2))
                droneTKE.append(vals[header.index("TKE[m2/s2]")])

            simSpeeds = [[] for k in globe.dirs]
            simTKEs = [[] for k in globe.dirs]

            for i in range(len(droneClock)):
                clock = droneClock[i]
                lat = droneLat[i]
                lon = droneLon[i]
                alt = droneAlt[i]
                message("Clock: " + str(clock), msglvl)
                msglvl += 1
                
                ic, jc = closestCell(lat, lon, lats, lons)
                for dir in globe.dirs:
                    idir = globe.dirs.index(dir)
                    message(globe.labels[idir], msglvl)
                    wrfFiles = select(listFiles(dir), globe.bases[idir])
                    nearestFiles, weights = selectNearestTimeFrames(wrfFiles, timeInstant(day +  ", " + clockFormat(clock)))
                   
                    #lats, lons, _, _ = loadGrid(nearestFiles[0])
                    #ic, jc = closestCell(lat, lon, lats, lons)
                    z = singleFrameVerticalProfile(nearestFiles[0], zMassName(), lat, lon, lats, lons, ic, jc, None)
                    ib, w = Interpolate1D(z, alt)
                    
                    speedProfile = 0
                    tkeProfile = 0
                    for timeFrame in nearestFiles:
                        speedProfile += weights[nearestFiles.index(timeFrame)] * singleFrameVerticalProfile(timeFrame, 
                            speedName(), lat, lon, lats, lons, ic, jc, None)
                        tkeProfile += weights[nearestFiles.index(timeFrame)] * singleFrameVerticalProfile(timeFrame, 
                            tkename(), lat, lon, lats, lons, ic, jc, None)

                    simSpeeds[idir].append(w * speedProfile[ib] + (1-w) * speedProfile[ib+1])
                    simTKEs[idir].append(w * tkeProfile[ib] + (1-w) * tkeProfile[ib+1])
                msglvl -= 1
            msglvl -= 1
            
            # export csv
            message("Exporting results to csv", msglvl)
            self.writecsv(flightIndex, droneClock, droneLat, droneLon, droneAlt, droneSpeed, droneTKE, simSpeeds, simTKEs, globe)

            message("Formatting figure", msglvl)
            iax = int((flightIndex - 1) * 2)
            # Left column for speed
            ax = axgs[iax]
            row, _ = rowCol(iax, ncols)
            curve, = ax.plot(droneLat, droneSpeed, color=colordict["OBS"], linestyle = styledict["OBS"],
                    label= "OBS", linewidth = linewidth)
            if ifile == 0: legendHandles.append(curve)
            for idir in range(len(globe.dirs)):
                label = globe.labels[idir]
                curve, = ax.plot(droneLat, simSpeeds[idir], color=colordict[label], linestyle = styledict[label],
                    label= label, linewidth = linewidth)
                if ifile == 0: legendHandles.append(curve)
            ax.text(0.91, 0.88, "(" + toLetters(iax) + ")", fontsize = labelsize, transform=ax.transAxes)
            ax.grid(False)
            ax.set_ylim(12,17)
            ax.set_yticks([12,13,14,15,16,17])
            ax.set_xticks([53.9,53.95,54,54.05,54.1,54.15,54.2,54.25])
            # ax.axvline(x = 53.95, color = "grey", linestyle = "dotted")
            # ax.axvline(x = 54.15, color = "grey", linestyle = "dotted")
            ax.axvspan(53.95, 54.15, alpha=0.2, color='grey')
            ax.tick_params(axis='both', labelsize=ticksize)
            ax.set_ylabel("U [$\mathrm{m\,s^{-1}}$]", fontsize = labelsize)
            if row == nrows-1:
                ax.set_xlabel("Lat [$^{\circ}$]", fontsize = labelsize)
            else:
                ax.xaxis.set_ticklabels([])

            # Right column for TKE
            ax = axgs[iax + 1]
            curve, = ax.plot(droneLat, droneTKE, color=colordict["OBS"], linestyle = styledict["OBS"],
                    label= "OBS", linewidth = linewidth)
            for idir in range(len(globe.dirs)):
                label = globe.labels[idir]
                curve, = ax.plot(droneLat, simTKEs[idir], color=colordict[label], linestyle = styledict[label],
                    label= label, linewidth = linewidth)
            ax.text(0.91, 0.88, "(" + toLetters(iax+1) + ")", fontsize = labelsize, transform=ax.transAxes)
            ax.set_ylabel("TKE [$\mathrm{m^2\,s^{-2}}$]", fontsize = labelsize)
            ax.tick_params(axis='both', labelsize=ticksize)
            ax.set_ylim(0,3)
            ax.set_yticks([0,0.5,1,1.5,2,2.5,3])
            ax.set_xticks([53.9,53.95,54,54.05,54.1,54.15,54.2,54.25])
            ax.axvspan(53.95, 54.15, alpha=0.2, color='grey')
            if row == nrows-1:
                ax.set_xlabel("Lat [$^{\circ}$]", fontsize = labelsize)
            else:
                ax.xaxis.set_ticklabels([])

        message("Saving figure", msglvl)
        fig.legend(handles = legendHandles, loc = [0.11, 0.94], ncol=7, prop={'size': legendsize})
        fig.subplots_adjust(bottom=0.05, top=0.93, left=0.1, right=0.95, wspace=0.225, hspace=0.2)
        fig.savefig(outdir + "/" + self.name + ".pdf")
        fig.clf()
