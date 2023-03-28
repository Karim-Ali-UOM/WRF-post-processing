from src.orders import order
from src.timeMod import selectNearestTimeFrames, selectTimeFrames, timeInstant, timeInterval
from numpy import array, min, max, arange, sum, linalg, cross, dot, zeros, divide
from wrf import interplevel, to_np
from src.utilities import createCrsPanel, createPanel, degToKm, durIdentifier, dxname, extractDirectVariable,\
    field_diff_percent, field_diff_percent2, fieldExtract, heightIdentifier, levelsIdentifier, listFiles, loadGrid,\
    markWindFarmCells, message, newline, pairsIdentifier, prepareloc, refWindFarmName, rowCol, select, speed10Name,\
        speedName, splitline, timeInstantIdentifier, toLetters, zoomIdentifier
from math import ceil

def wakeFolder():
    return "wake"

class wakeContourPlot(order):
    instant = timeInstant
    levels = []
    zoom = []
    height = 0

    def __init__(self, name, identValues, identNames):
        self.name = name
        self.folder = wakeFolder()
        for identifier in identNames:
            k = identNames.index(identifier)
            if identifier == timeInstantIdentifier():
                self.instant = identValues[k]
            elif identifier == levelsIdentifier():
                self.levels = identValues[k]
            elif identifier == zoomIdentifier():
                self.zoom = identValues[k]
            elif identifier == heightIdentifier():
                self.height = identValues[k]
                
    def execute(self, globe):

        from matplotlib.cm import get_cmap
        import cartopy.crs as crs
        from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter)
        import matplotlib.colors as clr

        msglvl = 0
        newline()
        message("Executing " + self.name, msglvl)
        msglvl += 1

        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)

        # These assisting boxes are hardwired to remove background noise
        # in regions away from the wind farms.
        # If other scenarios are to be exmained s=using this code, either
        # comment out these lines or change accordingly
        boxes = array([[7.5,100,-100,54.1],
                [7,100,-100,53.8],
                [-100,6.4,-100,53.8],
                [-100,7,54.7,100],
                [-100,100,-100,53.6]])

        maxwake = max(self.levels)
        minwake = min(self.levels)
        color = get_cmap("viridis")
        crss = crs.PlateCarree()
        plotsPerRow = 3
        nRows = ceil((len(globe.labels)-1) / plotsPerRow)  
        fig, axGlobal = createCrsPanel(plotsPerRow, nRows, len(globe.labels) - 1, 5.5, 3.2, crss)
        axGlobal = axGlobal.flatten()
        labels = []

        # find wind-farms footprint using the power field in any file
        files = select(listFiles(globe.dirs[0]), globe.bases[0])
        lats, lons, wfCells = markWindFarmCells(files[0])

        # Reference field
        iref = globe.labels.index(globe.ntLabel)
        files = select(listFiles(globe.dirs[iref]), globe.bases[iref])
        selectedTimeFiles, weights = selectNearestTimeFrames(files, self.instant)
        # since all cases will be compared to this refernece case, 
        # it is assumed that all cases have the same grid
        lats, lons, z, _ = loadGrid(selectedTimeFiles[0]) 
        
        refField = 0
        message("Reference field", msglvl)
        for file in selectedTimeFiles:
            if self.height == 10:
                refField += fieldExtract(file, speed10Name()) * weights[selectedTimeFiles.index(file)]
            else:
                refField += to_np(interplevel(fieldExtract(file, speedName()), z, self.height)) * weights[selectedTimeFiles.index(file)]

        iplot = -1
        for dir in globe.dirs:
            idir = globe.dirs.index(dir)
            if idir == iref: continue
            message(globe.labels[idir], msglvl)
            iplot += 1
            ax = axGlobal[iplot]
            labels.append(globe.labels[idir])
            files = select(listFiles(dir), globe.bases[idir])
            selectedTimeFiles, weights = selectNearestTimeFrames(files, self.instant)
            sol = 0
            for file in selectedTimeFiles:
                if self.height == 10:
                    sol += fieldExtract(file, speed10Name()) * weights[selectedTimeFiles.index(file)]
                else:
                    sol += to_np(interplevel(fieldExtract(file, speedName()), z, self.height)) * weights[selectedTimeFiles.index(file)]

            relDiff = field_diff_percent2(sol, refField)
            for j in range(lons.shape[0]):
                for i in range(lons.shape[1]):
                    if relDiff[j,i] > maxwake: relDiff[j,i] = maxwake
                    elif relDiff[j,i] < minwake: relDiff[j,i] = float("Nan")
                    for b in range(boxes.shape[0]):
                        if (boxes[b,0] <= lons[j,i] < boxes[b,1]) and (boxes[b,2] <= lats[j,i] < boxes[b,3]):
                            relDiff[j,i] = float("Nan")
            im = axGlobal[iplot].contourf(lons, lats, relDiff, self.levels, cmap=color, transform=crss)
            ax.contourf(lons, lats, wfCells, [0.9, 1], transform=crss, 
                cmap = clr.ListedColormap(['black', 'black']), zorder=20, alpha = 1) 
    
        message("Formatting figure", msglvl)
        iax = -1
        for ax in axGlobal:
            iax += 1
            row, col = rowCol(iax, plotsPerRow)

            ax.yaxis.tick_left()
            lonStep = 0.5
            latStep = 0.3
            ax.set_xticks([round(q,1) for q in arange(self.zoom[0]+lonStep, self.zoom[1], lonStep)], crs=crss)
            ax.set_yticks([round(q,1) for q in arange(self.zoom[2]+latStep, self.zoom[3] - latStep, latStep)], crs=crss)
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            ax.yaxis.set_major_formatter(lat_formatter)
            ax.tick_params(axis='both', which='major', labelsize=5)
            ax.grid(False)

            ax.text(5.6, 54.75, "(" + toLetters(iax) + ") " + labels[iax], fontsize = 5)
            ax.set_extent(self.zoom)
            if row != nRows-1: ax.xaxis.set_ticklabels([])
            if col != 0: ax.yaxis.set_ticklabels([])

        cb_ax = fig.add_axes([0.2, 0.06, 0.6, 0.02])
        cbar = fig.colorbar(im, cax=cb_ax, orientation='horizontal', format='%.0f%%')
        cbar.ax.tick_params(labelsize=5) 
        ticks = self.levels
        cbar.set_ticks(ticks)
        cbar.ax.set_title(" ")

        message("Saving figure", msglvl)
        fig.subplots_adjust(bottom=0.14, top=0.98, left=0.08, right=0.97, wspace=0.05, hspace=0.05)
        fig.savefig(outdir + "/" + self.name + ".png", dpi=500)
        fig.clf()

class windFarmWakeStats(order):
    duration = timeInterval
    height = 0
    windfarm = 0
    pairs = []

    def __init__(self, name, identValues, identNames):
        self.name = name
        self.folder = wakeFolder()
        for identifier in identNames:
            k = identNames.index(identifier)
            if identifier == durIdentifier():
                self.duration = identValues[k]
            elif identifier == heightIdentifier():
                self.height = identValues[k]
            elif identifier == refWindFarmName():
                self.windfarm = int(identValues[k]) -1
            elif identifier == pairsIdentifier():
                self.pairs = identValues[k]

    def calcPairs(self, data, globe, levels, tag):
        if len(self.pairs) == 0: return
        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)
        outfile = open(outdir + "/Pairs(" + tag + ")-" + self.name + ".csv", "w")
        dataholder = zeros([len(levels), len(self.pairs)])
        header = "Level,"
        for pair in self.pairs:
            splitdata = splitline(pair, ",")
            caseLabel = splitdata[0]
            refLabel = splitdata[1]
            header += caseLabel + "/" + refLabel + ","
            caseindex = globe.labels.index(caseLabel) 
            refindex = globe.labels.index(refLabel) 
            for kl in range(len(levels)):
                caseValue = data[caseindex, kl]
                refValue = data[refindex, kl]
                dataholder[kl,self.pairs.index(pair)] = (caseValue - refValue) / refValue * 100

        outfile.write(header+ "\n")
        for kl in range(len(levels)):
            line = str(levels[kl]) + ","
            for ipair in range(len(self.pairs)): line += str(dataholder[kl,ipair]) + ","
            outfile.write(line + "\n")
        outfile.close()
          
    def execute(self, globe):
        msglvl = 0
        newline()
        message("Executing " + self.name, msglvl)
        msglvl += 1
        
        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)

        # Style
        nrows = 2
        ncols = 1
        barwidth = 0.25
        labelsize = 7
        ticksize = 6
        legendsize = 6
        levelcolors = ["darkred", "darkgoldenrod", "tab:blue"]
        areatags = ["$\mathrm{\geq 20\%}$", "$\mathrm{15\% - 20\%}$", "$\mathrm{10\% - 15\%}$"]
        disttags = ["$\mathrm{20\%}$", "$\mathrm{15\%}$", "$\mathrm{10\%}$"]
        levels = [20,15,10]

        outfile = open(outdir + "/" + self.name + "-WakeStats.csv", 'w')
        csvheader = "info,"
        for label in globe.labels: 
            if label == globe.ntLabel: continue 
            csvheader += label + ","
        outfile.write(csvheader + "\n")

        farmarea = windFarmsArea(self.windfarm)
        wakebnds = wakebox(self.windfarm)
        _, wakecenter = processbox(wakebnds)
        fc = farmcenter(self.windfarm)
        normals = getNorms(wakebnds, wakecenter)
        minlon = min(wakebnds[:,0])
        maxlon = max(wakebnds[:,0]) 
        minlat = min(wakebnds[:,1])
        maxlat = max(wakebnds[:,1]) 
                
        fig, axs = createPanel(ncols, nrows, 2, 5.5, 3.7)
        axs = axs.flatten()
        areahandles = []
        disthandles= []

        iref = globe.labels.index(globe.ntLabel)
        files = selectTimeFrames(select(listFiles(globe.dirs[iref]), globe.bases[iref]), self.duration)
        # to make calculations more accurate, we calculate wake per time frame
        # hence, all simulations are assumed to have the same time frames
        timeframes = [splitline(file, "/")[-1] for file in files]
        # since all cases will be compared to this refernece case, 
        # it is assumed that all cases have the same grid
        lats, lons, z, _ = loadGrid(files[0])
        dx = extractDirectVariable(files[0], dxname()) / 1000

        # to save time, all cases are assumed to have the same grid.
        # hence, we store the points inside the wind-farm's wake box
        # and loop over them for each WFP. 
        # If grids are different, these points will need to be re-calculated
        # for each variation of a grid.
        stored_i = []
        stored_j = []
        for j in range(lons.shape[0]):
            for i in range(lons.shape[1]):
                if not ((minlon <= lons[j,i] <= maxlon) and (minlat <= lats[j,i] <= maxlat)): continue
                if not dual_inbox(array([lons[j,i],lats[j,i]]), wakebnds, normals): continue
                stored_i.append(i)
                stored_j.append(j)

        labels = []
        dists = zeros([len(globe.dirs)-1, len(levels)])
        areas = zeros([len(globe.dirs)-1, len(levels)])
        for frame in timeframes:
            message("Time frame: " + frame, msglvl)
            msglvl += 1
            iframe = timeframes.index(frame)

            maxdists = zeros([len(globe.dirs)-1, len(levels)])
            ncells = zeros([len(globe.dirs)-1, len(levels)])

            # get the reference simulation
            message("Reference field", msglvl)
            iref = globe.labels.index(globe.ntLabel)
            if self.height == 10:
                refField = fieldExtract(globe.dirs[iref] + "/" + frame, speed10Name())
            else:
                refField = to_np(interplevel(fieldExtract(globe.dirs[iref] + "/" + frame, speedName()), z, self.height))

            position = -1
            for dir in globe.dirs:
                idir = globe.dirs.index(dir)
                label = globe.labels[idir]
                if label == globe.ntLabel: continue
                message(label, msglvl)

                if iframe == 0: labels.append(label)
                position += 1
                if self.height == 10:
                    field = fieldExtract(dir + "/" + frame, speed10Name())
                else:
                    field = to_np(interplevel(fieldExtract(dir + "/" + frame, speedName()), z, self.height))
                    
                wake = field_diff_percent2(field, refField)

                for kstore in range(len(stored_j)):
                    i = stored_i[kstore]
                    j = stored_j[kstore]
                    if wake[j,i] < min(levels): continue
                    for kl in range(len(levels)):
                        if wake[j,i] >= levels[kl]:
                            d = dist(fc,array([lons[j,i],lats[j,i]]))
                            if d > maxdists[position, kl]: maxdists[position, kl] = d
                            ncells[position, kl] += 1
                            break

                areas[position] += ncells[position,:] * dx * dx / float(farmarea)
                dists[position] += degToKm(maxdists[position,:])
            msglvl -= 1

        avgAreas = array([divide(areas[k], len(timeframes)) for k in range(len(globe.labels)-1)])
        avgDists = array([divide(dists[k], len(timeframes)) for k in range(len(globe.labels)-1)])
        
        for kl in range(len(levels)):
            line1 = "dist(" + str(levels[kl]) + "%),"
            line2 = "area(" + str(levels[kl]) + "%),"
            for km in range(len(globe.labels)):
                if globe.labels[km] == globe.ntLabel: continue
                line1 += str(avgDists[km, kl]) + ","
                line2 += str(avgAreas[km, kl]) + ","
            outfile.write(line1 + "\n")
            outfile.write(line2 + "\n")
        outfile.close()

        message("Evaluating pairs", msglvl)
        self.calcPairs(avgAreas, globe, levels, "NWE")
        self.calcPairs(avgDists, globe, levels, "Lw")

        message("Creating figure", msglvl)
        xaxis = [1 + 4*barwidth*n for n in range(len(labels))]
        data = []
        for j in range(2):
            if j == 0:
                data = array(avgAreas)
            elif j == 1:
                data = array(avgDists)

            ax = axs[j]
            offset = zeros(len(labels))

            for i in range(len(levels)):
                localx = [x + (i-1)*barwidth for x in xaxis]
                if j == 0: 
                    hatch_ = ' '
                    tag = areatags[i]
                    ec = "none"
                elif j == 1:
                    if i == 0:
                        hatch_ = '//'
                    elif i == 1:
                        hatch_ = '//'
                    elif i == 2:
                        hatch_ = '//'
                    ec = "black"
                    tag = disttags[i]
                bar = ax.bar(localx, data[:,i], barwidth, bottom=offset, color=levelcolors[i], label = tag, hatch=hatch_, edgecolor=ec)
                
                if j == 0: areahandles.append(bar)
                elif j == 1: disthandles.append(bar)

            ax.set_xticks(xaxis)
            if j == 1:
                ax.set_xticklabels(labels, fontsize = labelsize, rotation = 45)
            else:
                ax.xaxis.set_ticklabels([])

            if j == 0:
                ax.set_ylabel("NWE [-]", fontsize = labelsize)
                ax.set_ylim([0,3.1])
                ax.set_yticks([0,0.5,1,1.5,2,2.5,3])
                ax.set_xticks([])
            elif j == 1:
                ax.set_ylabel("$\mathrm{\mathcal{L}_w}$ [km]", fontsize = labelsize)
                ax.set_ylim([0,130])
                ax.set_yticks([0,20,40,60,80,100,120])
        
            ax.tick_params(axis='both', labelsize=ticksize)

            if j == 0:
                hand = areahandles
            elif j == 1:
                hand = disthandles
            ax.legend(handles = hand, loc='upper center', ncol=3, prop={'size': legendsize})

        message("Saving figure", msglvl)
        fig.subplots_adjust(bottom=0.12, top=0.97, left=0.1, right=0.97, wspace=0.05, hspace=0.05)
        fig.savefig(outdir + "/" + self.name + ".pdf")
        fig.clf()

def wakebox(farm):
    if farm == 0:
        return array([
        [7.011409732004955, 54.12819548872181],
        [7.092151937405591, 53.98947368421052],
        [7.9804180315348185, 54.290601503759405],
        [7.916529472665864, 54.40902255639098]
        ])
    else:
        return float("Nan")

def windFarmsArea(farm):
    if farm == 0:
        return 108.34 # km2
    else:
        return float("Nan")

def processbox(bnds):
    center = array([sum(bnds[:,0]), sum(bnds[:,1])]) / float(bnds.shape[0])
    area = 0
    for k in range(bnds.shape[0]):
        kp1 = k+1
        if k == bnds.shape[0]-1: kp1 = 0
        area += triarea(center, array(bnds[k,:]), array(bnds[kp1,:]))
    return area, center

def triarea(p1,p2,p3):
    return 0.5*linalg.norm(cross(p2-p1,p3-p2))

def farmcenter(farm):
    if farm == 0:
        return array([6.9772059946706495, 54.036842105263155])
    else:
        return float("Nan")

def getNorms(bnds, center):
    normals = bnds.copy()
    for k in range(bnds.shape[0]):
        p1 = array(bnds[k,:])
        if k == bnds.shape[0]-1:
            p2 = array(bnds[0,:])
        else:
            p2 = array(bnds[k+1,:])
        v = p2 - p1
        vc = (p1+p2)/2.0
        n = array([-v[1],v[0]])
        if dot(vc-center, n) < 0:
            normals[k,:] = -n / linalg.norm(n)
        else:
            normals[k,:] = n / linalg.norm(n)
    return normals

def dist(p1,p2):
    return linalg.norm(p2-p1)

def dual_inbox(p, bnds, normals):
    mindist =10000
    iclose = 0
    for k in range(bnds.shape[0]):
        d = dist(p, array(bnds[k,:]))
        if d < mindist: 
            mindist = d
            iclose = k

    if iclose == bnds.shape[0]-1:
        va = 0
    else:
        va = iclose + 1

    if iclose  == 0:
        vb = bnds.shape[0]-1
    else:
        vb = iclose - 1

    a1 = array(bnds[vb,:]) - array(bnds[iclose,:])
    a2 = array(bnds[va,:]) - array(bnds[iclose,:])
    r1 = p - (array(bnds[vb,:]) + array(bnds[iclose,:]))/2.0
    r2 = p - (array(bnds[va,:]) + array(bnds[iclose,:]))/2.0

    c = a1[0]*a2[1] - a1[1]*a2[0]
    d1 = dot(r1,array(normals[vb,:]))
    d2 = dot(r2,array(normals[va,:]))

    if d1 > 0: s1 = 1
    elif d1 < 0: s1 = -1
    else: s1 = 0  

    if d2> 0: s2 = 1
    elif d2 < 0: s2 = -1
    else: s2 = 0

    if s1 == 0 or s2 == 0: return True

    if c > 0:
        if s1 == 1 and s2 == 1:
            return False
        else:
            return True
    elif c == 0:
        if s1 == 1: 
            return False
        else:
            return True
    else:
        if s1 == -1 and s2 == -1:
            return True
        else:
            return False
