from src.orders import order
from src.timeMod import selectNearestTimeFrames, timeInstant
from src.utilities import assistFilesIdentifier, capsIdentifier, createCrsPanel, fieldExtract, listFiles, loadGrid, markWindFarmCells, \
    maskField, maskName, message, newline, prepareloc, rowCol, select, speed10Name, splitline, timeInstantIdentifier, toLetters, zoomIdentifier
from math import ceil
from numpy import linspace, arange, min, max

def sarFolderName():
    return "SAR"

class sarPlot(order):
    sarfiles = []
    instant = timeInstant
    minWindSpeed = 0
    maxWindSpeed = 0
    zoom = []

    def __init__(self, name, identValues, identNames):
        self.name = name
        self.folder = sarFolderName()
        for identifier in identNames:
            k = identNames.index(identifier)
            if identifier == assistFilesIdentifier():
                self.sarfiles = identValues[k]
            elif identifier == timeInstantIdentifier():
                self.instant = identValues[k]
            elif identifier == zoomIdentifier():
                self.zoom = identValues[k]
            elif identifier == capsIdentifier():
                caps = identValues[k]
        self.minWindSpeed = min(caps)
        self.maxWindSpeed = max(caps)

    def execute(self, globe):

        import cartopy.crs as crs
        from netCDF4 import Dataset
        from wrf import to_np
        from matplotlib.cm import get_cmap
        from cartopy.mpl.ticker import (LongitudeFormatter, LatitudeFormatter)

        msglvl = 0
        newline()
        message("Executing " + self.name, msglvl)
        msglvl += 1

        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)

        # style
        nconts = 300
        plotsPerRow = 3
        nRows = ceil((len(globe.labels) + 1) / plotsPerRow)

        crss = crs.PlateCarree()
        fig, axGlobal = createCrsPanel(plotsPerRow, nRows, len(globe.labels) + 1, 5.5, 3.5, crss)
        axGlobal = axGlobal.flatten()

        # find wind-farms footprint using the power field in any file
        files = select(listFiles(globe.dirs[0]), globe.bases[0])
        lats, lons, wfCells = markWindFarmCells(files[0])

        # Plot SAR images
        for file in self.sarfiles:
            ifile = self.sarfiles.index(file)
            ax = axGlobal[0]

            sar1 = Dataset(file)
            sar_speed = to_np(sar1.variables["sar_wind"][:][:])
            sar_lon = to_np(sar1.variables["longitude"][:][:])
            sar_lat = to_np(sar1.variables["latitude"][:][:])
            mask = to_np(sar1.variables["mask"][:][:])
            for j in range(sar_lon.shape[0]):
                for i in range(sar_lon.shape[1]):
                    if mask[j,i] >= 0: sar_speed[j,i] = float("Nan")

            if ifile == 0:
                ax.contour(lons, lats, wfCells, 1, vmin = 0.5, vmax = 2, transform=crss, colors = "k",
                    zorder=20, alpha = 0.75, linewidths = 0.1) 
                   
            im = ax.contourf(sar_lon, sar_lat, sar_speed, linspace(self.minWindSpeed,self.maxWindSpeed,nconts),
                cmap=get_cmap("viridis"), transform=crss, zorder = 1)
        
        for dir in globe.dirs:
            idir = globe.dirs.index(dir)
            ax = axGlobal[idir + 1]
            files = select(listFiles(dir), globe.bases[idir])
            selectedTimeFiles, weights = selectNearestTimeFrames(files, self.instant)
            mask = fieldExtract(selectedTimeFiles[0], maskName())
            sol = 0
            localLats, localLons, _, _ = loadGrid(selectedTimeFiles[0])
            for file in selectedTimeFiles:
                sol += fieldExtract(file, speed10Name()) * weights[selectedTimeFiles.index(file)]
            # wind farm footprint
            ax.contour(lons, lats, wfCells, 1, vmin = 0.5, vmax = 2, transform=crss, colors = "k",
                    zorder=20, alpha = 0.75, linewidths = 0.1)
            sol = maskField(sol, mask, 2)
            # wind speed
            im = ax.contourf(localLons, localLats, sol, linspace(self.minWindSpeed,self.maxWindSpeed,nconts),
                cmap=get_cmap("viridis"), transform=crss, zorder = 1)

        fig.subplots_adjust(bottom=0.14, top=0.97, left=0.08, right=0.97, wspace=0.05, hspace=0.07)
        iax = -1
        for ax in axGlobal:
            iax += 1
            row, col = rowCol(iax, plotsPerRow)

            if iax == 0: name = "SAR"
            else: name = globe.labels[iax-1]
            ax.text(5.58, 54.77, "(" + toLetters(iax) + ") " + name, fontsize = 6)
            ax.set_extent(self.zoom)        
            ax.yaxis.tick_left()

            lonStep = 0.5
            latStep = 0.3 
            ax.set_xticks([round(q,1) for q in arange(self.zoom[0]+lonStep, self.zoom[1], lonStep)], crs=crss)
            ax.set_yticks([round(q,1) for q in arange(self.zoom[2]+latStep, self.zoom[3] - latStep, latStep)], crs=crss)
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            ax.yaxis.set_major_formatter(lat_formatter)
            ax.tick_params(axis='both', which='major', labelsize=6)
            ax.grid(False)
            if row != nRows-1: ax.xaxis.set_ticklabels([])
            if col != 0: ax.yaxis.set_ticklabels([])

        cb_ax = fig.add_axes([0.2, 0.06, 0.6, 0.02])
        cbar = fig.colorbar(im, cax=cb_ax, orientation='horizontal', format='%.0f')
        cbar.ax.tick_params(labelsize=6) 
        ticks = linspace(self.minWindSpeed,self.maxWindSpeed,10)
        cbar.set_ticks(ticks)
        cbar.ax.set_title("U $\mathrm{[m\,s^{-1}]}$", fontsize = 6, x=1.1, y=0.1, pad=0)
        
        fig.savefig(outdir + "/" + self.name + ".png", dpi=500)
        fig.clf()