from glob import glob
from src.orders import order
from src.utilities import boundsIdentifier, createFigure, dronPathIdentifier, dxname, extractDirectVariable, filename, kmToDeg, loadDefaultColors, \
    loadGrid, message, newline, prepareloc, refFileName, splitline, windFarmNameLocation, printedWindFarmsNames, zoomIdentifier
from numpy import min, max

def domainPlotFolder():
    return "domains"

class domainPlot(order):
    bounds = []
    zoom = []
    domainFiles = []
    dronPathsFiles = []
    folder = str

    def __init__(self, name, identValues, identNames):

        self.folder = domainPlotFolder()
        self.name = name
        for identifier in identNames:
            k = identNames.index(identifier)
            if identifier == boundsIdentifier():
                self.bounds = identValues[k]
            elif identifier == zoomIdentifier():
                self.zoom = identValues[k]
            elif identifier == refFileName():
                try: self.domainFiles =  glob(identValues[k])
                except: self.domainFiles = []
            elif identifier == dronPathIdentifier():
                try: self.dronPathsFiles =  glob(identValues[k])
                except: self.dronPathsFiles = []

    def execute(self, globe):

        import cartopy.feature as cfeature
        from matplotlib.pyplot import axes, scatter, locator_params, text
        import cartopy.crs as crs
        from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
        from mpl_toolkits.axes_grid1.inset_locator import mark_inset

        msglvl = 0
        newline()
        message("Executing " + self.name, msglvl)
        msglvl += 1

        # Style
        bndLineWidth = 1
        bndColor = "red"
        fontSize = 10
        txtNudge = 50
        turbSize = 0.05
        nameFontSize = 10

        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)

        # load global data
        crss = crs.PlateCarree()
        fig = createFigure(12,6)
        ax = axes(projection=crss)
        colors = loadDefaultColors()
        if bndColor in colors: colors[colors.index(bndColor)] = "darkolivegreen"

        land = cfeature.NaturalEarthFeature(category="physical", scale="50m",
                                    facecolor=cfeature.COLORS['land'],
                                    name="land",edgecolor='face')
        water = cfeature.NaturalEarthFeature(category="physical", scale="50m",
                                    facecolor=cfeature.COLORS['water'],
                                    name="ocean", edgecolor='face')
                        
        border = cfeature.NaturalEarthFeature(category="cultural", scale="50m",
                                    facecolor="none",
                                    name="admin_0_countries", edgecolor='face')

        ax.set_extent(self.bounds)

        for file in self.domainFiles:
            message("Plotting domain boundaries from " + filename(file), msglvl)
            bndLats, bndLons, _, _ = loadGrid(file)
            ax.plot(bndLons[0,:], bndLats[0,:],linewidth = bndLineWidth, color = bndColor)
            ax.plot(bndLons[-1,:], bndLats[-1,:],linewidth = bndLineWidth, color = bndColor)
            ax.plot(bndLons[:,0], bndLats[:,0],linewidth = bndLineWidth, color = bndColor)
            ax.plot(bndLons[:,-1], bndLats[:,-1],linewidth = bndLineWidth, color = bndColor)
            cornerLon = bndLons[-1,0]
            cornerLat = bndLats[-1,0]
            dx = extractDirectVariable(file, dxname()) / 1000
            if dx - int(dx) == 0:
                txt = str(dx) + " km"
            else:
                txt = str(round(dx,2)) + " km"
            ax.text(cornerLon + kmToDeg(60), cornerLat - kmToDeg(90), 
                txt, fontsize = fontSize, zorder=100, weight='bold')

        message("Plotting turbines location on main axes", msglvl)
        for kt in range(len(globe.turbLatDict)):
            scatter(float(globe.turbLonDict[kt+1]), float(globe.turbLatDict[kt+1]), s = turbSize,
            c = colors[int(globe.turbineIndexToWindFarmIndex(kt+1))-1], zorder=2)

        message("Formatting main axes", msglvl)
        ax.add_feature(land, linewidth=.25, edgecolor="black")
        #ax.add_feature(water, linewidth=.25, edgecolor="black")   
        ax.add_feature(border, linewidth=.25, edgecolor="black")
        ax.coastlines('50m', linewidth=0.8)
        gl = ax.gridlines(draw_labels=True, crs=crss, color = "lightgray", alpha = 0.75)
        gl.xlabels_bottom = False
        gl.ylabels_left = False
        gl.xlines = True
        gl.ylines = True
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        message("Creating secondary axes for zoom-in", msglvl)
        ax2 = axes([0.126, 0.11, 0.4, 0.35], projection=crss)  
        ax2.set_xlim(self.zoom[0], self.zoom[1])
        ax2.set_ylim(self.zoom[2], self.zoom[3])
        locator_params(axis='y', nbins=4)
        locator_params(axis='x', nbins=4)
        mark_inset(ax, ax2, loc1=2, loc2=4, fc="none", ec='k', zorder=200)     
        
        message("Plotting turbines location on secondary axes", msglvl)
        for kt in range(len(globe.turbLatDict)):
                scatter(float(globe.turbLonDict[kt+1]), float(globe.turbLatDict[kt+1]), s = turbSize,
                c = colors[int(globe.turbineIndexToWindFarmIndex(kt+1))-1], zorder=2)

        message("Labelling turbines location on secondary axes", msglvl)
        for j in range(globe.numWindFarms):
            loc = windFarmNameLocation(globe.windFarmName(j+1))
            text(loc[0], loc[1], printedWindFarmsNames(globe.windFarmName(j+1)), fontsize = nameFontSize, zorder=10*(j+1))

        message("Formatting secondary axes", msglvl)
        ax2.add_feature(land, linewidth=.25, edgecolor="black")
        #ax2.add_feature(water, linewidth=.25, edgecolor="black")   
        ax2.add_feature(border, linewidth=.25, edgecolor="black")
        ax2.coastlines('50m', linewidth=0.8)
        gl = ax2.gridlines(draw_labels=True, color = "lightgray", alpha = 0.75)
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = True
        gl.ylines = True
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        message("Plotting drone paths on secondary axes", msglvl)
        for file in self.dronPathsFiles:
            drone = open(file, "r")
            latdrone = []
            londrone = []
            while True:
                try: 
                    line = splitline(drone.readline(), ",")
                except:
                    break
                if len(line) == 0: break
                if line[0] == "Time": continue

                latdrone.append(float(line[1]))
                londrone.append(float(line[2]))
            ax2.plot(londrone, latdrone, color="green", linestyle = "--", linewidth = 0.7)
        
        # add FINO-1
        message("Adding the FINO-1 location to secondary axes", msglvl)
        ax2.scatter(6.58764, 54.01486, s = 5, c = "tab:red", zorder=2)
        ax2.annotate("FINO-1",
                xy=(6.58764, 54.01486), xycoords='data',
                xytext=(6.605, 54.35), textcoords='data',
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad = 0.3"), color = "tab:red")
                
        message("Saving the figure", msglvl)
        fig.savefig(outdir + "/" + self.name + ".pdf")
        fig.clf()