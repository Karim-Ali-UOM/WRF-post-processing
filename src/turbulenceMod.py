from src.orders import order
from src.timeMod import selectTimeFrames
from src.utilities import createPanel, durIdentifier, listFiles, locIdentifier, message, metersToIndex, newline, prepareloc, q2BudgetColor, \
    q2BudgetFields, q2BudgetTag, qtName, segmentToMeters, segmentsIdentifier, toLetters, select, qdName
from src.verticalMod import vertical
from numpy import zeros, arange, trapz, sum, array, max, min
from math import ceil

def turbulenceFolder():
    return "turbulence"

class qbudget(order):
    segments = []
    verticalProfile = vertical 
    fields = []

    def __init__(self, name, identValues, identNames):
        self.name = name
        self.folder = turbulenceFolder()
        self.fields = q2BudgetFields()
        for identifier in identNames:
            k = identNames.index(identifier)
            if identifier == locIdentifier():
                location = identValues[k]
            elif identifier == durIdentifier():
                duration = identValues[k]
            elif identifier == segmentsIdentifier():
                self.segments = identValues[k]
        self.verticalProfile = vertical(name, location, duration)
        
    def execute(self, globe):
        msglvl = 0
        newline()
        message("Executing " + self.name, msglvl)
        msglvl += 1
        
        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)
        # Style
        barwidth = 0.5
        ncols = 2
        labelsize = 6
        ticksize = 6
        legendsize = 7

        nrows = ceil(len(self.segments)/ncols)
        xaxis = arange(len(globe.labels)) + (1-barwidth)/2.0

        # Obtain limiting indices per segment and per directory
        # Grid is assumed not to be changing with time
        limIndices = zeros([len(globe.labels), len(self.segments), 2])

        # Convert segments to meters
        limMeters = zeros([len(self.segments), 2])
        includeLimits = zeros([len(self.segments), 2])
        message("Converting segments to meters", msglvl)
        for segment in self.segments:
            isegment = self.segments.index(segment)
            limMeters[isegment,0], limMeters[isegment,1], includeLimits[isegment,0], \
                includeLimits[isegment,1] = segmentToMeters(segment, globe)

        # Identify the closest grid cell
        message("Finding the closest grid-cell for all cases", msglvl)
        self.verticalProfile.closestGridCell(globe)

        for interval in self.intervals:
            message("Interval " + str(self.intervals.index(interval)), msglvl)
            msglvl += 1
            # Every interval have its own plot
            legendHandles = []
            kinterval = self.intervals.index(interval)
            avgFieldsHolder = zeros([len(self.segments), len(globe.dirs), len(self.fields)])
            dh = limMeters[-1, 1] - limMeters[0, 0]

            for dir in globe.dirs:
                idir = globe.dirs.index(dir)
                message(globe.labels[idir], msglvl)
                msglvl += 1
                files = selectTimeFrames(select(listFiles(dir), globe.bases[idir]), interval)

                # this is done for the first interval only
                if kinterval == 0: limIndices[idir, :, :] = metersToIndex(self.verticalProfile.locLat,
                    self.verticalProfile.locLon, files[0], self.verticalProfile.ic[idir], 
                    self.verticalProfile.jc[idir], limMeters, includeLimits)

                for field in self.fields:
                    ifield = self.fields.index(field)
                    zholder, dataholder = self.verticalProfile.execute(globe, field, interval, dir, None)
                    gain = 1
                    if field == qdName(): gain = -1
                    elif field == qtName() and globe.labels[idir] == globe.ntLabel: gain = 0
                    for segment in self.segments:
                        isegment = self.segments.index(segment)
                        i1 = int(limIndices[idir,isegment,0])
                        i2 = int(limIndices[idir,isegment,1])
                        avgFieldsHolder[isegment, idir, ifield] = gain* trapz(dataholder[i1:i2], zholder[i1:i2]) / dh
                msglvl -= 1
            msglvl -= 1

            # scale q2
            avgFieldsHolder = avgFieldsHolder * 1000

            # all vertical averages are calculated
            message("Formatting figure", msglvl)
            fig, axs = createPanel(ncols, nrows, len(self.segments), 6, 5)
            axs = axs.flatten()
            for segment in self.segments:
                isegment = self.segments.index(segment)
                ax = axs[isegment]
                row = isegment // ncols
                col = isegment - row * ncols
                yoffset = zeros(len(globe.dirs))
                yoffsetup = zeros(len(globe.dirs))
                yoffsetdown = zeros(len(globe.dirs))
                qtot = [sum(array(avgFieldsHolder[isegment, idir, :])) for idir in range(len(globe.dirs))]
                
                for ifield in range(len(self.fields)):
                    if ifield > 0:
                        for idir in range(len(globe.dirs)):
                            if avgFieldsHolder[isegment, idir, ifield] >= 0:
                                yoffset[idir] = yoffsetup[idir]
                            else:
                                yoffset[idir] = yoffsetdown[idir]

                    bar = ax.bar(xaxis, avgFieldsHolder[isegment, :, ifield], barwidth, bottom=yoffset,
                        color=q2BudgetColor(self.fields[ifield]), label = q2BudgetTag(self.fields[ifield]))

                    if isegment == 0: legendHandles.append(bar)

                    for idir in range(len(globe.dirs)):
                        if avgFieldsHolder[isegment, idir, ifield] >= 0:
                            yoffsetup[idir] += avgFieldsHolder[isegment, idir, ifield]
                        else:
                            yoffsetdown[idir] += avgFieldsHolder[isegment, idir, ifield]

                ax.set_xticks(xaxis)
                ax.set_xticklabels(globe.labels, fontsize = labelsize, rotation = 45)
                if row != nrows-1:
                    ax.xaxis.set_ticklabels([])
                ax.set_ylabel("$\\times 10^{-3}\, \mathrm{[cm^2\,s^{-2}]}$", fontsize = labelsize)
                ax.set_ylim([1.05 * min(yoffsetdown), 1.15 * max(yoffsetup)])

                title = '%g' %limMeters[isegment, 0] + " m $-$ " + '%g' %limMeters[isegment, 1] + " m"
                
                ax.text(0.6, 0.93, "(" + toLetters(isegment) + ") " + title, fontsize = labelsize, transform=ax.transAxes)        
                ax.tick_params(axis='both', labelsize=ticksize)
                ax.axhline(y=0, color='grey', linestyle='--', linewidth = 0.5)
                ax.plot(xaxis, qtot, color="k", marker='o', markersize = 2)

            message("Saving figure", msglvl)
            fig.legend(handles = legendHandles, loc = [0.32, 0.94], ncol=5, prop={'size': legendsize})
            fig.subplots_adjust(bottom=0.1, top=0.92, left=0.1, right=0.95, wspace=0.3, hspace=0.225)
            fig.savefig(outdir + "/" + self.name + "-I(" + interval.name + ").pdf")
            fig.clf()
