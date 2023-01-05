from src.orders import order
from src.utilities import durIdentifier, locIdentifier, message, nearest, newline, prepareloc, refWindFarmName, toLetters, uname, vname, createFigure
from src.verticalMod import vertical
import matplotlib.pyplot as plt
from numpy import gradient, array

def vgradFolder():
    return "SpeedVerticalGradient"

class vgradPlot(order):
    verticalProfile = vertical
    wfindex = 0
    folder = ""

    def __init__(self, name, identValues, identNames):
        self.folder = vgradFolder()
        self.name = name
        for identifier in identNames:
            k = identNames.index(identifier)
            if identifier == locIdentifier():
                location = identValues[k]
            elif identifier == durIdentifier():
                duration = identValues[k]
            elif identifier == refWindFarmName():
                self.wfindex = identValues[k]
        self.verticalProfile = vertical(name, location, duration)
        
    def execute(self, globe):
        msglvl = 0
        newline()
        message("Executing " + self.name, msglvl)
        msglvl += 1

        outdir = globe.outdir + "/" + self.folder
        prepareloc(outdir)

        #style
        linewidth = 0.9
        labelsize = 7
        ticksize = 6
        hubwidth = 0.5
        legendsize = 7

        colors = {"A100":"tab:red", "A90":"tab:red", "A80":"tab:red", "F100":"tab:blue", "F25":"tab:blue", "F100-NE":"tab:blue", 
                "V100":"tab:green", "V120":"tab:green", "V80":"tab:green", "PAN":"tan", "R100":"purple", "R25":"purple",
                "NT":"black", "EXP":"black"}
        styles = {"EXP":"--", "A100":"-", "A90":"--", "A80":":", "F100":"-", "F25":"--", "F100-NE":":", 
                "V100":"-", "V120":"--", "V80":":", "PAN":"-", "P25":"--", "R100":"-", "R25":"--", "NT":"-"}
        tags = ["$\partial u/\partial z$", "$\partial v/\partial z$", "$\mathcal{U}$"]
        
        legendHandles = []
        arrowloc = [230, 170, 75]
        arrowloc2 = [23, 23, 23]
        fig = createFigure(6,7)
        axu1 = plt.subplot2grid((10,3),(0,0),rowspan=7)
        axu2 = plt.subplot2grid((10,3),(8,0),rowspan=2)
        axv1 = plt.subplot2grid((10,3),(0,1),rowspan=7)
        axv2 = plt.subplot2grid((10,3),(8,1),rowspan=2)
        axS1 = plt.subplot2grid((10,3),(0,2),rowspan=7)
        axS2 = plt.subplot2grid((10,3),(8,2),rowspan=2)
        axu1t = axu1.twiny()
        axu2t = axu2.twiny()
        axv1t = axv1.twiny()
        axv2t = axv2.twiny()
        axS1t = axS1.twiny()
        axS2t = axS2.twiny()
        axs1 = [axu1, axv1, axS1]
        axs1t = [axu1t, axv1t, axS1t]
        axs2 = [axu2, axv2, axS2]
        axs2t = [axu2t, axv2t, axS2t]

        ugradHolder = [[] for k in globe.dirs]
        vgradHolder = [[] for k in globe.dirs]
        zholder = [[] for k in globe.dirs]
        
        message("Extracting vertical profiles", msglvl)
        for dir in globe.dirs:
            idir = globe.dirs.index(dir)
            message(globe.labels[idir], msglvl+1)
            z, u = self.verticalProfile.execute(globe, uname(), self.verticalProfile.intervals[0], dir, None)
            _, v = self.verticalProfile.execute(globe, vname(), self.verticalProfile.intervals[0], dir, None)

            zholder[idir] = z
            ugradHolder[idir] = array(gradient(u,z)) * 60 # convert to 1/min
            vgradHolder[idir] = array(gradient(v,z)) * 60 # convert to 1/min

        iref = globe.labels.index(globe.ntLabel)
        ugradRef = ugradHolder[iref]
        vgradRef = vgradHolder[iref]

        message("Plotting results", msglvl)
        for label in globe.labels:
            ilabel = globe.labels.index(label)
            ugrad = ugradHolder[ilabel]
            vgrad = vgradHolder[ilabel]
            z = zholder[ilabel]
            if label == globe.ntLabel:
                curve, = axu1t.plot(ugrad, z, color=colors[label], 
                    linestyle = styles[label], label=label, linewidth = linewidth)
                legendHandles.append(curve)
                axv1t.plot(vgrad, z, color=colors[label], linestyle = styles[label],
                    label=label, linewidth = linewidth)
                axS1t.plot((ugrad**2+vgrad**2), z, color=colors[label], linestyle = styles[label],
                    label=label, linewidth = linewidth)

                axu2t.plot(ugrad, z, color=colors[label], linestyle=styles[label], label=label,
                    linewidth = linewidth)
                axv2t.plot(vgrad, z, color=colors[label], linestyle=styles[label], label=label,
                    linewidth = linewidth)
                axS2t.plot((ugrad**2+vgrad**2), z, color=colors[label], linestyle=styles[label],
                    label=label, linewidth = linewidth)
            else:
                curve, = axu1.plot((ugrad-ugradRef), z, color=colors[label], linestyle=styles[label],
                    label=label, linewidth = linewidth)
                legendHandles.append(curve)
                axv1.plot((vgrad-vgradRef), z, color=colors[label], linestyle=styles[label],
                    label=label, linewidth = linewidth)
                axS1.plot((ugrad**2 + vgrad**2)-(ugradRef**2 + vgradRef**2), z, color=colors[label],
                    linestyle = styles[label], label=label, linewidth = linewidth)

                axu2.plot((ugrad-ugradRef), z, color=colors[label], linestyle=styles[label],
                    label=label, linewidth = linewidth)
                axv2.plot((vgrad-vgradRef), z, color=colors[label], linestyle=styles[label],
                    label=label, linewidth = linewidth)
                axS2.plot((ugrad**2 + vgrad**2)-(ugradRef**2 + vgradRef**2), z, 
                    color=colors[label], linestyle=styles[label], label=label, linewidth = linewidth)

        # Ref turbine hub and diameter
        turbType = globe.windFarmTurbineType(int(self.wfindex))
        hubheight = globe.turbTypeHH(turbType)
        diameter = globe.turbineTypeDiameter(turbType)

        message("Formatting figure", msglvl)
        # top plots
        # bottom axes
        iax = -1
        for ax in axs1:
            iax += 1
            ax.tick_params(axis='both', labelsize=ticksize)
            ax.text(0.88, 0.95, "(" + toLetters(iax) + ")", fontsize = labelsize, transform=ax.transAxes)
            ax.set_xlabel("$\mathrm{\Delta} ($" + tags[iax] + ") $[\mathrm{{min}^{-1}}]$", fontsize = labelsize)
            if iax == 0:
                ax.set_ylabel("z [m]", fontsize = labelsize)
                ax.set_xlim([-1,1.2])
                ax.set_xticks([-0.8,-0.4,0,0.4,0.8])
            elif iax == 1:
                ax.get_yaxis().set_visible(False)
                ax.set_xlim([-0.8,1.2])
                ax.set_xticks([-0.8,-0.4,0,0.4,0.8,1.2])
            elif iax == 2:
                ax.get_yaxis().set_visible(False)
                ax.set_xlim([-2.3,4])
                ax.set_xticks([-2,-1,0,1,2,3,4])
                ax.set_xlabel("$\mathrm{\Delta}$" + tags[iax] + " $[\mathrm{{min}^{-2}}]$", fontsize = labelsize)
            ax.set_ylim(bottom = 50, top = 400)
            ax.axhline(y=hubheight, color='grey', linestyle='--', linewidth = hubwidth)
            ax.axhline(y=hubheight+diameter/2.0, color='grey', linestyle='--', linewidth = hubwidth)
            ax.axhline(y=hubheight-diameter/2.0, color='grey', linestyle='--', linewidth = hubwidth)    
            ax.axvline(x=0, color='k', linestyle='--', linewidth = hubwidth) 

        # top plots
        # top axes
        iax = -1
        for ax in axs1t:
            iax += 1
            ax.set_xlabel(tags[iax] + " [$\mathrm{{min}^{-1}}$]", fontsize = labelsize)
            ax.tick_params(axis='both', labelsize=ticksize)
            #ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            if iax == 0:
                ax.set_xlim([-0.2,1.6])
                arrowxloc = ugradRef[nearest(zholder[iref],arrowloc[iax])]
            elif iax == 1:
                ax.set_xlim([-1.2,0.6])
                arrowxloc = vgradRef[nearest(zholder[iref],arrowloc[iax])]
            elif iax == 2:
                ax.set_xlim([0,2.5])
                arrowxloc = ugradRef[nearest(zholder[iref],arrowloc[iax])]**2 + vgradRef[nearest(zholder[iref],arrowloc[iax])]**2
                ax.set_xlabel(tags[iax] + " [$\mathrm{{min}^{-2}}$]", fontsize = labelsize)
            ax.annotate("", xy=(arrowxloc, arrowloc[iax] + 40), xycoords='data', 
                xytext=(arrowxloc, arrowloc[iax] + 5), textcoords='data', 
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3"), 
                color="k")

        # bottom plots
        # bottom axes
        iax = -1
        for ax in axs2:
            iax += 1
            ax.tick_params(axis='both', labelsize=ticksize)
            ax.text(0.88, 0.85, "(" + toLetters(iax+3) + ")", fontsize = labelsize, transform=ax.transAxes)
            ax.set_xlabel("$\mathrm{\Delta} ($" + tags[iax] + ") $[\mathrm{{min}^{-1}}]$", fontsize = labelsize)
            if iax == 0:
                ax.set_ylabel("z [m]", fontsize = labelsize)
                ax.set_xlim([-5,1])
            elif iax == 1:
                ax.get_yaxis().set_visible(False)
                ax.set_xlim([-4,1])
            elif iax == 2:
                ax.get_yaxis().set_visible(False)
                ax.set_xlim([-120,25])
                ax.set_xlabel("$\mathrm{\Delta}$" + tags[iax] + " $[\mathrm{{min}^{-2}}]$", fontsize = labelsize)
            ax.set_ylim(bottom = 0, top = 50)
            ax.axhline(y=hubheight-diameter/2.0, color='grey', linestyle='--', linewidth = hubwidth)    
            ax.axvline(x=0, color='k', linestyle='--', linewidth = hubwidth) 

        # bottom plots
        # top axes
        iax = -1
        for ax in axs2t:
            iax += 1
            ax.set_xlabel(tags[iax] + " [$\mathrm{{min}^{-1}}$]", fontsize = labelsize)
            ax.tick_params(axis='both', labelsize=ticksize)
            #ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
            wght = (arrowloc2[iax]-zholder[iref][1])/(zholder[iref][2]-zholder[iref][1])
            if iax == 0:
                ax.set_xlim([1,11])
                arrowxloc = wght*ugradRef[2] + (1-wght) * ugradRef[1]
            elif iax == 1:
                ax.set_xlim([0, 8])
                arrowxloc = wght*vgradRef[2] + (1-wght) * vgradRef[1]
            elif iax == 2:
                ax.set_xlim([0,160])
                arrowxloc = wght*(ugradRef[2]**2 + vgradRef[2]**2) + (1-wght)*(ugradRef[1]**2 + vgradRef[1]**2)
                ax.set_xlabel(tags[iax] + " [$\mathrm{{min}^{-2}}$]", fontsize = labelsize)
            ax.annotate("", xy=(arrowxloc, arrowloc2[iax] + 18), xycoords='data', 
                xytext=(arrowxloc, arrowloc2[iax]+2), textcoords='data', 
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3"), 
                color="k")

        message("Saving figure", msglvl)
        fig.legend(handles = legendHandles, loc = [0.11, 0.9], ncol=7, prop={'size': legendsize})
        fig.subplots_adjust(bottom=0.075, top=0.845, left=0.1, right=0.95, wspace=0.1, hspace=0.4)
        fig.savefig(outdir + "/" + self.name + ".pdf")
        fig.clf()