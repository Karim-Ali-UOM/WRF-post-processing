from src.domains import domainPlot
from src.sar import sarPlot
from src.temporal import temporalPlot
from src.transect import TransectFlightsPlot
from src.turbulenceMod import qbudget
from src.utilities import domainPlotName, message, qbudgetName, sarPlotName, temporalPlotName, transectPlotName, \
    verticalProfilesName, vgradPlotName, wakeContourName, wakeStatsNames, wfcapacityName
from src.powerMod import capacityFactorPlot
from src.verticalMod import verticalProfilesPlot
from src.vgrad import vgradPlot
from src.wake import wakeContourPlot, windFarmWakeStats

def register(name, orderType, identValues, identNames, msglvl):
    message("Registering order " + name, msglvl)
    
    if orderType == wfcapacityName():
        return capacityFactorPlot(name, identValues, identNames)
    elif orderType == qbudgetName():
        return qbudget(name, identValues, identNames)
    elif orderType == vgradPlotName():
        return vgradPlot(name, identValues, identNames)
    elif orderType == domainPlotName():
        return domainPlot(name, identValues, identNames)
    elif orderType == sarPlotName():
        return sarPlot(name, identValues, identNames)
    elif orderType == wakeContourName():
        return wakeContourPlot(name, identValues, identNames)
    elif orderType == verticalProfilesName():
        return verticalProfilesPlot (name, identValues, identNames)
    elif orderType == temporalPlotName():
        return temporalPlot (name, identValues, identNames)
    elif orderType == transectPlotName():
        return TransectFlightsPlot(name, identValues, identNames)
    elif orderType == wakeStatsNames():
        return windFarmWakeStats(name, identValues, identNames)
