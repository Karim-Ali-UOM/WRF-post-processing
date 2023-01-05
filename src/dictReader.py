from src.utilities import assistFilesIdentifier, boundsIdentifier, capsIdentifier, dronPathIdentifier,\
        durIdentifier, heightIdentifier, intervalsIdentifier, levelsIdentifier, locIdentifier, message, newline, \
        pairsIdentifier, refCaseKeyword, refFileName, refWindFarmName, segmentsIdentifier, splitline, dirsKeyword, \
        outkeyword, ordersKeyword, activeIdentifier, timeInstantIdentifier, typeIdentifier, \
        wfIndexIdentifier, turbLocKeyword, turbCoordIdentifier, turbCoordCode, \
        turbineTypesTablesIdentifier, turbTypeCode, wfnamesCode, windFarmsnamesIdentifier,\
        mainPacketKeyword, fataError, subpacketIdentifiers, zoomIdentifier
from src.register import register
from src.timeMod import timeInstant, timeInterval
from src.globe import globe

def extractLinesInLevel(lines, level):
    extractedLines = []
    openBrackets = 0
    for line in lines:
        if line == "(": openBrackets += 1
        elif line == ")": openBrackets -= 1
        elif openBrackets == level: extractedLines.append(line)

    return extractedLines

def refineLines(file):

    startedMainPacket = False
    foundMainPacketBracket = False
    openBrackets = 0
    lines = []
    dictfile = open(file,"r")
    alllines = dictfile.readlines()
    for currentline in alllines:
        try: line = splitline(currentline,"\n")
        except: continue
        ## check if it is an empty line
        if len(line) == 0 or line[0] in ["", "\n", " "] or line[0][0] == "#": continue    
        lines.append(line[0])
        if line[0] == mainPacketKeyword(): startedMainPacket = True
        if line[0] == "(" and startedMainPacket and not foundMainPacketBracket: foundMainPacketBracket=True
        if line[0] == "(": openBrackets += 1
        if line[0] == ")": openBrackets -= 1
        if foundMainPacketBracket and openBrackets == 0: break
    dictfile.close()
    return lines

def splitPacket(lines, keyword):
    openBrackets = 0
    foundKeyword = False
    foundKeywordBracket = False
    foundKeywordSecBracket = False
    recordNBackets = 0
    extractedLines = []
    leftoverLines = []
    for line in lines:
        if line == "(": openBrackets += 1
        if line == ")": openBrackets -= 1
        if line == keyword: 
            foundKeyword = True
            recordNBackets = openBrackets
            continue
        if line == "(" and foundKeyword and not foundKeywordBracket:
            foundKeywordBracket= True
            continue
        if foundKeyword and foundKeywordBracket and not foundKeywordSecBracket:
             if openBrackets > recordNBackets: extractedLines.append(line)
             else: foundKeywordSecBracket = True
        else:
            leftoverLines.append(line)
    return extractedLines, leftoverLines

def removePacket(lines, keyword):
    _, leftoverLines = splitPacket(lines, keyword)
    return leftoverLines

def extractPackets(lines, keyword):
    extractedLines, _ = splitPacket(lines, keyword)
    return extractedLines

def extractKeywords(orderPacket):
    return [splitline(line,"=")[0] for line in extractLinesInLevel(orderPacket, 0)]

def readDict(file):
    msglvl = 0
    newline()
    message("Reading input dictionary", msglvl)
    msglvl += 1
    lines = refineLines(file)

    for keyword in extractLinesInLevel(lines,1):
        ps = extractPackets(lines, keyword)
        if len(ps) == 0: fataError("Cannot find any entries for the keyword " + keyword + " in the " + mainPacketKeyword() +
             " dictionary.\nExiting..")
        if keyword == dirsKeyword():
            dirs, labels, bases = readDirectories(ps)
            message("Read directories", msglvl)
        elif keyword == outkeyword():
            outdir = ps[0]
            message("Read output directory", msglvl)
        elif keyword == refCaseKeyword():
            refcase = ps[0]
            message("Read reference case", msglvl)
        elif keyword == turbLocKeyword():
            domainIndex, turbDataFileLoc = readTurbLoc(ps)
            message("Read turbines data files", msglvl)
        elif keyword == ordersKeyword():
            orders = readOrders(ps, msglvl)

    return globe(dirs, labels, bases, refcase, outdir, domainIndex, turbDataFileLoc), orders

def readTurbLoc(lines):
    domainIndex = []
    turbDataFileLoc = []
    for line in lines:
        entry = splitline(line,"=")
        identifier = entry[0] 
        if identifier[0] == "d": domainIndex.append(int(identifier[1]))
        elif identifier == turbCoordIdentifier(): domainIndex.append(turbCoordCode())
        elif identifier == turbineTypesTablesIdentifier(): domainIndex.append(turbTypeCode())
        elif identifier == windFarmsnamesIdentifier(): domainIndex.append(wfnamesCode())
        turbDataFileLoc.append(entry[1])

    return domainIndex, turbDataFileLoc

def readOrders(lines, msglvl):
    ordersList = []
    for order in extractLinesInLevel(lines, 0):

        identValues = []
        identNames = []
        intervalsList = []
        orderPacket = extractPackets(lines, order)
        #keywords = extractKeywords(orderPacket)
        #activeIndex = keywords.index(activeIdentifier())
        #if splitline(entry, "=")

        for keyword in extractKeywords(orderPacket):
            if keyword in subpacketIdentifiers():
                subpacket = extractPackets(orderPacket, keyword)
                if keyword == intervalsIdentifier():
                    intervalsNames = extractLinesInLevel(subpacket, 0)
                    for name in intervalsNames:
                        intervalsList.append(timeInterval(name, extractPackets(subpacket, name)))
                    identValues.append(intervalsList)
                elif keyword in [pairsIdentifier(), segmentsIdentifier(), assistFilesIdentifier()]:
                    identValues.append(subpacket)
                elif keyword == durIdentifier():
                    duration = timeInterval(keyword, subpacket)
                    identValues.append(duration)
                identNames.append(keyword)
                orderPacket = removePacket(orderPacket, keyword)

        for entry in orderPacket:
            splitEntry = splitline(entry, "=")
            identifier = splitEntry[0]
            value = splitEntry[1]
            if identifier == activeIdentifier():
                orderActive = value.lower() == "yes"
            elif identifier == typeIdentifier():
                orderType = value
            elif identifier == wfIndexIdentifier():
                identValues.append([int(q) for q in splitline(value, ",")])
            elif identifier in [locIdentifier(), refWindFarmName(), refFileName(), 
                                dronPathIdentifier(), heightIdentifier()]:
                identValues.append(value) 
            elif identifier in [boundsIdentifier(), zoomIdentifier(), capsIdentifier(), levelsIdentifier()]:
                identValues.append([float(q) for q in splitline(value, ",")]) 
            elif identifier == timeInstantIdentifier():
                identValues.append(timeInstant(value))

            if not identifier in [activeIdentifier(), typeIdentifier()]: identNames.append(identifier)
        
        if orderActive: ordersList.append(register(order, orderType, identValues, identNames, msglvl))
    return ordersList

def readDir(dict):
    # read opening bracket
    line = dict.readline()

    # readlines until the closing bracket
    while True:
        try: line = splitline(dict.readline(),"\n")
        except: break
        if line[0] == ")": break
        if line[0][0] == "/": continue

        outdir = line[0]
    return outdir
            
def readDirectories(lines):
    
    # placeholders for directories, labels, and base name
    dirs = []
    labels = []
    bases = []
    for line in lines:
        splitdata = splitline(line,",")
        if int(splitdata[3]) == 1:
            dirs.append(splitdata[0])
            labels.append(splitdata[1])
            bases.append(splitdata[2])
    return dirs, labels, bases