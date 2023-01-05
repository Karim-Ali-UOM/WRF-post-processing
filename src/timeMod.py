from src.utilities import splitline, startIdentifier, endIdentifier

def clockFormat(clock):
    hrs = int(clock)
    leftmin = (clock - hrs)*60
    mins = int(leftmin)
    leftsec = (leftmin - mins)*60
    sec = int(leftsec)
    return str(hrs) + ":" + str(mins) + ":" + str(sec)

def fileToHour(file, start):

    s1 = splitline(splitline(file,"/")[-1], "_")
    s2 = splitline(s1[-2], "-")
    return timeInstant(s2[2] + "-" + s2[1] + "-" + s2[0] + ", " + s1[-1]) - start

def isLeap(year):
    # rule out centuries
    if year % 100 == 0:
        # except for years divisble by 400 
        if year % 400 == 0:
            return 1
        else:
            return 0
    elif abs(year-2000) % 4 == 0:
        return 1
    else:
        return 0

def daysPerMonth(month, year):
    if month in [1, 3, 5, 7, 8, 10, 12]:
        return 31
    elif month in [4, 6, 9, 11]:
        return 30
    if month == 2:
        return 28 + isLeap(year)

def daysApart(ins1, ins2):
    if ins1 < ins2:
        small = ins1.copy()
        big = ins2.copy()
        s = -1
    else:
        small = ins2.copy()
        big = ins1.copy()
        s = 1
    small.hour = 0
    small.min = 0
    small.sec = 0
    big.hour = 0
    big.min = 0
    big.sec = 0
    ndays = 0

    while small < big:
        ndays += 1
        small.addDay()
    return ndays * s

def selectTimeFrames(files, interval):
    return [file for file in files if interval.within(wrfoutToTimeInstance(file))]

def selectNearestTimeFrames(files, instant):
    diffBelow = 10000
    diffAbove = 10000
    filebelow = -1
    fileabove = -1
    for file in files:
        diff = wrfoutToTimeInstance(file) - instant
        if diff > 0:
            if diff < diffAbove:
                diffAbove = diff
                fileabove = file
        elif diff < 0:
            if abs(diff) < diffBelow:
                diffBelow = abs(diff)
                filebelow = file
        else:
            return [file], [1]

    dtot = wrfoutToTimeInstance(fileabove) - wrfoutToTimeInstance(filebelow)
    dabove = wrfoutToTimeInstance(fileabove) - instant
    dbelow = instant - wrfoutToTimeInstance(filebelow)

    return [filebelow, fileabove], [dabove/dtot, dbelow/dtot]

def wrfoutToTimeInstance(file):
    splitname = splitline(file,"_")
    time1 = splitline(splitname[-2], "-")
    time2 = splitline(splitname[-1], ":")
    return timeInstant(int(time1[0]), int(time1[1]), int(time1[2]),
                    int(time2[0]), int(time2[1]), int(time2[2]))

class timeInstant:
    year = int
    month = int
    day = int
    hour = int
    min = int
    sec = int

    def __init__(self,*args):
        if len(args) > 1:
            self.year = int(args[0])
            self.month = int(args[1])
            self.day = int(args[2])
            self.hour = int(args[3])
            self.min = int(args[4])
            self.sec = int(args[5])
        else:
            splitdata = splitline(args[0],",")
            time1 = splitline(splitdata[0],"-")
            time2 = splitline(splitdata[1],":")
            self.year = int(time1[2])
            self.month = int(time1[1])
            self.day = int(time1[0])
            self.hour = int(time2[0])
            self.min = int(time2[1])
            self.sec = int(time2[2])

    def __sub__(self, other):
        return daysApart(self, other) * 24 + \
               self.hour - other.hour + \
               (self.min - other.min) / 60.0 + \
               (self.sec - other.sec) / 3600.0

    def __gt__(self, other):
        if self.year != other.year: return self.year > other.year
        elif self.month != other.month: return self.month > other.month 
        elif self.day != other.day: return self.day > other.day 
        elif self.hour != other.hour: return self.hour > other.hour 
        elif self.min != other.min: return self.min > other.min 
        elif self.sec != other.sec: return self.sec > other.sec 
        elif self.month != other.month: return self.month > other.month 
        else: return False

    def __lt__(self, other):
        if self.year != other.year: return self.year < other.year
        elif self.month != other.month: return self.month < other.month 
        elif self.day != other.day: return self.day < other.day 
        elif self.hour != other.hour: return self.hour < other.hour 
        elif self.min != other.min: return self.min < other.min 
        elif self.sec != other.sec: return self.sec < other.sec 
        elif self.month != other.month: return self.month < other.month 
        else: return False

    def __eq__(self, other):
        return (self.year == other.year and
                self.month == other.month and
                self.day == other.day and
                self.hour == other.hour and
                self.min == other.min and
                self.sec == other.sec)

    def copy(self):
        return timeInstant(self.year,self.month,self.day,self.hour,self.min,self.sec)

    def after(self, year=0, month=0, day=0, hour=0, min=0, sec=0):
        other = self.copy()
        other.year += year
        other.month += month
        other.day += day
        other.hour += hour
        other.min += min
        other.sec += sec 
        return other

    def addDay(self):
        self.day += 1
        dpm = daysPerMonth(self.month, self.year)
        if self.day > dpm:
            # increment month and adjust day
            self.month += 1
            self.day -= dpm

        if self.month > 12:
            # increment year and adjust month
            self.year += 1
            self.month -= 12

class timeInterval:
    start = timeInstant
    end = timeInstant
    name = ""

    def __init__(self,name,*args):
        self.name = name
        if len(args) == 1:
            for entry in args[0]:
                splitdata = splitline(entry,"=")
                if splitdata[0] == startIdentifier():
                    self.start = timeInstant(splitdata[1])
                elif splitdata[0] == endIdentifier():
                    self.end = timeInstant(splitdata[1])
        elif len(args) == 2:
            if args[0] < args[1]:
                self.start = args[0]
                self.end = args[1]
            else:
                self.start = args[1]
                self.end = args[0]

    def within(self,instant):
        return (self.start < instant < self.end) or (self.start == instant) or (self.end == instant )

    def copy(self):
        return timeInterval(self.start,self.end,self.name)