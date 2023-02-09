'''
Code to model clam growth based on temperature and d18o

Author: (Original) Sawyer Hilt, (Modifications) Ashwin Lall, (Technical Assist) Dave Goodwin

'''

import pip
pip.main(["install","matplotlib"])
pip.main(["install","pandas"])
import matplotlib.pyplot as pyplot
import math
from statistics import NormalDist
import pandas as pd
import random
import datetime

def periodicFunction(a,b,c,d,t):
    """
    Purpose: this function simulates temperature variations using periodic function shell
    inputs:
    a= amplitude
    b= perioid
    c= phase shift
    d= mean value
    t= number of iterations of periodic function / time period
    output:
    Y =  list of periodic variations
    """
    
    Y = []
    for i in range(len(t)):
        y = (a*(math.sin(b*(i+c)))+d)
        Y.append(y)
    return Y

def genHrFromDay(d, temp_amplitude_day = 2.5):
    """
    Purpose: this function takes daily simulated temperatures and creates 24 hours of temps from them, using
    periodic function and daily simulated temp as mean temp for each 24 hr increment
    
    input:
    d= list of mean daily temperatures
    temp_amplitude_day = temperature amplitude over a day
    
    output:
    hourData = list of mean hourly temperatures
    """
    hourData = []
    for i in range(len(d)):
        out = periodicFunction(temp_amplitude_day, math.pi/12, 3, d[i],list(range(24)))
        hourData = hourData + out

    return hourData

def genHourlyTemps(noYears, temp_mean = 20, temp_amplitude_year = 7.5, temp_amplitude_day = 2.5):
    """
    Purpose: use genHrFromDay function to create hourly temperatures for specified time period

    input:
    noYears = number of years to run 'experiment'
    temp_mean = mean temperature for the year
    temp_amplitude_year = temperature amplitude over the year
    temp_amplitude_day = temperature amplitude over a day

    output:
    d18oShellVPDB = hourly d18o values for clam shells, properly converted to d18o values on VPDB scale
    """

    meanDailyTemps = periodicFunction(temp_amplitude_year, (2*math.pi)/365, -90.75, temp_mean, list(range(0,365))) * noYears
    hourlyTemps = genHrFromDay(meanDailyTemps, temp_amplitude_day)
    return hourlyTemps

def genHourlyWaterd18o(noYears, mean, amplitude):
    """
    Purpose: Generates hourly d18o for water using a sine curve with peak in the summer of maxd18o

    input:
    noYears = number of years to run 'experiment'
    mean, amplitude: mean and amplitude of sine curves

    output:
    d18oWater = hourly d18o values for water
    """
    return periodicFunction(amplitude, (2*math.pi/(365*24)), -90.75*24, mean, list(range(365 * 24))) * noYears

def genHourlyd18o(temps, d18oWater):
    """
    Purpose: this function takes a given list of temperatures and converts them to d18o values of modelled clam
    shells using Dettman et al's 1999 paleotemperature equation.
        ** function assumes that the d18o value of the water in VSMOW is = 0 for each sample

    input:
    Temps = list of temperatures to be converted

    output:
    d18oshellVPDB = temperature lists converted to list of d18o values
    """
    #d18oWater = 0 # assuming a constant d18o value of water of 0 (on VSMOW scale)
    d18oShellVPDB = []
    for i in range(len(temps)): # for loop that converts temperature for each item in 'Temps' list
        tempK = temps[i] + 273.15 
        dettman = (2.559*((10**6)*(tempK)**(-2))+0.715) #using dettman's 1999 equation to find alpha lna100 given temp of water
        lna = dettman / 1000 # dividing by 1000 to get ln(a)
        alpha = math.exp(lna) # converting ln(a) to a (alpha)
        d18oVSMOW = alpha*(1000 + d18oWater[i]) - 1000 # converting to d18o in VSMOW
        d18oVPDB = ((d18oVSMOW - 30.91)/ (1.03091)) # converting d18o from VSMOW to VPDB scale
        d18oShellVPDB.append(d18oVPDB) # appending to list 

    return d18oShellVPDB

def hourlyIncWidth(temps, ogt, sd, max_growth, shutdowns, growthRate, startJulian, endJulian, springIncrease):
    """
    Purpose: calculates the width of each hourly increment of shell, based on the optimal temperature.
    Ex: An hour that has a temperature equal to the optimal growth temperature will produce the greatest
    width of incremental growth. Increment widths for temperatures which are not equal to the optimal
    growth temperature will be smaller. 
    
    Inputs:
    temps = list of temps for all hours of simulated sampling
    ogt = optimal growth temperature
    sd = standard deviation
    max_growth = maximum growth... growth which occurs if temp = ogt
    shutdowns = lengths of time that growth shuts down each year; first must always be 0; after the last
    year in the list, all subsequent shutdowns are the same as the last
    growthRate = fraction of the maximum possible growth (compared to the first year) for each year; first
    must always be 0; after the last year in the list, all subsequent rates are the same as the last
    Parameters for nutrience abundance in the spring:
    startJulian = start of spring
    endJulian = end of spring
    springIncrease = factor by which growth is higher in the spring

    Output:
    the width of each hourly increment of shell
    """
    hourlyWidth = []
    OGT_width = NormalDist(mu = ogt, sigma = sd).pdf(ogt)
    for i in range(len(temps)):
        calc = NormalDist(mu = ogt, sigma = sd).pdf(temps[i])
        weighted = calc / OGT_width
        width = weighted * max_growth

        # always leave the first half year alone and start counting years from there
        if i < 365 * 24 // 2:
            year = 0
        else:
            year = (i - 365 * 24 // 2) // (365 * 24) + 1

        if year < len(shutdowns):
            shutdown = shutdowns[year]
        else:
            shutdown = shutdowns[-1]  # use the last shutdown for each subsequent year

        # shutdown growth at the beginning and end of the year
        if i % (365 * 24) < shutdown/2 * 30 * 24 or i % (365 * 24) >  365 * 24 - shutdown/2 * 30 * 24:
            width = 0

        # account for slowdown in growth rate
        calendarYear = i // (365 * 24)
        if calendarYear >= len(growthRate):
            width = width * growthRate[-1]
        else:
            width = width * growthRate[calendarYear]

        # account for abundance of nutrients in the spring
        julianDay = i % (365 * 24) // 24
        if julianDay >= startJulian and julianDay <= endJulian:
            width = width * springIncrease
        
        hourlyWidth.append(width)

    return hourlyWidth
  

def genPotentialSample(hourlyIncWidth, hourlyd18o, sampleSize):
    """ Edited function description 08.22.22
    Purpose: function creates samples from given hourly increment widths, d18o values,
    and sample sizes
    - sample sizes are equal to the smallest number of increment widths whose sum was greater
    than or equal to the required increment width.
    

    inputs:
    hourlyIncWidth = list of hourly increment widths
    hourlyd18o = list of hourly d18o values from shell
    sampleSize = how large should sample be? ranges from 50-150 micrometers

    outputs:
    sampleEndHour = hour in which sampling ends. 
        Ex: if increment widths for hours 1-4 are 100, 105, 110, and 95 micrometers, and the sample size is 400 micrometers,
        this value will be equal to 4.
        (Ideally, the hour catalogued for the sample is weighted using the different increment widths which contributed to
        to each sample)
    d18oSamples = d18o values of samples of specified size taken from shells 
    """
    samples = []
    d18oSamples = []
    sampleEndHour = []
    EndHours = []
    hours = 0
    sampleWidth = []
    j = 0

    while j < len(hourlyIncWidth):
        
        sampleWidth = 0
        k = j 
        
        while sampleWidth < sampleSize and j < len(hourlyIncWidth):
            sampleWidth = sampleWidth + hourlyIncWidth[j]
            j = j + 1

        hours = j - k
        
        if sampleWidth < sampleSize:
            return d18oSamples, sampleEndHour
        
        d18oSample = (sum([hourlyIncWidth[i] * hourlyd18o[i] for i in range(k,j)])) / sampleWidth # weighted sample

        if sampleWidth >= sampleSize:
            d18oSamples.append(d18oSample)
            sampleEndHour.append(j)
    
        

    return d18oSamples, sampleEndHour
            
            
                                                                              
def sampleShell(noShells, temps, mu_ogt, sd_ogt, mu_ogt_sd, sd_ogt_sd, maxHrWidth, shutdowns, growthRate, startJulian, endJulian, springIncrease, d18oWater):
    """
    Purpose: sampleShell uses avg hourly increment widths calculated with hourlyIncWidth function to 

    
    Inputs:
    allShells = empty list
    noSample = number of shells to be sampled
    temps = list of temperatures
    mu_ogt = mu for the OGT mean 
    sd_ogt = sd for the OGT mean
    mu_ogt_sd = mu for the OGT sd
    sd_ogt_sd = sd for the OGT sd
    maxHrWidth = maximum hourly growth
    shutdowns = lengths of time that growth shuts down each year; first must always be 0; after the last
    year in the list, all subsequent shutdowns are the same as the last
    growthRate = fraction of the maximum possible growth (compared to the first year) for each year; first
    must always be 0; after the last year in the list, all subsequent rates are the same as the last
    Parameters for nutrience abundance in the spring:
    startJulian = start of spring
    endJulian = end of spring
    springIncrease = factor by which growth is higher in the spring
    d18oWater = d18o of Water

    Outputs:
    allShells = list of lists of growth per hour for each shell
    allTimes = list of lists of timestamps for each hour for each shell
    OGT_mu = the actual OGTs for each shell
    OGT_sd = the sd for the OGT for each shell
    """
    allShells = []
    allTimes = []
    OGT_mu = NormalDist(mu_ogt, sd_ogt).samples(noShells)         # generate the OGT for all specimen
    OGT_sd = NormalDist(mu_ogt_sd, sd_ogt_sd).samples(noShells)   # generate the OGT sd for all specimen
    for i in range(noShells):
        hourlyWidth = hourlyIncWidth(temps, OGT_mu[i], OGT_sd[i], maxHrWidth, shutdowns, growthRate, startJulian, endJulian, springIncrease) 
        hourlyd18o = genHourlyd18o(temps, d18oWater)
        samples, sampleEndHour = genPotentialSample(hourlyWidth, hourlyd18o, 300)
        allShells.append(samples)
        allTimes.append(sampleEndHour)
      
    return allShells, allTimes, OGT_mu, OGT_sd

def subSample(samples, times, windowSize, windowSD):
    """
    Purpose: Take a number of complete samples and subsamples one per window,
    picking samples near the center of the window with higher chance using a
    normal distribution
    
    samples = a list of completely sampled specimen
    times = timestamps for each sample
    windowSize = the size of a window
    windomSD = the standard deviation for sampling near the center of the window. A low value
    means that it will almost always pick near the center. A large value is closer to uniform.
    """

    subsamples = []
    subtimes = []
    for i in range(len(samples)):
        sample = samples[i]
        time = times[i]
        subsample = []
        subtime = []
        for i in range(0, len(sample), windowSize):
            # select the sample within this window assuming normally distributed weights across the window
            # with a mean at the center of the window and sigma equal to windowSD
            sampledIndex = random.choices(list(range(i, i + windowSize)), [NormalDist(mu = i + windowSize//2, sigma = windowSD).pdf(i) for i in range(windowSize)])[0]
            if sampledIndex < len(sample):  # may fall off on the last window
                subsample.append(sample[sampledIndex])
                subtime.append(time[sampledIndex])
        subsamples.append(subsample)
        subtimes.append(subtime)
    return subsamples, subtimes


def plotShells(samples, subsamples, times, noYears, sampleEndHour, hourlyTempsC, shutdowns, growthRate, startJulian, endJulian, springIncrease, d18oWater):
    """
    Purpose: plots all shell values onto one plot

    input:
    samples =  list of d18o samples of every shell
    subsamples = list of subsampled d18o values for each of the samples
    times = timestamps for each subsample
    noYears = number of yrs to run model for
    SampleEndHour = list of last hour of temp recording included in sample, for each sample
    hourlyTempsC = list of hourly temps for as long as experiment runs for
    shutdowns = lengths of time that growth shuts down each year; first must always be 0; after the last
    year in the list, all subsequent shutdowns are the same as the last
    growthRate = fraction of the maximum possible growth (compared to the first year) for each year; first
    must always be 0; after the last year in the list, all subsequent rates are the same as the last
    Parameters for nutrience abundance in the spring:
    startJulian = start of spring
    endJulian = end of spring
    springIncrease = factor by which growth is higher in the spring
    d18oWater = d18o of water
    
    output:
    shows two plots, one which displays ALL sampled shell values and last hr of sample taken on x-axis, and
    another plot which shows temperature over time 
    """
    pyplot.subplot(6, 1, 1)
    pyplot.plot(range(24*365*noYears), hourlyTempsC)
    pyplot.xlabel('Hour')
    pyplot.ylabel('Hourly Temperature (deg C)')

    pyplot.subplot(6, 1, 2)
    pyplot.plot(range(24*365*noYears), genHourlyd18o(hourlyTempsC, d18oWater))
    pyplot.xlabel('Hour')
    pyplot.ylabel('d18o')

    pyplot.subplot(6, 1, 3)
    pyplot.plot(range(24*365*noYears), d18oWater)
    pyplot.xlabel('Hour')
    pyplot.ylabel('d18o Water')

    pyplot.subplot(6, 1, 4)
    pyplot.plot(range(24*365*noYears), hourlyIncWidth(hourlyTempsC, 25, 6, 5, shutdowns, growthRate, startJulian, endJulian, springIncrease))
    pyplot.xlabel('Hour')
    pyplot.ylabel('Hourly growth')

    pyplot.subplot(6, 1, 5)
    pyplot.plot(samples[0])
    pyplot.xlabel('Time stamp')
    pyplot.ylabel('Real d18o value (VPDB)')

    pyplot.subplot(6, 1, 6)
    pyplot.plot(subsamples[0])
    pyplot.xlabel('Subsample number')
    pyplot.ylabel('Sample d18o value (VPDB)')
    
    pyplot.show()



def getStatsByYear(sample, times, BLANK = ""):
    """
    Purpose: compute statistics by year for a given sample

    input:
    samples =  list of d18o samples of a shell
    times = timestamps for each sample
    BLANK = entry for a missing field
    
    output:
    a list of statistics for each year; each element of the list is a tuple with:
    min, max, range
    """
    dataByYear = []
    for i in range(len(sample)):
        year = times[i] // (24 * 365)
        while year >= len(dataByYear):
            dataByYear.append([])
        dataByYear[year].append(sample[i])

    #print(sample)
    #print(times)

    stats = []
    for i in range(len(dataByYear)):
        data = dataByYear[i]
        #print("Year", i, ":", data)
        if len(data) > 0:
            stats.append((min(data), max(data), max(data) - min(data)))
        else:
            stats.append((BLANK, BLANK, BLANK))
    return stats
                     

def saveData(samples, times, OGT_mu, OGT_sd, timestamp, SEP = ","):
    """
    Purpose: Dump data into files

    input:
    samples =  list of d18o samples of a number of shells
    times = timestamps for each sample
    OGT_mu = OGT for the shells
    OGT_sd = sd of the OGT for the shells
    SEP = separator for the CSV files
    timestamp = timestamp of this run
    
    output:
    a file with the OGT stats followed by the sample readings for each shell
    a file with some stats per year for each shell
    """

    # dump data for OGT stats and raw sample readings
    raw_data = open(timestamp + "-raw_data.csv", "w")   
    raw_data.write("OGT" + SEP + "sd\n")
    for i in range(len(samples)):
        raw_data.write(str(OGT_mu[i]) + SEP + str(OGT_sd[i]))
        for j in range(len(samples[i])):
            raw_data.write(SEP + str(samples[i][j]))
        raw_data.write("\n")
    raw_data.close()

    # dump statistics per year for each specimen
    yearly_min = open(timestamp + "-yearly_min.csv", "w")   # this clobbers any existing file in the folder
    yearly_max = open(timestamp + "-yearly_max.csv", "w")   # this clobbers any existing file in the folder
    yearly_range = open(timestamp + "-yearly_range.csv", "w")   # this clobbers any existing file in the folder
    all_stats = []
    for i in range(len(samples)):
        stats = getStatsByYear(samples[i], times[i])
        all_stats.append(stats)
    num_years = max([len(i) for i in all_stats])
    yearly_min.write("Iteration" + SEP + SEP.join(["Year" + str(i) for i in range(1, num_years + 1)]) + "\n")
    yearly_max.write("Iteration" + SEP + SEP.join(["Year" + str(i) for i in range(1, num_years + 1)]) + "\n")
    yearly_range.write("Iteration" + SEP + SEP.join(["Year" + str(i) for i in range(1, num_years + 1)]) + "\n")
    for i in range(len(samples)):
        stats = all_stats[i]
        #print(stats)
        mins = SEP.join([str(i[0]) for i in stats])
        maxs = SEP.join([str(i[1]) for i in stats])
        ranges = SEP.join([str(i[2]) for i in stats])
        yearly_min.write(str(i + 1) + SEP + mins + "\n")
        yearly_max.write(str(i + 1) + SEP + maxs + "\n")
        yearly_range.write(str(i + 1) + SEP + ranges + "\n")
    yearly_min.close()
    yearly_max.close()
    yearly_range.close()


def dumpShell(hourlyGrowth, timestamp):
    '''
    This dumps the shell's values to a CSV
    Ignores all entries that are 0

    input:
        hourlyGrowth: hourly growth of the shell
        timestamp = timestamp of this run
    '''
    file = open(timestamp + "-shell.csv", "w")
    for val in hourlyGrowth:
        file.write(str(val) + "\n")
    file.close()

def main():
    noYears = 10 # no. of years to run model
    noShells = 5 # no. of shells to be sampled

    temp_mean = 20               # mean temperature for the year
    temp_amplitude_year = 7.5    # temperature amplitude over the year
    temp_amplitude_day = 2.5     # temperature amplitude over a day

    # timestamp used for all the files created by this run
    timestamp = datetime.datetime.now().replace(microsecond=0).isoformat().replace("T", "-").replace(":", "-")

    #ogt = 25 # optimal growth temperature
    #sd = 5   # standard deviation
    mu_ogt = 25      # mu for the OGT mean 
    sd_ogt = 2       # sd for the OGT mean
    mu_ogt_sd = 5    # mu for the OGT sd
    sd_ogt_sd = 1    # sd for the OGT sd
    
    maxHrWidth = 5 # maximum growth per unit time

    # parameters for sampling from a window
    windowWidth = 11 # size of the window. all samples are divided into windows with this width.
    windowSD = 5     # standard deviation for sampling from this window. Small=> concentrated at center. High=>closer to uniform

    # shutdown is the number of months that the organism shuts down for during each winter period
    # the first value should always be 0 as you will not want to shut down in the first year
    shutdowns = [0]  # no shutdown at all
    shutdowns = [0, 0, 3, 4, 5, 6]
    #shutdowns = [0, 6, 8, 10, 12] # extreme shutdown for testing

    # growthRate is the fraction of the maximum possible growth (compared to the first year) for each year
    # the first value should always be 1
    growthRate = [1, 0.8, 0.6, 0.4, 0.2, 0.1]

    # parameters for spring nutrient abundance
    startJulian = 90
    endJulian = 150
    springIncrease = 2.0

    # water d18o parameter
    meanWaterd18o = 0.25
    amplitudeWaterd18o = 0.25
    d18oWater = genHourlyWaterd18o(noYears, meanWaterd18o, amplitudeWaterd18o)

    # generate hourly temps 
    hourlyTempsC = genHourlyTemps(noYears, temp_mean, temp_amplitude_year, temp_amplitude_day) # creating list of temperatures for time period
    
    # list which contains hourly temps and hr of last increment width for each sample
    samples, times, OGT_mu, OGT_sd = sampleShell(noShells, hourlyTempsC, mu_ogt, sd_ogt, mu_ogt_sd, sd_ogt_sd, maxHrWidth, shutdowns, growthRate, startJulian, endJulian, springIncrease, d18oWater) 

    # compute a subsample of each sample
    subsamples, times = subSample(samples, times, windowWidth, windowSD)

    #print(len(samples), len(samples[0]), len(subsamples[0]))

    # dump one specific sample into a file
    dumpShell(hourlyIncWidth(hourlyTempsC, 25, 6, 5, shutdowns, growthRate, startJulian, endJulian, springIncrease), timestamp)

    # print the stats for each sample, up to the first 5
    #for i in range(min(noShells, 5)):
    #    print(getStatsByYear(subsamples[i], times[i]))

    # dump the results into a file
    saveData(subsamples, times, OGT_mu, OGT_sd, timestamp)
    
    # plot all samples and compare against temp over time plot for same study period
    plotShells(samples, subsamples, times, noYears, list(range(len(samples))), hourlyTempsC, shutdowns, growthRate, startJulian, endJulian, springIncrease, d18oWater) 
    
    

    
main()
