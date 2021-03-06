# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 19:33:37 2019

@author: myran
"""
import matplotlib.pyplot as plt
from seabird import fCNV
import numpy as np
import sys
import os
import cartopy.crs as ccrs
from shapely.geometry.polygon import LinearRing
import cartopy.feature as cfeature
sys.path.insert(0, '/Users/myran/Documents/Oceanography/Ocean 340')
from Myranda_Dang_a6_340b_cross_sections import twoDemData, calcDistance

def allProfDict(profList, varList): 
    '''
    This function creates a "master dictionary" and returns a dictionary containing all the 
    stations and their values only during the downcast. Also gets list of max depths of each
    stations
    '''
    ProfDict = {}
    for i in range(len(profList)):
        prof = getDownCastData(profList[i], varList)    #one station
        ProfDict['Station ' + str(i + 1)] = prof
    return(ProfDict)

def getDownCastData(dataset, varList):
    '''
    This function takes in the filepath name of the station's cnv file, creates a dictionary of 
    the station variables and values only in the CTD downcasting, and returns that dictionary
    '''
    station = fCNV(dataset)
    newProf = {}
    lat, lon = latLong(station)
    newProf['LATITUDE'] = lat
    newProf['LONGITUDE'] = lon
    for var in station.keys():
        newProf[var] = []
    depth = station['DEPTH']
    maxNum = depth.max()
    end = list(depth).index(maxNum)
    minNum = depth.min()
    start = list(depth).index(minNum)
    for var in varList:    #Look at each variable
        for i in range(start, end):
            newProf[var].append(station[var][i])
    return(newProf)
            
def latLong(stationDict):
    '''
    This function takes in a dicionary containing a station and its variables and values, and 
    returns the latitude and longitude (in that order) of where the profile was taken with rounded 
    decimals of four digits after the decimal place
    '''
    head = stationDict.attributes
    lat = float(head['LATITUDE'])
    lat = round(lat, 4)
    long = float(head['LONGITUDE'])
    long = round(long, 4)
    return lat, long

def plotCrossSection(depths2d, dist2d, var2d, figNum, xlabel, ylabel, title, barLabel, maxDepth):
    '''
    This function takes in a 2-D array of depth, distance,
    variable data, title, and a figure number, and create a filled contour plot of the 
    variable data array with respect to the depths and distance arrays given
    '''
    plt.figure(num = figNum)
    plt.gca().invert_yaxis()
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title, fontsize = 15)
    plt.contourf(dist2d, depths2d, var2d)
    #plt.fill_between(dist2shore,bathy,1000,color='k') #Mapping bathymetry
    cbar = plt.colorbar()
    cbar.set_label(barLabel)
    plt.show()
    
def plotVar(staDict, varName, figNum, yVar = 'DEPTH'):
    """This function takes in a profile of a variable in the form of a list, another variable
    in a list (most likely depth), the figure number, and the label for the x-axis. It then plots 
    the profile on a specific figure, adding labels, a title, inverts the y-axis for proper 
    depth profiles, and caption.
    """
    plt.figure(figNum, figsize = (15,20))
    for station in staDict.keys():
        var1 = staDict[station][varName]
        var2 = staDict[station][yVar]
        plt.plot(var1, var2, label = station)
    plt.xlabel(varName)
    plt.ylabel(yVar)
    if 'DEPTH' in yVar:
        plt.gca().invert_yaxis()
        plt.title(varName + ' Profile Plot')
    else:
        plt.title(varName + ' vs '+yVar)
    plt.legend()
    plt.show()
    #Optional saving profiles
    #file = xlabel + '.png'
    #plt.savefig(file)
    
def mapStations(lats, longs, figSize=(20,20), res='10m', cmap='rainbow'):
    '''
    This function takes in lists of longitude and latitude and plots 
    them on a map, as well as taking in optional parameters such as
    figure size, default of (20, 20), map resolution (as string), 
    default of '10m', and colormap (string), default of 'rainbow'
    '''
    fig = plt.figure(1,figsize=figSize)
    fig.suptitle('Station Locations', fontsize=25)
    map1 = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree())
    map2 = fig.add_subplot(1,2,1,projection=ccrs.PlateCarree())
    map1.coastlines(resolution=res)
    map1.set_extent([-124, -122, 47, 48.5]) #Puget Sound area 
    pslons = [-122.10, -122.10, -122.40, -122.40]
    pslats = [47.89, 48.07, 48.07, 47.89]
    ring = LinearRing(list(zip(pslons, pslats)))
    map2.coastlines(resolution=res)
    map2.set_extent([-122.40, -122.10, 47.89, 48.07], crs=ccrs.PlateCarree())    #Posession Sound
    land = cfeature.NaturalEarthFeature('physical', 'land', res, facecolor='darkseagreen')
    map1.add_feature(land)
    map2.add_feature(land)
    map1.add_geometries([ring],ccrs.PlateCarree(),facecolor='none',edgecolor='darkred',linewidth=4)
    color = eval('plt.cm.'+cmap+'(np.linspace(0,1, len(lats)))')
    for i in range(len(lats)):
        c = color[i]
        lat = lats[i]
        long = longs[i]
        profNum = i+1 
        l = "Station" + str(profNum)
        map2.plot(long,lat,'o',color=c, label = l, markersize=14, transform=ccrs.PlateCarree())
    map2.legend()
    map1.set_title('Puget Sound', fontsize=14)
    map2.set_title('Possession Sound', fontsize=14)
    plt.show

def mapDepthCrossSection(lats,longs,varName,stationName,varData,barLabel):
    '''
    This function takes in a list of latitudes and longitudes, the variable name and data that 
    are being plotted and shown as changes in each depth layer for each station, label containing
    the units of the inputted variable, and the figure number.
    
    '''
    plt.figure(figsize=(20,20))
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_title(varName + ' m', fontsize=25)
    land = cfeature.NaturalEarthFeature('physical', 'land', '10m', 
                                        facecolor='darkseagreen')
    ax.set_extent([-122.40, -122.10, 47.89, 48.07], crs=ccrs.PlateCarree())
    ax.coastlines(resolution='10m')
    ax.add_feature(land) 
    plt.scatter(longs,lats,c=varData,s=200,transform=ccrs.PlateCarree())
    #Need to check to see if SAL/AATT = 0, then pt would be black 'k'
    for i in range(len(lats)):
        lat = lats[i]
        long = longs[i]
        if varData[i] == 0: #Check if station has data or not
            plt.plot(long,lat,'o',color='darkgray', markersize=14, transform=ccrs.PlateCarree())
        plt.text(longs[i], lats[i]+0.005, 'Sta. '+ str(stationName[i]),
                fontsize=14,transform=ccrs.PlateCarree())
    cbar = plt.colorbar()
    cbar.set_label(barLabel)
    plt.show()
    
def depthRange(depth, station, varName):
    '''
    This function takes in depth(as a string) that samples were collected at, a list of the station, 
    their variables and values, and also the variable name of which the user
    wishs to calculate; calculates an average variable over the passed in depth and returns average
    as a float
    '''
    avgList = []
    count = 0
    for i in range(len(station['DEPTH'])):
        num = station['DEPTH'][i]
        if depth in str(num):
            if num < (float(depth) + 2) and num >= (float(depth) - 2):
                value = list(station['DEPTH']).index(num)
                avgList.append(station[varName][value])
                count += 1
    result = 0
    for j in range(len(avgList)):
        result += avgList[j]
    if result != 0 and count != 0:
        result = result/count
    else:
        result = 0
    return result
  
def createLayersDict(varName, stationDict, layerDepths):
    '''
    This function takes in a variable name, a dictionary containing all the stations and their
    values, and a list of depths that are considered shallow, mid, and deep layers. The function
    creates a dictionary of variable values at those depth layers and returns created dictionary.
    '''
    layerVar = {}
    for layer in layerDepths:
        layerVar[layer] = []
        for station in stationDict.keys():
            if station not in 'Everett Docks':
                layerData = depthRange(str(layer),stationDict[station],varName)
                layerVar[layer].append(layerData)
    return layerVar

def main(): 
    #Creating list of profile filenames for Stations 1-9
    fp = '/Users/myran/Documents/Oceanography/Ocean Senior Thesis/DATA/rc0015_CTD_MD/PROCESSED/Stations 1-9/'
    fileList = os.listdir(fp)
    profList = []
    for item in fileList:
        if '.cnv' in item:    #Looks to see if file is .cnv file
            openFile = fp + item
            profList.append(openFile)
    
    #List of variables 
    varList = ['PSAL', 'CStarAt0', 'CStarTr0', 'DEPTH', 'flECO-AFL']

    #Attenuation and Salinity labels
    attLabel = 'Attenuation (1/m)'
    salLabel = 'Salinity (PSU)'
    
    #Create dictionary of all profile
    d = allProfDict(profList, varList)
    
    #Create list of max depths of all stations
    maxDepth = []
    for station in d.keys():
        maxNum = max(d[station]['DEPTH'])
        maxDepth.append(maxNum)

    #Append Everett Dock file to master dictionary
    dockName = '/Users/myran/Documents/Oceanography/Ocean Senior Thesis/DATA/rc0015_CTD_MD/PROCESSED/RC0015_MD_dockTest.cnv'
    dockFile = getDownCastData(dockName,varList)
    d['Everett Docks'] = dockFile
    
    # Create numpy array of new depths
    newDepth = np.arange(10, 60, 2)
    
    #Getting stations that are in a line (Stations 2 to 9) for crossections
    statLine = {}
    stationList = list(d.keys())
    for station in stationList[1:9]:
        statLine[station] = d[station]
    #Create array of latitude and longitude
    lats = []
    longs = []
    for station in d.keys():
        lat = d[station]['LATITUDE']
        lon = d[station]['LONGITUDE']
        lats.append(lat)
        longs.append(lon)
        
    #Get lats/longs of stations to for transect
    crossLats =[]
    crossLongs = []
    crossDepth = []
    for i in range(4):
        latCross1 = lats[1 + i]
        depth1 = lats[1 + i]
        crossLats.append(latCross1)
        crossDepth.append(depth1)
        longCross1 = longs[1 + i]
        depth2 = longs[1 + i]
        crossLongs.append(longCross1)
        crossDepth.append(depth2)
        latCross2 = lats[8 - i]
        crossLats.append(latCross2)
        longCross2 = longs[8 - i]
        crossLongs.append(longCross2)
     
    # Get list of distances between datapoints
    dists = calcDistance(crossLats, crossLongs)
    
    # Create 2-D array of data salinity and attenuation
    attData = twoDemData('CStarAt0', statLine, newDepth)
    salData = twoDemData('PSAL', statLine, newDepth)
    
    # Create grid of depths and distances
    depths_2d, dists_2d = np.meshgrid(newDepth, dists)
    
    #Plot lat/long pair of all stations
    mapStations(lats, longs)

    #Plotting cross section
    plotCrossSection(depths_2d,dists_2d,attData,2,'Distance (km)','Depth(m)',
                     'Attenuation CrossSection',attLabel,crossDepth)
    plotCrossSection(depths_2d,dists_2d,salData,3,'Distance (km)','Depth(m)',
                     'Salinity CrossSection', salLabel,crossDepth)
    

    #Creating "shallow, mid, deep" data of Attenuation and Salinity 
    layerDepth = [10,20,30,40,50,60,70,80,90]
    attLayers = createLayersDict('CStarAt0', d, layerDepth)
    salLayers = createLayersDict('PSAL', d, layerDepth)

    #Getting lats/longs for depth layers
    depthCoorLats = lats
    depthCoorLongs = longs
    depthCoorLats.pop()
    depthCoorLongs.pop() 
    
    #Plotting depth interval layers 
    count = 4
    for layer in attLayers.keys():
        mapDepthCrossSection(depthCoorLats,depthCoorLongs,'Attenuation: '+str(layer),
                             [1,2,3,4,5,6,7,8,9],np.array(attLayers[layer]),attLabel)
        count += 1
        mapDepthCrossSection(depthCoorLats,depthCoorLongs,'Salinity: '+str(layer),
                             [1,2,3,4,5,6,7,8,9],np.array(salLayers[layer]),salLabel)
        count += 1
    '''
    #Creating depth profiles
    for var in varList: #plotting profiles for each variable
        plotVar(d, var, count)
        count += 1
    plotVar(d, 'PSAL', count, yVar = 'CStarAt0')
    '''   
    
if __name__ == "__main__":
    main()
    
    
    
