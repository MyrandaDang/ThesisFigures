# Myranda Dang
# Problem Set #6: Cross Sections & Interpolation
# OCEAN 340B, Winter 2019

"""In this OCEAN 340B assignment we will examine how to create cross section plots. Part of this
process will involve interpolation, which is an estimation of a value between two given values.
Interpolation is very common in oceanographic computer programming, and cross sections are a very 
common way to represent oceanographic data. You will also be using the function given to you to
calculate the distance between two points (specifically within question 3). Please read the 
docstring and examine the code of this function to learn more about it.

The data you will be working with after you complete the functions (and should keep in mind while
writing the functions) is CNV data found on Canvas. These are profiles from around the Puget Sound.
You will read in the data using the fCNV function in your main function, which is created in 
question 5.

____________________________________________________________________________________________________                                                               
_________________________________________ QUESTIONS ________________________________________________
____________________________________________________________________________________________________ 

Question 1: This question asks you to write a function that interpolates data. Please use the
following example, courtesy of SciPy documentation:
    
https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html


Create a function that takes in an profile of variable data (ex. temperature),
a profile of depths, and a profile of new depths. All of these profiles should be
numpy arrays. This function should then interpolate the variable data relative to the new depths 
given, and return this new interpolation of the variable data in the form of a numpy array.

*** 10 points ***


Question 2: This function will create a 2-D array of vertical profiles.

Create a function that takes in a variable name, a list of fCNV datasets (NOT filepaths), and a  
numpy array of depths for data to be interpolated to. This function returns a 2-D numpy array of 
data from the given datasets. The data returned should be that of the given variable name, 
interpolated to the depths given, and returned in the form of a numpy array.

*** 20 points ***


Question 3: Create a function that takes in a list of latitudes and a list of longitudes. This 
function should then calculate the distance between two consecutive latitudes and longitudes
in the given lists. The output is a numpy array with the first value being zero, and the following
values being the distances between the given latitudes and longitudes.

Example:
    latList = [40, 40.5, 41]
    lonList = [120, 120.2, 120.4]
    calculate distance between point 1 (40, 120) and point 2 (40.5, 120.2)...
    calculate distance between point 2 (40, 120) and point 3 (40.5, 120.2)...
    OUTPUT:[  0., 58.13064838, 116.22461932]
    
*** 20 points ***


Question 4: Create a function that takes in a 2-D depths array, a 2-D distance array, a 2-D
variable data array, and a figure number. This function should create a filled contour plot of the 
variable data array with respect to the depths and distance arrays given. 
Be sure to use the figure number when plotting.

*** 15 points ***


Question 5: Create a main function that plots the given data for this assignment by doing the
            following:
    1. Open all of the data given for this assignment using the fCNV function from the seabird
        package and put the datasets in a list.
    2. Create a numpy array of new depths to interpolate to. This should be, for this problem set,
        the depths 5 - 100, inclusive. HINT: Try using the np.arange function for this.
    3. Create 2-D array of variable data, in this case both salinity and temperature, using 
        createTwoDimData function.
    4. Create a list of all latitudes from the CNV datasets and a list of all longitudes from the
        CNV datasets.
    5. Get the distance between each latitude and longitude pair. HINT: Use function from question 3.
    6. Create a 2-D array of distances and a 2-D array of new depths. HINT: Use np.meshgrid
    7. Plot the cross sections using contourf of both salinity and temperature from each of the
        three given CNV files.

*** 20 points ***

Please complete the following problems. There are 100 points total, with 15 points being used to 
grade code style/logic.
"""


import matplotlib.pyplot as plt
import numpy as np
from seabird import fCNV
from scipy import interpolate
from math import asin, sin, cos, sqrt, radians


def getDistBetween(lat1, lon1, lat2, lon2):
    """Given a latitude and longitude of a point, and a latitude and longitude
    of another point (all in the form of floats), this function returns the distance
    in kilometers between the two points.
    """
    # Radius of earth in kilometers
    radius = 6371.0
    
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 

    km = radius * c
    return km

def interpretData(depthProf, varProf, newDepth):
    '''
    This function takes in an profile of variable data, a profile of depths, and profile 
    of new depths, interpolate the variable data relative to the new depths 
    given, and return this new interpolation of the variable data in form of an array
    '''
    f = interpolate.interp1d(depthProf, varProf, fill_value="extrapolate")
    newVarProf = f(newDepth)
    return newVarProf

def twoDemData(varName, setDict, depth): #Depth what we want to interpolate over
    '''
    This function takes in a variable name, a list of fCNV datasets, and depth for data to be 
    interpolated and returns data of the variable name interpolated to the given depths in form 
    of a numpy array
    '''
    output = []
    for dataset in setDict.keys():
        #get data
        varProf = setDict[dataset][varName]
        depthProf = setDict[dataset]['DEPTH']
        #interpret new depth
        newVarProf = interpretData(depthProf, varProf, depth)
        output.append(newVarProf) #make cross section
    return output

def calcDistance(lats, lons):
    '''
    This function takes in a list of latitudes and a list of longitudes, calculate the 
    distance between two consecutive latitudes and longitudes, and returns an array of distance
    with the first value as zero
    '''
    output = [0]
    for i in range(len(lons) - 1):
        lat1 = lats[i]
        lon1 = lons[i]
        lat2 = lats[i + 1]
        lon2 = lons[i + 1]
        dist = getDistBetween(lat1, lon1, lat2, lon2)
        output.append(dist + output[i])
    return np.array(output)

def plotCrossSection(depths2d, dist2d, var2d, figNum):
    '''
    This function takes in a 2-D array of depth, distance, and
    variable data, and a figure number, and create a filled contour plot of the 
    variable data array with respect to the depths and distance arrays given
    '''
    plt.figure(num = figNum)
    plt.gca().invert_yaxis()
    plt.contourf(dist2d, depths2d, var2d, cmp=plt.cm.coolwarm)
    plt.colorbar()
    plt.show()

def main():
    # Create and put all CNV datasets into a list
    prof1 = fCNV("prof1.cnv")
    prof2 = fCNV("prof2.cnv")
    prof3 = fCNV("prof3.cnv")
    profList = [prof1, prof2, prof3]
    
    # Create numpy array of new depths
    newDepth = np.arange(5, 101, 5)
    
    # Create 2-D array of data salinity and temperature using createTwoDimData function
    tempData = twoDemData('TEMP', profList, newDepth)
    salData = twoDemData('PSAL', profList, newDepth)

    # Get lat and longitude arrays, data from each cnv dataset
    lats = [prof1.attrs['LATITUDE'], prof2.attrs['LATITUDE'], prof3.attrs['LATITUDE']]
    lons = [prof1.attrs['LONGITUDE'], prof2.attrs['LONGITUDE'], prof3.attrs['LONGITUDE']]
    
    # Get list of distances between datapoints
    dists = calcDistance(lats, lons)
    
    # Create grid of depths and distances
    depths_2d, dists_2d = np.meshgrid(newDepth, dists)
    
    # Plot cross section of depths, distances and data on separate figures
    plotCrossSection(depths_2d, dists_2d, tempData, 0)
    plotCrossSection(depths_2d, dists_2d, salData, 1)


if __name__ == "__main__":
    main()