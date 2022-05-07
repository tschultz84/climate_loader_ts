#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 14:35:20 2020

@author: scs-rd
"""

import numpy as np
import scipy as sc
from scipy import stats
import matplotlib.pyplot as plt
import math
import array as arr
import pandas as pd
import sys
import scipy.constants
import descartes
import geopandas as gpd
from shapely.geometry import Polygon, Point, MultiPolygon
from matplotlib import pyplot

#import shapefile


#%%
#Loading the station data



#This laods station data from .dly files
#stationset1 is just a list of IDs of stations
#elements are the elements extracted, in an array, like TMIN, PRCP, etc.


class StationData :
    def __init__(self,stationset1,elements=["TMAX","TMIN"],years=np.arange(1880,2020),display=True,dropflags=True):
       
        self.flags=dropflags
        
        self.displayer1 = display
        self.keepels = elements
        
        #This flag just prompts the calcualtion of average temperature
        #Defaults to false
        self.tavgcalc= False
        
        if "TMEAN" in elements:
            self.keepels = ["TMAX","TMIN"]
            self.tavgcalc = True
            
            print("You requested the calculation of the average temperature.")
            print("This is the average of TMAX and TMIN for each day in the dataset.")
            print("Only TMAX and TMIN are being loaded. The output will be calculated based on this data.")
        
        
        self.keepyears = years
        if (type(elements) != list)   :
            sys.exit("Please specify a correct format for elements which are to be retained. it should be a list of strings.")
            
        
        #This below records only the unique stations, removing any redundancies
        #self.stationset = np.asarray(stationset1['ID'].unique())
        self.stationset = np.asarray(stationset1['ID'].unique())
        #self.stationset = np.asarray(stationset1[['ID']]).reshape([-1])
        self.shortstationset = self.stationset
        #self.datapath = '/Users/scs-rd/Desktop/OneDrive-Tobias/OneDrive - SCS Global Services/Reference/Temperature-Datasets/ghcnd_all/ghcnd_all/'
        self.datapath = '/Users/scs-rd/Desktop/Tobias/Github/climate/Datasets/ghcnd_all/ghcnd_all/'
        
        #DATAPATH is just where all of the information is located.
        #self.stDATAPATH = '/Users/scs-rd/Desktop/OneDrive-Tobias/OneDrive - SCS Global Services/Reference/Temperature-Datasets/'
        self.stDATAPATH = '/Users/scs-rd/Desktop/Tobias/Github/climate/Datasets/'
        
        
        #This loads all of teh station meta information, like lat/lon.
        self.stationlist = pd.read_csv(self.stDATAPATH+'ghcnd-stations.txt',sep='\t')
        
        self.stationdata = None #This is initalized as None before being populated.
        
        print("Initializing "+str(elements)+" from "+str(len(self.shortstationset))+ " stations:")
        print(self.stationset)
        
        
        if display:
        
            #This retrioeves all lat/lon information
            latlonmat = []
            for i in range(0,len(self.stationset)):
                inter = self.stationlist[self.stationlist["ID"]==self.stationset[i]][["LATITUDE",'LONGITUDE']]          
                latlonmat.append(np.asarray(inter)[0])
            
            self.lats = np.asarray(latlonmat)[:,0]
            self.longs = np.asarray(latlonmat)[:,1]
                    
            self.showonmap(self.stationset)
           
        
    #This displays a set of stations, in a array containing their ID,
    #on the world map  
    #stations_shown is a list of stations, by ID number
    #if  nont a number, it defautls to all stations in thie object
    def showonmap(self,stations_shown=np.nan):
        #This loads the world map, necessary to display station locations
        world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
        if type(stations_shown) == float:
            station = self.stationset
        if type(stations_shown) == np.ndarray:
            station=stations_shown

        #This shows the global scale map. 
        with plt.style.context(("seaborn", "ggplot")):
            world.plot(figsize=(18,10),
                       color="white",
                       edgecolor = "grey");
        
            
            plt.scatter(self.longs,self.lats, s=5, color="red", alpha=0.3)
            plt.xlabel("Longitude")
            plt.ylabel("Latitude")
        #This shows a map constrained just to the region where the stations re present.     
        with plt.style.context(("seaborn", "ggplot")):
            world.plot(figsize=(18,10),
                       color="white",
                       edgecolor = "grey");
        
            
            plt.scatter(self.longs,self.lats, s=5, color="red", alpha=0.3)
            plt.xlabel("Longitude")
            plt.ylabel("Latitude")
            plt.xlim([np.min(self.longs),np.max(self.longs)])
            plt.ylim([np.min(self.lats),np.max(self.lats)])
    
    #This loads the data for the stations listed in the IDs.
    #randomsize is the number of stations about which data is loaded (if a random sample is desired).
#this is saved into self.statinodata for oeasy reference by other functions.         
    def load_stats(self,randomsize=None):
        #This resets shortstationset, in case  it has been previously defined
        self.shortstationset = self.stationset
        #If random size is not specified, then all stations are loaded.
        if randomsize != None:
            self.shortstationset = np.random.choice(self.stationset,int(randomsize))
        dropflags=self.flags    
        #If dropflags=True, then all the flag columns are dropped.
        #If False, then everything is loaded.

        #This is the reader function for the temperature station data, which is organized
        # in a very idiosyncratic fashion. It generates the column number dividers 
        #needed to read the file.
        #See https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/readme.txt
        
        #The following are the column dividers and names, pre-loaded
        #The difference in the reading algorithm is how many columns are read in
        if dropflags == True:
            self.colspecs1=[(0, 11),(11, 15), (15, 17),  (17, 21),  (21, 26),  (29, 34),(37, 42), (45,50),  (53, 58),
                   (61, 66),(69, 74),(77, 82),  (85, 90), (93, 98),  (101, 106),  (109, 114),  
                     (117, 122), (125, 130),(133, 138), (141, 146), (149, 154), (157, 162), 
                     (165, 170),(173, 178),(181, 186),(189, 194), (197, 202), (205, 210), 
                    (213, 218),(221, 226),  (229, 234), (237, 242), (245, 250),(253, 258), (261, 266)]
            self.columnnames=np.array(['Station_ID', 'Year', 'Month', 'Element', 'Day1', 
                                        'Day2',  'Day3','Day4','Day5','Day6','Day7', 'Day8',  'Day9', 
                                        'Day10', 'Day11',  'Day12', 'Day13', 'Day14',  'Day15',  
                                        'Day16','Day17', 'Day18',  'Day19','Day20',  'Day21',
                                         'Day22','Day23', 'Day24', 'Day25','Day26', 'Day27', 
                                         'Day28','Day29','Day30','Day31']) 
        if dropflags != True:            
            self.colspecs1 = [(0, 11), (11, 15), (15, 17), (17, 21), (21, 26), (26, 27), (27, 28),
                        (28, 29), (29, 34), (34, 35), (35, 36), (36, 37), (37, 42), (42, 43),
                        (43, 44), (44, 45), (45, 50), (50, 51), (51, 52), (52, 53), (53, 58),
                        (58, 59), (59, 60), (60, 61), (61, 66), (66, 67), (67, 68), (68, 69),
                        (69, 74), (74, 75), (75, 76), (76, 77), (77, 82), (82, 83), (83, 84),
                        (84, 85), (85, 90), (90, 91), (91, 92), (92, 93), (93, 98), (98, 99),
                        (99, 100), (100, 101), (101, 106), (106, 107), (107, 108), (108, 109),
                        (109, 114), (114, 115), (115, 116), (116, 117), (117, 122),
                        (122, 123),  (123, 124),  (124, 125),  (125, 130),  (130, 131),  (131, 132),  (132, 133),
                        (133, 138),  (138, 139),  (139, 140),  (140, 141),  (141, 146),  (146, 147),
                        (147, 148),  (148, 149),  (149, 154),  (154, 155),  (155, 156),  (156, 157),
                        (157, 162),  (162, 163),  (163, 164),  (164, 165),  (165, 170),  (170, 171),
                        (171, 172),  (172, 173),  (173, 178),  (178, 179),  (179, 180),  (180, 181),  (181, 186),
                        (186, 187),  (187, 188),  (188, 189),  (189, 194),  (194, 195),  (195, 196),
                        (196, 197),  (197, 202),  (202, 203),  (203, 204),  (204, 205),  (205, 210),  (210, 211),  (211, 212),  (212, 213),  (213, 218),  (218, 219),  (219, 220),  (220, 221),
                        (221, 226),  (226, 227),  (227, 228),  (228, 229),  (229, 234),  (234, 235),  (235, 236),  (236, 237),  (237, 242),  (242, 243),  (243, 244),  (244, 245),  (245, 250),  (250, 251),
                        (251, 252),  (252, 253),  (253, 258),  (258, 259),  (259, 260),  (260, 261),  (261, 266),  (266, 267),  (267, 268),  (268, 269)]

            self.columnnames = np.array(['Station_ID', 'Year', 'Month', 'Element', 'Day1', 'MFLAG1',
               'QFLAG1', 'SFLAG1', 'Day2', 'MFLAG2', 'QFLAG2', 'SFLAG2', 'Day3',
               'MFLAG3', 'QFLAG3', 'SFLAG3', 'Day4', 'MFLAG4', 'QFLAG4', 'SFLAG4',
               'Day5', 'MFLAG5', 'QFLAG5', 'SFLAG5', 'Day6', 'MFLAG6', 'QFLAG6',
               'SFLAG6', 'Day7', 'MFLAG7', 'QFLAG7', 'SFLAG7', 'Day8', 'MFLAG8',
               'QFLAG8', 'SFLAG8', 'Day9', 'MFLAG9', 'QFLAG9', 'SFLAG9', 'Day10',
               'MFLAG10', 'QFLAG10', 'SFLAG10', 'Day11', 'MFLAG11', 'QFLAG11',
               'SFLAG11', 'Day12', 'MFLAG12', 'QFLAG12', 'SFLAG12', 'Day13',
               'MFLAG13', 'QFLAG13', 'SFLAG13', 'Day14', 'MFLAG14', 'QFLAG14',
               'SFLAG14', 'Day15', 'MFLAG15', 'QFLAG15', 'SFLAG15', 'Day16',
               'MFLAG16', 'QFLAG16', 'SFLAG16', 'Day17', 'MFLAG17', 'QFLAG17',
               'SFLAG17', 'Day18', 'MFLAG18', 'QFLAG18', 'SFLAG18', 'Day19',
               'MFLAG19', 'QFLAG19', 'SFLAG19', 'Day20', 'MFLAG20', 'QFLAG20',
               'SFLAG20', 'Day21', 'MFLAG21', 'QFLAG21', 'SFLAG21', 'Day22',
               'MFLAG22', 'QFLAG22', 'SFLAG22', 'Day23', 'MFLAG23', 'QFLAG23',
               'SFLAG23', 'Day24', 'MFLAG24', 'QFLAG24', 'SFLAG24', 'Day25',
               'MFLAG25', 'QFLAG25', 'SFLAG25', 'Day26', 'MFLAG26', 'QFLAG26',
               'SFLAG26', 'Day27', 'MFLAG27', 'QFLAG27', 'SFLAG27', 'Day28',
               'MFLAG28', 'QFLAG28', 'SFLAG28', 'Day29', 'MFLAG29', 'QFLAG29',
               'SFLAG29', 'Day30', 'MFLAG30', 'QFLAG30', 'SFLAG30', 'Day31',
               'MFLAG31', 'QFLAG31', 'SFLAG31'])

        
        self.nostations = len(self.shortstationset)
        print("Loading data  from "+str(len(self.shortstationset))+ " stations.")
        print("Loading in process......")
        
        #This defines the pre-loading data, i.e. which columns are used to skip rows.
        shortcols = [ (11, 15), (17, 21)]
        shortnames =  np.array([ 'Year', 'Element']) 

        #counter is just used to help dispaly the number of stations so far successfully loaded
        counter = 0
        
        dayscount=0
        monthscount=0
        
        #This initializes the dataframe which will be returned. 
        stationer1  = pd.DataFrame(columns=self.columnnames)


        for j in np.arange(self.nostations):
            #Where the station data files are
            path = self.datapath+str(self.shortstationset[j])+'.dly'
            
            #This reads in the column which is defined such to skip rows .
            prelim = pd.DataFrame()    
            prelim=pd.read_fwf(path,colspecs=shortcols)
            prelim.columns=shortnames
            
            #This generates the rows to skip.   
            #skiprows1 = prelim[~prelim['Element'].isin(self.keepels)].index
            skiprows1 = prelim[(~prelim['Element'].isin(self.keepels)) | (~prelim['Year'].isin(self.keepyears))].index
            
            
            inter0 = pd.DataFrame()      
            #This then Reads out the next file.  
            inter0 = pd.read_fwf(path,colspecs=self.colspecs1,skiprows=skiprows1+1)
            inter0.columns=self.columnnames
            
            #here, if the taverage calculation is initialized, it doe the average
            if self.tavgcalc == True:
                
                #This averages across TMAX and TMIN values for each day.
                inter0 = inter0.groupby(["Year","Month"], as_index=False).mean()
                inter0['Element']='TMEAN'
                inter0["Station_ID"]=self.shortstationset[j]
                #This extends the output dataframe.
                stationer1 = stationer1.append(inter0)
             
            #OR, if not average temperature is calculated, jsut reads out the data    
            if self.tavgcalc != True:
                #This extends the output dataframe.

                stationer1 = stationer1.append(inter0)
            
        
            
            #The following lines just print to screen the progress. 
            counter = counter +1
            #This increments by adding the days worth of data.
            dayscount = dayscount+inter0.count().iloc[4:].sum()
            
            #and the number of months
            monthscount = monthscount+inter0.shape[0]
            if len(self.shortstationset)>=500:
                resetter = int(len(self.shortstationset)/50)
            if len(self.shortstationset)<500:  
                resetter=25
            if counter == resetter:
                
                print(str(j)+" stations loaded...Loading still in process, "+str(int(100*j/len(self.shortstationset)))+"% completed.")
                counter=0
           
             
        #this replaces the no data values with  not a number  values   
        #The value -9999 denotes no data fields. This allows the
         #function which takes the mean to work properly. 
        
        stationer1 = stationer1.replace(-9999,np.nan) 
        self.stationdata=stationer1    
        
        liststations = self.stationdata['Station_ID'].unique()
        print("I loaded data from these  "+str(len(liststations))+ " stations:")
        print(liststations)
        print("It includes "+str(monthscount)+" months of data")
        print(" and "+str(dayscount)+" days of data across all elements.")

        print("Not all stations were necessarily loaded, because not all have the ")
        print(" elements which are specified.")
        
        #This provides a summary of the data, by element, if 
        #display is set to True.
        if self.displayer1:
            eldat = pd.DataFrame(columns=['Year']+self.keepels)
            
            
            for j in np.arange(0,len(self.keepyears)):
                yearint = self.keepyears[j]
                inter0 = stationer1[stationer1['Year']==yearint]
                #Initializes the dictionary which is used to create the Pandas Dataframe
                h5 = {}
                h5["Year"]=yearint
                emp=[]
                for i in np.arange(0,len(self.keepels)):
                    inter = 0
                    #IT sets different parameters, depending on if taverage is calculated
                    if self.tavgcalc ==True:
                        ele0 = "TMEAN"
                    if self.tavgcalc !=True:
                        ele0 = self.keepels[i]
                    #This counts all values with numbers for this element, for this year.
                    inter = inter0[inter0['Element']==ele0].count().iloc[4:].sum()
                    h5[str(ele0)]=inter
                    
                    
                eldat=eldat.append(h5,ignore_index=True)
            eldat=eldat.sort_values('Year')
            self.filldata = eldat           
            eldat.plot(x='Year',kind='line',
                       title='Number of days per year of elment day (total for '+str(self.nostations)+' stations)')
                       
        return stationer1
    

#%%   
#All analysis functions for station data are in this class.
#It is passed an object that is the output of the StationData class. 
#The argumetns are : mind is the minimum number of days needed in every nmonth of a 
#dataset to include; byears and ayears are the baseline and assessment period years  
#ele, if not jsut a stadnard operato, is one of the followign:
#HEAT : Selects TMAX from the summer months (analyzing summer heat waves)
#COLD : Selects TMIN from teht winter months (analyzing winter cold snaps)
#TMID : Measures the average of TMAX and TMIN for every day in the dataset           
class AnalyzeStationData :
    def __init__(self,stationdataobj,mind=25,ele="TMAX",byears=[1880,1920],ayears=[2000,2020],noby=20,noay=10,show=False):
        #Initializes entire object.
        self.stobj = stationdataobj
        
        #This field determines if varioous charts and printed features are actually displayed
        self.displayon=show
        
        print("Loading in data from "+str(self.stobj.stationdata['Station_ID'].unique().shape[0])+" stations.")
        
        
        self.element=ele
        
        #This initiatilizes a temporary variable loading stationdata
        stdataint = self.stobj.stationdata
        
        #Strips out unused data.
        self.stdata=stdataint[stdataint['Element']==self.element]
        
        #Removes years greater than or less than the earliest/latest baseline/assessment period years
        self.stdata=self.stdata[self.stdata['Year']>=byears[0]]
        self.stdata=self.stdata[self.stdata['Year']<=ayears[1]]
        #And strip out intermediate years
        self.stdata=self.stdata[(self.stdata['Year']<=byears[1]) | (self.stdata['Year']>=ayears[0])]
        
        self.stationlist = stdataint['Station_ID'].unique()
        
        if ele=="HEAT":
            stdataint = stdataint[stdataint['Element']=="TMAX"]
            stdataint = stdataint[stdataint['Month'].isin([6,7,8])]
        
        if ele=="COLD":
            stdataint = stdataint[stdataint['Element']=="TMIN"]
            stdataint = stdataint[stdataint['Month'].isin([12,1,2])]
            

        if (ele != "SUMMER") and (ele!="TMID") and (ele!="COLD"):
        
            stdataint = stdataint[stdataint['Element']==self.element]
        
            
        self.basey = byears
        self.assey = ayears
        #mindays is the minimum number of days needed in every mnonth to be lincuded.
        self.mindays = mind
        
        
        
        #This specifies the minumum number of years  needed in t he b aseline and assessment.
        self.minbyears = noby
        self.minayears = noay
        
        print("Have loaded in station data containing  "+str(self.stdata.shape[0])+" rows.")
        print("The requirements for assessment are as follows.")
        print("Every month of each year must have "+str(mind)+" days of data.")
        print("The baseline period includes "+str(byears[1]-byears[0])+" years, of which "+str(noby)+" years must have data.")
        print("The baseline period runs from "+str(byears[0])+" to "+str(byears[1]))
        print("The assessment period includes "+str(ayears[1]-ayears[0])+" years, of which "+str(noay)+" years must have data.")
        print("The assessment period runs from "+str(ayears[0])+" to "+str(ayears[1]))
   
    #This checks if the year is complete.
    #It must have all months.
    #And the number of days in eveyr month must be higher than self.mmindays
    def check_year_complete(self,stationID,year):
        check = False
        
        #Loads the data for the selected station. 
        data = self.stdata[self.stdata['Station_ID']==stationID]
        
        data1 = data[data['Year']==year]
        months = np.sort(np.asarray(data1['Month']))
        #Checks if all months are present. If a month is missing, the dataset is not used.
        if np.alltrue(months == np.array([ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12])):
            #This counts the values in each row. It returns a dataframe  containing a list of numbers
            #First, drop extraneous columns.
            data5=data1.drop(["Station_ID","Element","Year","Month"],axis=1)
            #Then, count the number of cells which are NOT NaN using .count
            #Then see if the number in each row is greatter tha the set minimum number of days
            count=np.asarray(data5.count(1).ge(self.mindays))
            if np.alltrue(count):
                check = True
        return check
            
   #This, for a single station, lists the eyars which re compelte.
    def list_complete_years(self,station):
     #Loads the data for the selected station. 
        data = self.stdata[self.stdata['Station_ID']==station]
        #Generates list of the years in the baseline and assessment range. 
        listerb = np.arange(self.basey[0],self.basey[1]+1)
        listera = np.arange(self.assey[0],self.assey[1]+1)
        
        yearsb = []
        yearsa = []
        
        #This removes from the list,  years which do not have sufficient number of days from the baseline.
        for i in np.arange(0,len(listerb)):
            if (self.check_year_complete(station,listerb[i])):
                yearsb.append(listerb[i])
        
        #This removes from the list, years which do not have sufficient number of days from the assessment.
        for i in np.arange(0,len(listera)):
            if (self.check_year_complete(station,listera[i])):
                yearsa.append(listera[i])
                
        return [yearsb,yearsa]
    
    #This lists the stations which have a complete record. 
    def list_useful_stations(self,printer=True):
        
        #Lists all stations. 
        allstations=self.stdata["Station_ID"].unique()
        nostations = allstations.shape[0]
        
        print("Reviewing usefulness of the "+str(nostations)+" stations in the data.")

        returner = []
        for i in np.arange(0,nostations):

            anyb = self.list_complete_years(allstations[i])[0]
            anya = self.list_complete_years(allstations[i])[1]
            
            if (anyb != []) and (anya != []):
                returner.append(allstations[i])
        if printer == True:
            print("Here is the list of useful stations :")
            print(returner)
        return   returner
 
    #This grabs data for a specific period of time from a set of stations.    
    #stations is a station ID in the dataset object to be included
    #baseline and assyears are the baseline range and assessment range of years
    #element1 is the element to be evaluated, which can be different from the object
    def perioddata(self,station):
                
         #Loads the data for the selected station. 
        data = self.stdata[self.stdata['Station_ID']==station]

        #Thsi function p ulls the years with complete data for this station.
        yearsl = self.list_complete_years(station)
        yearsb = yearsl[0]
        yearsa = yearsl[1]

        #It turns only data for this years where enough data is present.         
        base = data[data["Year"].isin(np.asarray(yearsb))]
        asse = data[data["Year"].isin(yearsa)]
        
        return base,asse
  
    
      #This puts out the data in a format which is sutiable for a historgram. 
    def data_frequency(self,station):
        #This generates a lsit of eyars which have adequate data. 
        based1 , assed1 = self.perioddata(station)
        proceed = True
        
        printer=self.displayon
        
        
        based , assed = based1, assed1


        #This lines drop unecessary columsn. 
        
        based = based.drop(["Station_ID","Element","Year","Month"],axis=1)
        assed = assed.drop(["Station_ID","Element","Year","Month"],axis=1)
        
        #Then flatten. 
        based=np.asarray(based.values.tolist())
        based = based.flatten()
        
        assed=np.asarray(assed.values.tolist())
        assed = assed.flatten()

        
        #This modifies the variable if its' in temperature
        if self.element == "TMIN" or "TMAX" or "TAVG":
        
            based = based/10*9/5+32
            assed = assed/10*9/5+32
            label = "Farenheit Degrees"
        
        if printer == True:            
            print("---------")
            print("Element shown: "+str(self.element))
            print("Station: "+str(station))
            print("Baseline Years: "+str(based1['Year'].unique()))
            print("Assessment Years: "+str(assed1['Year'].unique()))
                
        return [based,assed]
    
        #this does an anlaysis, creating a statistical compendium of
    #differenes in the baseline period and tthe assessment period for the station
    #  noted. 
    #It  rerturns a dataframe summarzigin statistical summary
    def stats(self,stations):
        
        printer=self.displayon
        
        #This generates a lsit of eyars which have adequate data. 
        if printer==True:
            print("Analyzing statistical data for "+str(stations))
        information = self.data_frequency(stations)
        based = information[0]
        assed = information[1]
        
        
        #DFirst, check if there's enough data. 
        if (based.size!=0) and (assed.size!=0):
            proceed = True
        if (based.size==0) or (assed.size==0):    
            proceed= False
        
        if proceed ==False:
            if self.displayon==True:
                print("Error! Inadequate years in the baseline or assessment.")
                
                print("This is for: "+str(self.element))
                print("Inadequate data for Station : "+str(stations))

        
        #ONLY PROCEED IF THERE IS.
        if proceed == True:
        
             #Removes nan, so that statistical values can be assessed
            based = based[np.logical_not(np.isnan(based))]  
            assed = assed[np.logical_not(np.isnan(assed))] 
             
            #This evaluates key statistical metrics
            #Average
            basedmean = np.average(based)
            assmean = np.average(assed)
            #T-test
            tscore,pvalue = sc.stats.ttest_ind(based,assed,equal_var=False,nan_policy='omit')
            
            #Median
            basedmedian = np.median(based)
            assmedian = np.median(assed)
            #Standard error
            basedstd = np.std(based)
            asstd = np.median(assed)
            #85 percentile
            based85 = np.percentile(based,85)
            assed85 = np.percentile(assed,85)
            #15 percentile
            based15 = np.percentile(based,15)
            assed15 = np.percentile(assed,15)
            #number of days in sample
            basedno=len(based)
            assedno=len(assed)
            
            #This defines the return data
                
            dataret =  {"Station_ID":[stations],
                        "Element":[self.element],
                        "Mean, Baseline":[basedmean],
                        "Mean, Assessment":[assmean],
                        
                        "Median, Baseline":[basedmedian],
                        "Median, Assessment":[assmedian],
                        
                        "Std Dev, Baseline":[basedstd],
                        "Std Dev, Assessment":[asstd],
                        
                        "No. of Days in Sample, Baseline":[basedno],
                        "No. of Days in Sample, Assessment":[assedno],
                        
                        "85% Percentile, Baseline":[based85],
                        "85% Percentile, Assessment":[assed85],
                        
                        "15% Percentile, Baseline":[based15],
                        "15% Percentile, Assessment":[assed15]
                        
                        
                        }
            
            if self.element == "TMIN" or "TMAX" or "TAVG"  or "TMEAN":
                #Percent Days over 95
                basedhot = np.count_nonzero(based >= 95)/len(based)
                assedhot = np.count_nonzero(assed >= 95)/len(assed)
                #Percent Days beneath 20
                basedcold = np.count_nonzero(based <= 20)/len(based)
                assedcold = np.count_nonzero(assed <= 20)/len(assed)
                
                dataret['Fraction of Days >95, Baseline']=[basedhot]
                dataret['Fraction of Days >95, Assessment']=[assedhot]
                
                dataret['Fraction of Days <20, Baseline']=[basedcold]
                dataret['Fraction of <20, Assessment']=[assedcold]
            
            if printer == True:
                
                print("Statistical Summary:")
                
                print("Baseline number of days:"+str(len(based)))
                print("Assessment number of days:"+str(len(assed)))
                
                print("Baseline Average: "+str(basedmean))
                print("Assessment Average: "+str(assmean))
                
                print("T-Score of means:"+str(tscore))
                print("P-Value of difference:"+str(pvalue))
                print("that is, you would expect to see a difference between the sample means ")
                print("as extreme as "+str(basedmean-assmean)+" only "+str(pvalue*100)+"% of the time")
                
                
                print("Baseline Median: "+str(basedmedian))
                print("Assessment Median: "+str(assmedian))
            
                print("Baseline Standard Error: "+str(basedstd))
                print("Assessment Standard Error: "+str(asstd))
            
                print("Baseline 85% Percentile : "+str(based85))
                print("Assessment 85% Percentile : "+str(assed85))            
                
                        
                print("Baseline 15% Percentile : "+str(based15))
                print("Assessment 15% Percentile : "+str(assed15))
                if self.element == "TMIN" or "TMAX" or "TAVG":
                    print("Baseline Percent days over 95 degrees : "+str(basedhot))
                    print("Assessment Percent days over 95 degrees : "+str(assedhot))
                    print("Baseline Percent days beneath 20 degrees : "+str(basedcold))
                    print("Assessment Percent days beneath 20 degrees : "+str(assedcold))
                print("Analysis complete for "+str(stations)+"!")
                
            #Returns a Dataframe
            return pd.DataFrame(dataret)

        

        
        
    
    #Creates a pandas dataframe containing all summary datapoints for all statioons
    def all_stats(self):
        allstations=self.stationlist
        returner=pd.DataFrame()
        
        account=0
        counter=0
        for i in range(0,len(allstations)):
            
            returner=returner.append(self.stats(allstations[i]))
            counter=counter+1
            if counter>25:
                account=account+25
                counter=0
                print("Analyzed "+str(account)+" stations so far.")
                print("That is "+str(int(100*account/len(allstations)))+"% complete.")
                print(str(len(allstations)-account)+" remain.")
        return returner
    
    #this does an anlaysis, creating a statistical compendium of
    #differenes in the baseline period and tthe assessment period for the station
    #  noted. 
    def compdata(self,stations,hist=True,num_bins=20):
        #This generates a lsit of eyars which have adequate data. 
     
        information = self.data_frequency(stations)
        based = information[0]
        assed = information[1]
        
        
        #DFirst, check if there's enough data to even show a histogram. 
        if (based.size!=0) and (assed.size!=0):
            proceed = True
        if (based.size==0) or (assed.size==0):    
            proceed= False
            
        #ONLY PROCEED IF THERE IS.
        if proceed == True:
        
            if self.element == "TMIN" or "TMAX" or "TAVG":
                
                label = "Farenheit Degrees"
            
            if hist == True:
                pyplot.hist(based, bins=num_bins, alpha=0.5, label='Baseline',normed=1)
                pyplot.hist(assed, bins=num_bins, alpha=0.5, label='Assessment',normed=1)
                pyplot.legend(loc='upper right')
                pyplot.ylabel('Frequency of Data Points')
                pyplot.xlabel(label)
                pyplot.title('Climate Variable Distribution')
                pyplot.show()
                
            self.stats(stations)
            
            return based,assed
        
        if proceed ==False:
            print("Error! Inadequate years in the baseline or assessment.")
            
            print("This is for: "+str(self.element))
            print("Station not shown: "+str(stations))

        print("---------")

  
