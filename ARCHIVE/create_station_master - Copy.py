#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 07:11:53 2020

@author: scs-rd
"""
#This cell imports the right packages.

import pandas as pd


#DATAPATH is just where ghcnd-stations and ghcnd-inventory are located on your hard drive.
stDATAPATH = "C:\\ts_big_data_files\\"

#This loads the file containing all of the station meta information, including lat/lon.
allstations = pd.read_fwf(stDATAPATH+'ghcnd-stations.txt',colspecs=[
                                                                          #ID#
                                                                        (0,11),
                                                                        #Lat/Lon
                                                                        (11,20),(21,30),
                                                                        #Eelvation, tength of meters
                                                                        (31,37),
                                                                        #Station name
                                                                        (38,72),
                                                                        #Other
                                                                        (73,82),(83,85)
                                                                        ])
allstations.columns=['ID','Lat','Lon','Elevation (0.1m)','Name','Other1','Other2']
#%%
#This pulls the list of countries. THe first two eltters of each weather station ID
# is the country.
countrylist = pd.read_fwf(stDATAPATH+'ghcnd-countries.txt',colspecs=[(0,3),(3,49)])
countrylist.columns = ['Code','Name']
#%%
#allstations = pd.read_csv(stDATAPATH+'ghcnd-stations.txt',sep='\t')
#allstations.columns=['ID','Latitude','Longitude','Elevation','Statname']

#This loads the datafile containing the time periods of analysis
stationperiods = pd.read_fwf(stDATAPATH+'ghcnd-inventory.txt',colspecs=[
                                                                         #Country code#
                                                                        (0,2),                                                                       
                                                                        #ID#
                                                                        (0,11),
                                                                        #Lat/Lon
                                                                        (11,20),(21,30),
                                                                        #Others
                                                                        (31,35),
                                                                        #First year, Last Year
                                                                        (36,40),(41,45)
                                                                        ])
stationperiods.columns=['Country','ID','Lat','Lon','Element','Firstyear','Lastyear']

#%%
#This merges the data into one DataFrame containing all station related data, including all columns, and all elements, smashed together into one enormous
#dataframe (witih over 700,000 roles). 
allstdata = pd.merge(allstations,stationperiods)

#%%
#tempdata=allstdata[allstdata['Element'].isin(["TAVG","TMAX","TMIN"])]
#This creates a dataframe that focus on the columns and variables of interest. 
sttempdata=allstdata[allstdata['Element'].isin(["TAVG","TMAX","TMIN","PRCP","SNOW"])]
#STILL NEED TO DROP SOME COLUMSN.
sttempdata=sttempdata[['ID','Lat', 'Lon', 'Element', 'Firstyear', 'Lastyear']]