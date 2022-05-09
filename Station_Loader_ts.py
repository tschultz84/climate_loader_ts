# -*- coding: utf-8 -*-


import pandas as pd
from math import radians, cos, sin, asin, sqrt
import time
import sys
#import yaml
import numpy as np
from datetime import date
import requests
import warnings
warnings.filterwarnings("ignore")


"""
This object has functions to identify the station nearest to a reference lat, lon.
It then downloads the data from the climate station closest to this lat lon.
It then runs several checks on this data to see if it is acceptably complete, which you can infer from the input variables. 
It also calculates the TMID values, which is the average of TMAX and TMIN in every day.

#The outputs of interest are as follows:
    self.refpoint is the reference
    self.station_information gives meta information, like the station name, ID, lat/lon, etc.
    self.station_filters summarizes all the filters you entered as inputs. 
    self.closest_stations gives a list of the  nearest stations
    self.station_data is a dataframe ahving loaded all the statin's data, a numpy array 
     with these columns: [Year, Month, Day, Day of Year, TMAX, TMIN, TMID]

    The inputs are as follows:
    *point is a list, in the form [latitude,longitude], which is the reference point 
        (data will be loaded for the station closest to these coordinates where enough data
        is available.)
    *printupdate, is only used to display status updates for the script.If false, no updates are printed
    *min_days_per_mo is the minimum number of dasys of data which must be present for every month of a year for 
        a dataset to be included. If not, then the year is excluded from dataset.
    *search_radius is the lat/lon degrees around the reference piont that are searched. Defaults to 1, meaning 
        a searched is performed within 1 degree north, south, east, and west, of point.
    *firstyear is the earliest year which must be present in the dataset record for the weather station
        to be loaded. (Note: The station must have data for the 2nd most recent year, i.e., in 2022, for 2021.)
    *basenoyears is the number of years in the baseline period.
    *min_recent_years is the minimum number of years before he present which must be present
        for the statin to be loaded.
    *required_trend_years minimum number of years in the last 30 required to calculate this trend.
        Weather stations without this many years present in the last 30 years of data will be discarded.
    *lastbaseyear is the last year in which the baseline is calculated (not the same as last year)
    
"""    
class LoadStation :
    def __init__(self,point,printudpate=False,min_days_per_mo=15,search_radius=1,firstyear=1890,basenoyears=30,
            min_recent_years=5,required_trend_years =20,lastbaseyear=1955):
        #Basic initializations - define your print variable and ensure the 'point' variable is in the right format. 
        self.display=  printudpate   
        #A couple of checks. 
        #If a point (lat, lon) is passed, then create a point. 
        if (type(point) != str) and (len(point)==2): self.refpoint=point
        if required_trend_years >=30: required_trend_years = 30 #This by definition should not be more than 30. 
        #Initiate all the variables which are used later on, and load them into a Series for easy review.
        self.station_filters=pd.Series(data={
            "Reference Point Latitude":point[0],
            "Reference Point Longitude":point[1],
            "Earliest Year Required in Dataset":firstyear,
            "Min Days in Every Month":min_days_per_mo,
            "Search Radius (deg)":search_radius,
            "Min Number of Years before Present Year":min_recent_years,
            "Min Years in Last 30 Years for Trend":required_trend_years,
            "Last Year of Baseline":lastbaseyear,
            "Required Number of years in Baseline Period":basenoyears    })
        
        #Initialize the station_information Series with nan values, which are replaced if you find a station. 
        self.station_information = pd.Series(data={
                'Weather Station Name':np.nan, #Note station Name.
                'Station ID':np.nan, #Note station ID.
                'Earliest year in station record':np.nan,
                'Most recent year in record':np.nan,
                'Number of years in record meeting quality requirements':np.nan,
                'Station Latitude':np.nan, #Define station coordinates.
                'Station Longitude':np.nan,
                "Reference Point Latitude":point[0], #Add in the original reference point for reference.
                "Reference Point Longitude":point[1],
                'Miles from Reference Point':np.nan
                })

        #Files to load.
        self.STATION_META_FILENAME = "data\\ghcnd_station_master_ts_tmax.csv" #containing all the lat / lon information for the stations.
        self.NOAA_URL = 'https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/' #url to the directory containing the DLY files. 

        #First, the list of stations is generated. This function finds all the stations within search_radius nearest to the point, without applying any filters.
        self.closest_stations = self.nearest_station(point)
      
        #Now, all of the stations are checked for conformance to the filters. 
        #Add two flag variables. 
        keep_going = True # add a flag variable to keep the loop running. 
        isgood = False #Initialize this flag variable if a station is found to satisfy all filters.
        for index,station in self.closest_stations.iloc[1:].iterrows():
            if keep_going == True:
                self.run_this_baby(station['ID']) #Run the download from the NOAA URL. This defines self.station_data for this station.
                isgood = self.StationDataCheck(self.station_data) #Now, run a check on the downloaded station data. 
                #This flag stops the loop, and assigns all variables, if checks are complete.
                if isgood == True: 
                    keep_going=False #Kill the loop. 
                    #Output all the station meta information as a pandas Series. 
                    self.station_information = pd.Series(data={
                        'Weather Station Name':station['Name'], #Note station Name.
                        'Station ID':station['ID'], #Note station ID.
                        'Earliest year in station record':np.unique(self.station_data[:,0])[0],
                        'Most recent year in record':np.unique(self.station_data[:,0])[-1],
                        'Number of years in record meeting quality requirements':np.unique(self.station_data[:,0]).shape[0],
                        'Station Latitude':station['Latitude'], #Define station coordinates.
                        'Station Longitude':station['Longitude'],
                        "Reference Point Latitude":point[0], #Add in the original reference point for reference.
                        "Reference Point Longitude":point[1],
                        'Miles from Reference Point':station['Miles_from_Ref']
                        })
                    self.st_latlon_str = f"{round(station['Latitude'],2)} latitude {round(station['Longitude'],2)} longitude" #Lat Longitude string, used in an output table.
                                                         
                    
                    print(f"Station ID# {station['ID']}, called {station['Name']} is complete. It's good to use.")
                    print(f"This station is {self.station_information['Miles from Reference Point']} miles from the reference point.")
                        
                #If the check is False, keep the loop going. 
                if isgood == False:
                    keep_going=True #Keep going, since this station's data was incomplete. 
                    
                    print(f"Station ID# {station['ID']}, called {station['Name']} is incomplete. Do not use it.")
        if isgood==False: #If the loop is over, and no stations are found, throw an error.
            print(f"I checked { len(self.closest_stations) } weather stations within {self.station_filters['Search Radius (deg)']} degrees of {point} in all directions.")
            print("None have data which is adequate for my use.")
            print("Please enter a different point, or change the years you must include in the load.")        
                  
    #This funtion runs all functions.
    def run_this_baby(self,point):
        if self.display: startTime = time.time() 
      
        self.station_data=self.load_station(point)
        if(self.display): 
            print(f"Time to download this station data: {time.time() - startTime} seconds.")
            startTime = time.time() 
        
        self.station_data=self.StationDataCleaner(self.station_data)
                
        if(self.display): 
            print(f"Time to clean up this data: {time.time() - startTime} seconds.")
            
    """    
    #This function as a default takes the LoadStation object's
    #reference point of interest and searches for the nearest point in our list of points,
    #It returns a list of the nearest stations, ordered by distance.
    #This is returned in the foramt of a dataframe.
    #The first row of the returned dataframe is the reference point. 

    """
    def nearest_station(self,point):
       pt=point
       #This initializes the actual data to search, from the list of stations
       #It only reads in the columns which are necessary to search by distance. 
       #It also searches only for the TMAX values. 
       df_all = pd.read_csv(self.STATION_META_FILENAME,
                        dtype={'Firstyear': np.int64,'Lastyear': np.int64})[['ID','Name','Latitude','Longitude','Firstyear','Lastyear']]
               
       #This steps strips out lat and lon values that are not nearby, reducing the number
       #of distance computations required.
       bar = self.station_filters['Search Radius (deg)']
       #Create a mask - only these rows are kept (i.e., those within the search radius. )
       mask = (df_all['Longitude']>point[1]-bar) & (df_all['Longitude']<point[1]+bar) & (df_all['Latitude']>point[0]-bar) & (df_all['Latitude']<point[0]+bar)
       #Create the filtered dataframe. 
       df = df_all[mask].copy(deep=True)

       #These lines strip out stations where there is no recent data (from within the last year), 

       recentyear=date.today().year-1
       df = df[df['Lastyear']>=recentyear] 
       #and then strips out stations for which data is only very very recent. 
       baseyear=self.station_filters["Earliest Year Required in Dataset"]
       df = df[df['Firstyear']<=baseyear]
       if len(df) == 0:
           print(f"There are no stations within {bar} degrees of {point} in any direction that have data as far back as the year {baseyear}.")
           print(f"Please enter a more realistic base year for later than {baseyear} or expand your search radius.")
           
           #This creates the frist row of the return dataframe.
           refdf= pd.DataFrame(
               {'ID':['REFPT'],
                'Name':["Reference Location"],
                'Latitude':np.nan,
                'Longitude':np.nan,
                'Miles_from_Ref':np.nan}
               )
           
           return refdf
       #Now drops the year info, since we don't need it.
       df=df[['ID','Name','Latitude','Longitude']]

        #Prints the number of stations being searched.
       
       print(f"Searching closest station among {len(df)} stations within {self.station_filters['Search Radius (deg)']} degrees of the reference.")
       print(f"Which have data for {recentyear} and at least as early as {self.station_filters['Earliest Year Required in Dataset']}.")
       
       #This creates a list, where every element is the distance from the Reference Point
       #to each of the weather stations which were found within search radius of Reference Point.
       d1er=[] #Initialize this list. 
       for index, station in df.iterrows():
           d1er.append(self.ts_latlon_distance([point[0],point[1]],station[['Latitude','Longitude']]))
       
       #This appends the data series to the final return dataframe
       df['Miles_from_Ref']=d1er
       
       #and then sorts by miles from ref, to find the closest stations       
       returner=df.sort_values(by='Miles_from_Ref')
       
       #Then strip out to include only the closest 50 stations.
       returner=returner[0:50]
       
       #This creates the frist row of the return dataframe.
       refdf= pd.DataFrame(
           {'ID':['REFPT'],
            'Name':["Reference Location"],
            'Latitude':[point[0]],
            'Longitude':[point[1]],
            'Miles_from_Ref':[0]}
           )
       final = refdf.append(returner)   
       return final
    
    #This function returns the miles between two points which were given in decimal lat/lon. 
    def ts_latlon_distance(self,latlon1,latlon2):
        # radians which converts from degrees to radians.
        lon1 = radians(latlon1[1])
        lon2 = radians(latlon2[1])
        lat1 = radians(latlon1[0])
        lat2 = radians(latlon2[0])  
        # Haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat / 2)**2 + cos(lat1) * cos(lat2) * sin(dlon / 2)**2 
        c = 2 * asin(sqrt(a))
        
        r = 3956    # Radius of earth in miles. Use 3956 for miles 
        # calculate the result, returned in miles
        return(c * r)
         
    #This loads the station.
    #As default, station = None, and the id_closest_st for thsi object is run.
    #The station value must be a string. 
    def load_station(self,station):
        #First, assigning the URL Name. This is standard and based on NOAA conventions. 
        csv_url = self.NOAA_URL+str(station)+'.csv'
        if self.display == True: 
            #This gets the filesize from NOAA.
            info=requests.head(csv_url)
            filesize = int(info.headers['Content-Length'])/1000000
            print(f"Beginning to download {csv_url} ({filesize} MB). Please wait while the file is transferred.")
        #Then, checking if it exists. the script breaks if the file does not exist.
        response = requests.get(csv_url)
        if response.status_code!= 200:
            print(f"Bad Station ID. The file called {csv_url} does not exist.")
            print(f"Please check your station ID {station} and re-submit.")
            sys.exit("Break Error in load_station of StationReader: Bad Station ID.")
        
        #Then, downloads the file from NOAA and returns it.          
        return pd.read_csv(csv_url)   

    """
      This takes in data for single weather station - in the format of output from
      self.load_station -- and cleans it up. 
      The only input is stationd, which must be in the format of self.station_data
     It returns a numpy array version, where each row contains:
         Year, Month, Day, Day of Year, TMAX, TMIN, TMID
    """
    def StationDataCleaner(self,stationd):
       
        #Does the swap of -9999 for np.nan (-9999 are no data values)
        returner = stationd.replace(-9999,np.nan) 
        
        #Calculate TMID. (Note that they provide a TAVG value, but even if it exists, it is often incomplete. So we always need TMID.)
        returner['TMID'] = (returner['TMAX']+returner['TMIN'])/2
        
        #Now, drop extraneous columns.
        returner = returner[['DATE','TMAX','TMIN','TMID']]
                    
        #Converting frmo tenths of degrees C to C.
        threetemps=["TMAX",'TMIN',"TMID"]
        returner[threetemps] = returner[threetemps]/10
        
        #Convert from Celsius to Farnehit
        returner[threetemps] = (returner[threetemps]*9/5)+32
        
        #Create dates columns appropriately.
        returner['DATE']=pd.DatetimeIndex(returner['DATE'])
        returner['Year']=returner["DATE"].dt.year
        returner['Month']=returner["DATE"].dt.month
        returner['Day']=returner["DATE"].dt.day
        returner['DOY'] = returner['DATE'].dt.day_of_year
        
        #Removes the DATE column, which is now redudant.
        returner = returner.drop(['DATE'],axis=1)
        
        #Re-order the columns for cleanliness.
        returner1=returner.reindex(columns=['Year','Month','Day','DOY','TMAX','TMIN',"TMID"])
        #Change to numpy.
        returner1 = np.asarray(returner1)
        
        #THESE STEPS SCRUB OUT ALL YEARS WHERE THERE IS INSUFFICIENT DATA.
        #First, list out the unique years.
        uniqueyears = np.unique(returner1[:,0])
            
        showme = self.display

        #This initializes the "to keep" list -- the list of data to be included.
        tokeep = np.empty([0,1],dtype=int)
        #Loop over all years.
        for year in uniqueyears:
            #As a default, the year is not returned. This flag therefore starts as FALSE.
            return_this_year = False    
            #Select only the index values where the year is equal to year.
            index=np.where(returner1[:,0]==year)
            #Then, selects data for this year
            #All columns are included.
            yearsubset=returner1[index]
            
            #Finding all unique months in this year.
            uniquemonths = np.unique(yearsubset[:,1])
            #Coutns the number of months.
            norows = len(uniquemonths)
            #Print an update, if the year is removed.
            if (showme) and (norows < 12) and (year != uniqueyears[-1]):
                print(f"I dropped {year} for having only {norows} months of data.")
            #Checks that 12 months are actually present in the data.
            if norows == 12:
                #This becomes True for now, but becomes false if there are insufficient days
                #in any one month.
                return_this_year=True
                #Then, checks there are adequate data in each month for this year.
                for month in uniquemonths:
                    #Creates a list just including these months. 
                    monthssubset = yearsubset[np.where(yearsubset[:,1]==month)]
                    #Finds all numbered values.
                    listoftmax = ~np.isnan(monthssubset[:,4])
                    listoftmin = ~np.isnan(monthssubset[:,5])
                    listoftmid = ~np.isnan(monthssubset[:,6])
                    #Counts the number of numbers. 
                    tmax_number = np.count_nonzero(listoftmax)
                    tmin_number = np.count_nonzero(listoftmax)
                    tmid_number = np.count_nonzero(listoftmax)
                    #Then checks if each is over the require dminium.
                    #Every data point must be present at the right threshold.
                    tmax_enough = (tmax_number >= self.station_filters['Min Days in Every Month'])
                    tmin_enough = (tmin_number >= self.station_filters['Min Days in Every Month'])
                    tmid_enough = (tmid_number >= self.station_filters['Min Days in Every Month'])
                    
                    #If any condition is false -- there are not enough TMID, TMAX, or TMIN
                    #values in the month - then the whole year is dropped out.
                    if (tmax_enough == False) or (tmin_enough==False) or (tmid_enough ==False):
                        return_this_year = False
                        if showme: print(f"I dropped {year} for having less than {self.station_filters['Min Days in Every Month']}s days in month # {month}.")
            #Only if return_this_year, the flag field, is True, is the year added.
            #there is only one exception: the most recent year, which definitionally
            #will not have a complete reord of data, which is OK.
            if (return_this_year ==True) or (year == uniqueyears[-1]):
                #If everything is true; then add it to the list of rows to keep. 
                tokeep=np.append(tokeep,index)
        #Filter to just the years meeting all the conditions above.
        returner1 =returner1[tokeep]          

        return returner1
    
    #This function reviews the station data and checks that it is complete,.
    #using a few common-sense checks.    
    def StationDataCheck(self,stationd):
        if self.display: 
            print("-------------------------------------")
            print("Checking the completeness of the data.")
        #List unique years in the dataset.
        listyears = np.sort(np.unique(stationd[:,0]))
        #This is a flag field which defaults to TRue, but remains true
        #only if the dataset passes all checks.
        itsgood=True
        
        thisyear = date.today().year
        rec_years =int(self.station_filters['Min Number of Years before Present Year'] )  #Create a shortened variable name.
        #Ensure the years are sorted, because they do not work otherwise. 
        recent_years_in_data=pd.Series(listyears[-rec_years-1:-1]).astype(int) #List the recent years in the dataset.
        recentyearslist=pd.Series(np.arange(thisyear-rec_years,thisyear)).astype(int) #The list of recent years that should be present. 
        if len(recent_years_in_data) != len(recentyearslist): 
            itsgood = False 
            print(f"Flag: You specified that I need {rec_years} years of data before {thisyear}, but there are only {len(recent_years_in_data)} years.")
        
        #Then, check if there are sufficient years in the recent trend.
        min_trend_years = self.station_filters['Min Years in Last 30 Years for Trend']
        no_recent = np.shape(np.where(listyears>=thisyear-min_trend_years))[1] #Count the number of recent years. 
        if no_recent <= self.station_filters['Min Years in Last 30 Years for Trend']:
            itsgood = False
            
            print("Flag: Insufficent years to calculate recent trend.") 
            print(f"You specified that {self.station_filters['Min Years in Last 30 Years for Trend']} years of data in the last 30 be present to calculate a trend, but only {no_recent} are present.")
        #Then, that there are enough years to calculate a baseline.
        #First, limit it to the years in the baseline range. 
        early_index = np.where(listyears<=self.station_filters['Last Year of Baseline'])
        early_index = np.where(listyears[early_index]>=self.station_filters["Earliest Year Required in Dataset"])
        #Then, count the years.
        no_early = np.shape(early_index)[1]
        
        if no_early <= self.station_filters['Required Number of years in Baseline Period']:
            itsgood = False
            
            print(f"Flag: Insufficent years before {self.station_filters['Last Year of Baseline']} to calculate mean for baseline period.")
            print(f"There are {no_early} years in the period {self.station_filters['Earliest Year Required in Dataset']} to {self.station_filters['Last Year of Baseline']}, but I need {self.station_filters['Required Number of years in Baseline Period']} to set an accurate baseline.")
                
        return itsgood
              





