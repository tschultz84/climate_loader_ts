# -*- coding: utf-8 -*-


import pandas as pd
from math import radians, cos, sin, asin, sqrt
import time
import sys
#import yaml
import numpy as np
from datetime import date
import requests

"""
This object has functions to identify the station nearest to a reference lat, lon.
It then downloads the data from the climate station closest to this lat lon.
It then runs several checks on this data to see if it is acceptably complete, which you can infer from the input variables. 
It also calculates the TMID values, which is the average of TMAX and TMIN in every day.

#The outputs of interest are as follows:
    #self.refpoint is the reference
    #self.id_closest_station gives the ID of the nearest station
    #self.name_closest_statioon gives its name
    #self.miles_from_ref gives its distance from teh reference point
    #self.closest_stations gives a list of the  nearest stations
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
    *firstyear,lastyear are the first and last years which must be present in the dataset record for the weather station
        to be loaded. 
    *basenoyears is the number of years in the baseline period.
    *min_recent_years is the minimum number of years before he present which must be present
        for the statin to be loaded.
    *required_trend_years minimum number of years required to calculate this trend.
        Weather stations without this many years present will be discarded.
    *lastbaseyear is the last year in which the baseline is calculated (not the same as last year)
    
"""    
class LoadStation :
    def __init__(self,point,printudpate=False,min_days_per_mo=15,search_radius=1,firstyear=1890,lastyear=2020,basenoyears=30,
            min_recent_years=5,required_trend_years =20,lastbaseyear=1955):
        
        #Initiate all the variables which are used later on, and load them into a Series for easy review.
        self.station_filters=pd.Series(data={
            "Reference Point Latitude":point[0],
            "Reference Point Longitude":point[1],
            "Earliest Year Required in Dataset":firstyear,
            "Latest Year Required in Dataset":lastyear,
            "Min Days in Every Month":min_days_per_mo,
            "Search Radius (deg)":search_radius,
            "Min Number of Years before Present Year":min_recent_years,
            "Min Years for Trend":required_trend_years,
            "Last Year of Baseline":lastbaseyear,
            "Required Number of years in Baseline Period":basenoyears    })
        
        self.display=  printudpate   
        #If a point (lat, lon) is passed, then create a point. 
        if (type(point) != str) and (len(point)==2):
            self.refpoint=point
        
        #Files to load.
        self.GHCND_STATIONS_FILENAME = "data\\ghcnd-stations.txt" #The NOAA file containing station meta info (lat, lon, elevation.)
        self.STATION_MASTER_ALL_FILENAME = "data\\ghcnd_station_master_ts_all.csv" #contains all stations and all data about each.
        self.STATION_META_FILENAME = "data\\ghcnd_station_master_ts_tmax.csv" #containing all the lat / lon information for the stations.
        self.NOAA_URL = 'https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/' #url to the directory containing the DLY files. 

        #First, the list of stations is generated.
        self.closest_stations = self.nearest_station(point)
      
        #Now, the station data is all loaded up and then all stations are run through to find one which is complete.
        keep_going = True# add a flag variable
        for index,station in self.closest_stations.iloc[1:].iterrows():
            if keep_going == True:
                self.run_this_baby(station['ID']) #Run the download.
            
                #A check is performed tos ee if the data is compelte.
                isgood = self.StationDataCheck(self.station_data)
                if isgood == True: #Set flag to stop the loop, if the data is good.
                    keep_going=False
                    #SEtting variables to be extracted. 
                    self.name_closest_station=station['Name']
                    self.station_information = pd.Series(data={'Weather Station Name':station['Name']})
                    self.id_closest_station=station['ID']
                    self.station_information = self.station_information.append(pd.Series(data={'Station ID':station['ID']}))
                    #Lat Longitude.
                    self.st_latlon_str = f"{round(station['Latitude'],2)} latitude {round(station['Longitude'],2)} longitude"
                    
                    #Latitude, Longitude of the weather station.
                    self.st_latlon=[station['Latitude'],station['Longitude']]
                    self.station_information = self.station_information.append(pd.Series(data={'Station Latitude':station['Latitude']}))
                    self.station_information = self.station_information.append(pd.Series(data={'Station Longitude':station['Longitude']}))
                    
                    #Miles from reference point.
                    self.miles_from_ref = station['Miles_from_Ref']
                    #Add in the original reference point for reference.
                    self.station_information = self.station_information.append(
                        pd.Series(data={
                            "Reference Point Latitude":point[0],"Reference Point Longitude":point[1]}))

                   
                    self.station_information = self.station_information.append(pd.Series(data={'Miles from Reference Point':station['Miles_from_Ref']}))
                    if self.display:  
                        
                        print(f"Station ID# {station['ID']}, called {station['Name']} is complete. It's good to use.")
                        print("This station is "+str(self.miles_from_ref)+" miles from the reference point.")
                        
                #If it's not, we report that. ANd keep the flag TRue to keep goign. 
                if isgood == False:
                    keep_going=True
                    if self.display: 
                        print(f"Station ID# {station['ID']}, called {station['Name']} is incomplete. Do not use it.")
        if isgood==False: #If the loop is over, and no stations are found, throw an error.
            print(f"I checked { len(self.closest_stations) } weather stations within {self.station_filters['Search Radius (deg)']} degrees of {point} in all directions.")
            print("None have data which is adequate for my use.")
            print("Please enter a different point, or change the years you must include in the load.")
            #Set variables to return, which all are not numbers, since no reference station was found.
            self.id_closest_station = np.nan
            self.name_closest_station = np.nan
            self.miles_from_ref = np.nan
            self.st_latlon = np.nan
                    
                    
                  
    #This funtion runs all functions.
    def run_this_baby(self,point):
        if self.display: startTime = time.time() 
      
        self.station_data=self.load_station(point)
        if(self.display): 
            print("Time to load this station data:" + str(time.time() - startTime))
            startTime = time.time() 
        
        self.station_data=self.StationDataCleaner(self.station_data)
                
        if(self.display): 
            print("Time to clean up this data:" + str(time.time() - startTime))
            
    """    
    #This function as a default takes the LoadStation object's
    #reference point of interest and searches for the nearest point in our list of points,
    #searching distance by miles distance
    #It returns a list of the nearest stations, ordered by distance.
    #This is returned in the foramt of a dataframe.
    #The first row of the returned dataframe is the reference point. 

    """
    def nearest_station(self,point):
       pt=point
       #This initializes the actual data to search, from the list of stations
       #It only reads in the columns which are necessary to search by distance. 
       #It also searches only for the TMAX values. 
       df1 = pd.read_csv(self.STATION_META_FILENAME,
                        dtype={'Firstyear': np.int64,'Lastyear': np.int64})[['ID','Name','Latitude','Longitude','Firstyear','Lastyear']]
               
       #This steps strips out lat and lon values that are not nearby, reducing the number
       #of distance computations required.
       #bar = self.search_radius
       bar = self.station_filters['Search Radius (deg)']
       
       df=df1[df1['Longitude']>point[1]-bar]
       df=df[df['Longitude']<point[1]+bar]
       df=df[df['Latitude']>point[0]-bar]
       df=df[df['Latitude']<point[0]+bar]
       
       #These lines strip out stations where there is no recent data (from within the last year), 

       #recentyear=date.today().year-1
       recentyear=self.station_filters['Latest Year Required in Dataset']
       df = df[df['Lastyear']>=recentyear] 
       #and then strips out stations for which data is only very very recent. 
       baseyear=self.station_filters["Earliest Year Required in Dataset"]
       df = df[df['Firstyear']<=baseyear]
       if len(df) == 0:
           print(f"There are no stations within {bar} degrees of {point} in any direction that have data as far back as the year {baseyear}.")
           print(f"Please enter a more realistic base year for later than {baseyear} or expand your search radius.")
           #sys.exit("Base year is too early.")
           #Set variables to return, which all are not numbers, since no reference station was found.
           self.id_closest_station = np.nan
           self.name_closest_station = np.nan
           self.miles_from_ref = np.nan
           self.st_latlon = np.nan
           
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
       if self.display: 
           print("Searching closest station among "
                              +str(len(df))+" stations within "+str(self.station_filters["Search Radius (deg)"])+" degrees of the reference.")
           print("which have more recent data than "+str(recentyear)+" and at least as early as "+str(self.station_filters["Earliest Year Required in Dataset"]))
       
       #This creates a dataseries which calculates the distance from the ref "pt"
       #to all of the nearest stations.
       d1er=[]
       for i in range(0,len(df)):
           d1er.append(self.ts_latlon_distance([point[0],point[1]],df.iloc[i][['Latitude','Longitude']]))
       
       #This appends the data series to the final return dataframe
       df['Miles_from_Ref']=d1er
       
       #and then sorts by miles from ref, to find the closest stations       
       returner=df.sort_values(by='Miles_from_Ref')
       
       #Then return the closest 50 stations.
       returner=returner[0:50]
       
       #This sets a variable which is the ID Of the nearest station
       self.id_closest_station = returner.iloc[0,0]
       #Then a separate variable that is the name. 
       self.name_closest_station = returner.iloc[0,1]
       #Then a separate variable that is the miles distance from the ref point. 
       self.miles_from_ref = returner.iloc[0,4]
       
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
    
    #This function returns the miles between two poitns
    #def ts_latlon_distance(lat1, lat2, lon1, lon2):
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
        
        # Radius of earth in miles. Use 3956 for miles
        r = 3956     
        
        # calculate the result, returned in miles
        return(c * r)
         
    #This loads the station.
    #As default, station = None, and the id_closest_st for thsi object is run.
    #The station value must be a string. 
    def load_station(self,station=None):

        #If statino is none, then defaults to the station closest to the ref point. 
        if station == None:
            station=self.id_closest_station
       
        #Downloading the CSV file straight rom NOAA.
        #First, assigning the URL Name.
        csv_url = self.NOAA_URL+str(station)+'.csv'
        
        if self.display == True: 
            #This gets the filesize.
            info=requests.head(csv_url)
            filesize = int(info.headers['Content-Length'])/1000000
            print("Beginning to download "+csv_url+' ('+str(filesize)+"MB). Please wait while the file is transferred.")
        #Then, checking i it exists. the script breaks if the file does not exist.
        response = requests.get(csv_url)
        if response.status_code!= 200:
            print("Bad Station ID. The file called "+csv_url+" does not exist.")
            print("Please check your station ID "+str(station)+" and re-submit.")
            sys.exit("Break Error in load_station of StationReader: Bad Station ID.")
        
        #Then, downloads the file and returns it. 
        inter0 = pd.read_csv(csv_url)            
        return inter0
 
    #This function finds the Day of year, assuming every month has 31 days
    #taking in the day of month (dom) and month.
    def find_doy(self,dom,month):
        doy = dom+(month-1)*31
        return doy
    """
      This takes in data for single weather station - in the format of output from
      load_station and cleans it up. 
      The only input is stationd, which must be in the format of self.station_data
     It returns a numpy array version, where each row contains:
         Year, Month, Day, Day of Year, TMAX, TMIN, TMID
    """
    def StationDataCleaner(self,stationd):
       
        #Does the swap of -9999 for np.nan (I am not sure this is relevant for CSV)
        returner = stationd.replace(-9999,np.nan) 
        
        #Calculate TMID.
        #Note that you cna't just load TAVG.
        #EVen if it exists, it is often incomplete.
        returner['TMID'] = (returner['TMAX']+returner['TMIN'])/2
        
        #Now, drop extraneous columns.
        returner = returner[['DATE','TMAX','TMIN','TMID']]
                    
        #Converting frmo tenths of degrees C to C.
        threetemps=["TMAX",'TMIN',"TMID"]
        returner[threetemps] = returner[threetemps]/10
        
        #Convert from Celsius to Farnehit
        returner[threetemps] = (returner[threetemps]*9/5)+32
        
        #Create dates columns appropriately.
        returner['Year']=returner['DATE'].str[:4]
        returner['Month']=returner['DATE'].str[5:7]
        returner['Day']=returner['DATE'].str[8:10]
        
        returner['Year']=pd.to_numeric(returner.Year)
        returner['Month']=pd.to_numeric(returner.Month)
        returner['Day']=pd.to_numeric(returner.Day)
        #Creates a column for the DAy of the Year.
        returner['DOY'] = self.find_doy(returner['Day'],returner['Month'])
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
        #Loop oer all years.
        for year in uniqueyears:
            #As a default, the year is not returned. This flag therefore starts as FALSE.
            return_this_year = False    
            #Select only the index values wehre the year is equal to year.
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
                print("I dropped "+str(year)+" for having only "+str(norows)+" months of data.")
            #Checks that 12 months are actually present in the data.
            if norows == 12:
                #This becomes True for now ,but becomes false if there are insufficeitn days
                #in any one month.
                return_this_year=True
                #Then, checks there are adequate data in each month for this year.
                for month in uniquemonths:
                    #Creates a list jsut including tehse months. 
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
        #Pull the very first year in the dataset.
        self.veryfirstyear=np.unique(returner1[:,0])[0]
        return returner1
    
    #This function reviews the station data and checks that it is complete,.
    #using a few common-sense checks.    
    def StationDataCheck(self,stationd):
        if self.display: 
            print("-------------------------------------")
            print("Checking the completeness of the data.")
        #List unique years in the dataset.
        listyears = np.unique(stationd[:,0])
        #This is a flag field which defaults to TRue, but remains true
        #only if the dataset pasess al checks.
        itsgood=True
        
        """
        thisyear = date.today().year
        #Check the latest few years are present.
        if (listyears[-1]!=thisyear) or (listyears[-2]!=thisyear-1) or (listyears[-3]!=thisyear-2):
            itsgood = False
            if self.display: print("Flag: I need every year of ("+str(thisyear-3)+" to "+str(thisyear)+") - but they are not all not present.")
        """
        #Ensure the years are sorted, because they do not work otherwise. 
        listyears=np.sort(listyears)
        
        thisyear = date.today().year

        #If MIN_RECENT_YEARS = 0, then the check is skipped entirely.
        if self.station_filters['Min Number of Years before Present Year'] >0:
            #This is the value that should be recent in this part of the array, if all years are present. 
            threshyear=thisyear-self.station_filters['Min Number of Years before Present Year']
            
            #But this finds out if the threshyear is actually the right value in the array
            #This was throwing errors I think because the number of yeras was incorrect, need to return to this line.
            #arbiter = (listyears[-self.station_filters['Min Number of Years before Present Year']-1] - threshyear >= 0)
            
            #Check the latest few years are present.
            #If not, then itsgood is false and the year is bad. 
            #if arbiter == False:
            #    itsgood = False
            #    print("Flag: I need every year of ("+str(threshyear)+" to "+str(thisyear)+") - but they are not all not present.")

        
        #Then, check if there are sufficient years in the recent trend.
        #This is just the number of trend years plus 5.
        min_trend_years = self.station_filters['Min Years for Trend']
        no_recent = np.shape(np.where(listyears>=thisyear-min_trend_years))[1]
        if no_recent <= self.station_filters['Min Years for Trend']:
            itsgood = False
            if self.display: 
                print("Flag: Insufficent years to calculate recent trend.") 
                print("I need more than "+str(self.station_filters['Min Years for Trend'])+" years to calculate a trend, but only "+str(no_recent)+" available.")
        #Then, that there are enough years to calculate a baseline.
        #First, limit it to the years in the baseline range. 
        early_index = np.where(listyears<=self.station_filters['Last Year of Baseline'])
        early_index = np.where(listyears[early_index]>=self.station_filters["Earliest Year Required in Dataset"])
        #Then, count the years.
        no_early = np.shape(early_index)[1]
        
        if no_early <= self.station_filters['Required Number of years in Baseline Period']:
            itsgood = False
            if self.display: 
                print("Flag: Insufficent years before "+str(self.station_filters['Last Year of Baseline'])+" to calculate mean for baseline period.")
                print(f"There are {no_early} years in the period {self.station_filters['Earliest Year Required in Dataset']} to {self.station_filters['Last Year of Baseline']}, but I need {self.station_filters['Required Number of years in Baseline Period']} to set an accurate baseline.")
                
        return itsgood
              

bzdata=LoadStation([45.647256643331126,-111.04060494981753],True,15,1,1920,2020) #Loading in Bozeman, MT coordinates
print(bzdata.station_information)



