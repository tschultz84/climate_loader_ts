# -*- coding: utf-8 -*-

#This readme file describes the parameters in load_stats_static.yaml

#NOAA_URL is the url to the directory containing the DLY files.
    NOAA_URL : 'https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/'

#==============
    #BASEYEAR is the year before which station data must be availbale.
    #stations with no data prior to this year are not loaded.
    BASEYEAR : 1970

    #FIRST_BASE_YEAR is the first year for which data must be available.
    #stations without data going back to this year are not loaded.
    FIRST_BASE_YEAR : 1920

    #LAST_BASE_YEAR is the last year in the baseline for which data must be available.
    #stations without data going back to this year are not loaded.
    LAST_BASE_YEAR : 1955

    #BASENOYEARS is the number of years included in the baseline period.
    BASENOYEARS : 30

    #RECENT_TREND_YEARS is the number of years included (starting from the final year in the dataset).
    #in the annual temperature average trendline.  The linear trend is calculated only over this period.
    RECENT_TREND_YEARS : 30

    #REQUIRED_TREND_YEARS is minimum number of years required to calculate this trend.
    #Weather stations without this many years present will be discarded.
    REQUIRED_TREND_YEARS : 25
    
    #MIN_RECENT_YEARS is the minimum number of years before he present which must be present
    #for the statin to be loaded.
    #if this is 5, and the year is 2022, then every year since 2017 must be present for the weather station to load.
    #set this to zero to turn this off. 
    MIN_RECENT_YEARS : 0

    #MIN_DAYS_PER_MO is the minimum number of dasys of data which must be present
    #for every month of a year for a dataset to be included.
    # If not, then the year is excluded from dataset.
    MIN_DAYS_PER_MO : 20
    
    #ALPHA is the alpha value cmoparing TMID in the reference and baseline period. 
    #Basically if the P Value from the corresponding t-test exceeds this value,
    #then a change is deteced.
    ALPHA : 0.01
    
    
    #This is the lat/lon distance within which weather stations are serached.
    #This is used to streamline how fast the "closest station" function works. 
            #stations are only sorted by distance if they rae within +/- coord_radius
            #of the reference point. 
    SEARCH_RADIUS : 3


#=============
#REQUIRED SPECIFICS TO LEAD CSV FILES FROM NOAA WEBSITE
    #Columns to keep from CSV files.
    KEEP_COLS : ['DATE','TMAX','TMIN','TMID']


...