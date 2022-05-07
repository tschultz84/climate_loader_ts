# README
This file describes key environment variables, how to set up the package, and the needed program
inputs and outputs.

# OVERVIEW
These Python modules allow you to search for the land-based weather station nearest to you
and review historical temperature data at that site. You provide the programs with a latitude and longitude,
and a time period of interest; it automatically finds the closest weather station and provides an analysis for you. 

# About GHCND
This entire analysis is based upon the NOAA Dataset called Global Historical Climatology Network daily (GHCNd). To quote NOAA:
"The Global Historical Climatology Network daily (GHCNd) is an integrated database of daily climate summaries from land surface stations across the globe. 
GHCNd is made up of daily climate records from numerous sources that have been integrated and subjected to a common suite of quality assurance reviews.
GHCNd contains records from more than 100,000 stations in 180 countries and territories."

NOAA provides GHCNd data through the National Centers for Environmental Information (NCEI). See 
https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily

# ENVIRONMENT
Make sure the following files are all in your root directory:
* Station_Loader_ts.py and Station_analyzer_ts.py contain the Python modules (see their descriptions below)
* filesnames.yaml and load_stats_static.yaml contain, respectively, filenames to underlying data files and parameters used by StationAnalyzer. 
Files in the "/data" directory are vital to the correct operation of programs, especially ghcnd_station_mast_ts_tmax.csv. StationLoader will not work if
data/ghcnd_station_mast_ts_tmax.csv does not exist.
* examples.ipynb shows you how to run the Python modules. You don't technically need it, but I would recommend taking a look and then making a copy of this
to run other weather station analyses.
* create_station_master_1.pynb was used only to generate ghcnd_station_mast_ts_tmax.csv from underlying metadata in the GHCND database. You can ignore it. 


You will also need to ensure your Python programs know where to look, either by setting the current working directory using os.cwd() to the root directory or adding a new .pth file in your site-packages
directory that points there. 


# StationLoader and StationAnalyzer Classes
There are two modules in two Python files, each containing a Class:

StationLoader class (in Station_Loader_ts.py)
StationLoader takes in two input fields: 
* "point" in the form [float,float] which is the latitude and longitude of interest to you, and 
* "printupdate", which is TRUE or FALSE, and just indicates whether processing updates are returned to your Python console. 
* There are several parameters in load_stats_static.yaml of importance: 
    #BASEYEAR is the year before which station data must be availbale. Weather stations with no data prior to this year are not loaded.  (default is 1970)
    #BASENOYEARS is the number of years included in the baseline period (i.e., the baseline period includes BASENOYEARS after the earliest available year of data). 
    #MIN_DAYS_PER_MO is the minimum number of days of data which must be present for every month of a year for a dataset for the year to be included.  If not, then the entire year is excluded from dataset.
    #SEARCH_RADIUS this is the lat/lon distance from point within which weather stations are serached. This is used to streamline how fast the "closest station" function works. Stations are only sorted by distance if they rae within +/- coord_radius of point. 
    

StationAnalyzer class (in Station_analyzer_ts.py)
StationAnalyzer takes in four input fields:
* stationdata is a numpy array that matches the format of StationLoader.station_data (see below). (StationAnalyzer is designed
to receive inputs directly from StationLoader.)
*refst (e.g. '2020-01-31') is a string, in format YYYY-MM-DD, that is the beginning of the reference period you want analyzed.
*refend (e.g., '2020-12-31')  is a string, in format YYYY-MM-DD, that is the end of the reference period you want analyzed.
*display is TRUE or FALSE. If True, it automatically displays all the outputs of StationAnalyzer (tables and charts). If you want to scan
results for a single weather station, have display=TRUE. If you want to analyze a large number of weather stations, set it to False. 

#  Attributes of StationLoader 
StationLoader is a class. After you run StationLoader, the following object attributes are available.
* StationLoader.closest_stations is a list of the weather stations closest to the point lat ,lon, values, along with their index numbers, names, locations, and distances from point.
* StationLoader.name_closest_station and self.id_closest_station are strings that are respectively the name and ID of the closest weather station meeting data requirements (see below) to point. 
* StationLoader.st_latlon contains the latitude and longitude of the weather station as al ist: [lat,lon]
* StationLoader.miles_from_ref is the distance, in miles, of the closest weather station meeting data requirements
* StationLoader.station_data is a Nx7 numpy array that is the key returned attribute of interest. 

# Notes on StationLoader.station_data and Data Requirements for Selected Weather Station
StationLoader.station_data contains the data for name_closet_station. But this station is only chosen if it satisfies some criteria, evaluated as follows.
* After loading, all years that have less than 12 months of data are deleted (this would heavily bias data in these years).
* Additionally, every single month in a year must have more than yaml.MIN_DAYS_PER_MO days. If it has less than that, then the entire year's dataset is removed.
* CRITERIA 1: For a weather station to be used, it must have at least yaml.BASENOYEARS worth of data prior to the baseline year yaml.BASEYEAR.  
* CRITERIA 2: The weather station must have the latest 3 years of data complete and present (e.g., if it's 2022, complete years of data must be available for 2021, 2020, and 2019.)
* CRITERIA 3: There must be enough data in the past few years to calculate a trend. Specifically, the last 30 years must have at least yaml.['REQUIRED_TREND_YEARS'] worth of good data.
If any of these 3 criteria are not met, the weather station is not used, and the next closest weather station is instead loaded. 

# Structure of StationLoader.station_data
Each row of StationLoader.station_data corresponds to a single day, and contains TMAX, TMIN, and TMID, values. The columns are as follows:
[YEAR,MONTH, DAY OF MONTH, DAY OF YEAR, TMAX, TMIN, TMID].
The Year, Month, and Day of Month, values, are just as you expect.  
Day of Year is a derived value, and it is not the true "Day of Year", counting upwards from Day 1 (January 1). Rather, Day of Year is equal to (Day of Month)+(Month-1)*31. 
It essentially conts every month as having 31 days, for computational simplicity. 

TMAX, TMIN, TMID
StationLoader.station_data loads directly the TMAX and TMIN values for every day. These are the actually measured temperature data values in the GHCND. "TAVG" values are actually only derived
products that are equal to (TMAX + TMIN )/2 for every day. But TAVG values are provided only sporadically. So StationLoader also calculates the "TMID" value (equivalent to the normal
T average value), based upon TMAX and TMIN values. Therefore every single entry in self.station_data will have a TMAX, TMIN, and TMID value, and TMID will be equal to (TMAX+TMIN)/2
for every row. 
It also calculates the TMID values, which is the average of TMAX and TMIN in every day.

#  Attributes of StationAnalyzer 

StationAnalyzer.all_years_mean is an Nx2 numpy array. THe firts column is the year. The second is the average of TMID in that year for the days in the chosen reference period.
StationAnalyzer.tmid_ref_data and StationAnalyzer.tmid_base_data are the daily data values for just the reference and baseline range.
StationAnalyzer.key_metrics_table is a pandas dataframe containing all the climate metrics of interest across All Time and in the reference and baseline periods.
* Average, Variance, Maximum, and Minimum Temperatures, over the All-Time Record, Baseline, and Reference periods
* Dates of the Maxium and Minimum Temperatures for All-Time, Baseline, and Reference Periods
StationAnalyzer.key_stats is a pandas dataframe containing the statistics of significance about climate change. 
* The difference between the average temperatures in the Baseline and Reference Periods
* The recent trend in temperatures over a time period set by the YAML file variable RECENT_TREND_YEARS
(This trend is calculated using LinearRegression() from sklearn.linear_model)
* The Statistical Significance of these differences, and accompanying Alpha Values

# Calculation of Statistical Significance
Two statistical calculations are completed:
* Comparing the average temperatures in the Baseline and Reference Period. The statistical signifiance
is calculated using a two-sided t-test (using Python function stats.ttest_ind). 
The difference is also visualized in a histogram showing the
temperatures in the baseline and reference periods.
* The trend in temperatures in the reference period (days) over the most recent RECENT_TREND_YEARS (in the YAML file).
The significance is calculated using a Pearson's correlation using Pearon's r from scipy.stats.stats.
This is visualized on the chart of annual change, as a trend line.

For all calculations of statistical signifiance, ALpha is controlled using ALPHA in the YAML file.

#charts
Three charts are shown
* Frequency histogram of TMID values in the Reference and Baseline periods.
* ALl-time Monthly average temperatures.
* Averages across the Reference period days in every year in the temperature record. ON this chart
are superimposed the baseline and reference years, and a recent trend over YAML RECENT_TREND_YEARS