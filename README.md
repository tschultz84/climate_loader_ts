# Package Author and Use
This package was originally written by Tobias Schultz in 2022. First relase on May 10, 2022.

Feel free to use as you see fit, but do attribute it to me if you ever publish or write anything. Thanks!

# README
This file describes key environment variables, and the needed program
inputs and outputs.

# OVERVIEW
This Python module allows you to search for the land-based weather station nearest to you. You provide the programs with a latitude and longitude,
and a time period of interest; it automatically finds the closest weather station. In addition, it filters the weather stations to find ones meeting
certain quality criteria, which you specify as inputs to the class. The output is a numpy array containing daily temperature data, meant to be analyzed by other packages. 

# About GHCND
This data is all sourced from the NOAA Dataset called Global Historical Climatology Network daily (GHCNd). To quote NOAA:
"The Global Historical Climatology Network daily (GHCNd) is an integrated database of daily climate summaries from land surface stations across the globe. 
GHCNd is made up of daily climate records from numerous sources that have been integrated and subjected to a common suite of quality assurance reviews.
GHCNd contains records from more than 100,000 stations in 180 countries and territories."

NOAA provides GHCNd data through the National Centers for Environmental Information (NCEI). See 
https://www.ncei.noaa.gov/products/land-based-station/global-historical-climatology-network-daily

# ENVIRONMENT
It should work with Python 3.9.7 and above. 

Make sure the following files are all in your root directory:
* **Station_Loader_ts.py** containa the Python module
* Files in **/data** directory are vital to the correct operation of programs, especially **ghcnd_station_mast_ts_tmax.csv**. StationLoader will not work if
**data/ghcnd_station_mast_ts_tmax.csv** does not exist.
* **examples.ipynb** shows you how to run the Python modules. You don't technically need it, but I would recommend taking a look and then making a copy of this
to run other weather station analyses.
* **create_station_master_1.pynb** was used only to generate ghcnd_station_mast_ts_tmax.csv from underlying metadata in the GHCND database. You can ignore it. 

You will also need to ensure your Python programs know where to look, either by setting the current working directory using os.cwd() to the root directory or adding a new .pth file in your site-packages
directory that points there. 

You also need to ensure these packages are available in your environment: *pandas, math, time, sys, numpy, datetime, requests, warnings*

# StationLoader Class
There is one module in one Python (.py) file, which contains one Class, Station_Loader_ts (in Station_Loader_ts.py)

StationLoader takes in these input fields: 
1. **point** in the form [float,float] which is the latitude and longitude of interest to you. 
2. **printupdate**, which is TRUE or FALSE, and just indicates whether detailed processing updates are returned to your Python console. 
3. **min_days_per_mo** is the minimum number of days of data which must be present for every month of a year for a dataset for the year to be included.  If not, then the entire year is excluded from dataset.
4. **search_radius** is the lat/lon distance from **point** which is searched for weather stations. This is used to streamline how fast the "closest station" function works. Stations are only reviewed for quality if they are within +/- **search_radius** of **point**. 
5. **firstyear** is the earliest year which must be present in the dataset record for the weather station
        to be loaded. If the weather station does not have some data present going back to this year, then it is rejected and another weather station will be reviewed instead.
6. **lastbaseyear** is the last year in which the baseline will be calculated. The baseline period is defined as the period between **firstyear** and **lastbaseyear**. This variable is important for now because there must be **basenoyears** of data available in the baseline period. 
7. **basenoyears** is the number of years required in the baseline period. For example, if **basenoyears == 30**, then there must be 30 complete years of data between **firstyear** and **lastbaseyear**.
8. **min_recent_years** is the minimum number of years before the present which must be present for the statin to be loaded.
9. **required_trend_years** is the minimum number of years in the last 30 required to calculate a trend.   Weather stations without this many years present in the last 30 years of data will be discarded.
    
    
#  Attributes of StationLoader 
StationLoader is a class. After you run StationLoader, the following object attributes are available.
1) **StationLoader.station_data** is a Nx7 numpy array that is the actual temperature data for the selected weather station. Each row is a day; and it has these 7 columns: *[Year, Month, Day, Day of Year, TMAX (max. temperature), TMIN (min. temperature), TMID (midpoint temperature between TMAX and TMIN)]*.
2) **StationLoader.station_information** is a pandas DataFrame containing meta information for the selected weather station, like the station name, ID, lat/lon, etc.  
3) **StationLoader.closest_stations** is a pandas DataFrame which lists the weather stations closest to **point** (closest based on number of miles between) values, along with their index numbers, names, locations, and distances from point.
4) **StationLoander_station_filters** is a pandas DataFrame which summarizes all the filters you entered as inputs (e.g., firstyear, search_radius, etc).

# Notes on StationLoader.station_data and Data Requirements for Selected Weather Station
**StationLoader.station_data** contains the data for a single weather station. It found this weather station via this process:
1) Listing all weather stations within **search_radius** degrees, then sorting them by miles distance. *(Note: This list of stations sorted by distance is **StationLoader.closest_stations**.)*
2) The data for the closest station in the list is downloaded from the NOAA website (CSV format).
3) The weather station is read using pandas, then cleaned, adding date fields and calculating TMID for each day. 
4) One key cleaning step is then completed: Every year must have at least 12 months of data in its record. All years that have less than 12 months of data are deleted, as this would heavily bias data in these years towards the seasons where data is present. *(For example: If you only had data on July on August, and not December and January, that would  make the year seem extra hot.)*
5) With the cleaning and step #4 filter completed, now the dataset is checked against every one of your input "filter" parameters. Every single filter must be met, or the station is disqualified - and the next closest station in **StationLoader.closest_stations** is searched.  
6) If every station in **StationLoader.closest_stations** is searched and none qualify, you get a notice to this effect -- and no data is loaded. 

# Notes on TMAX, TMIN, TMID
**StationLoader.station_data** loads directly the TMAX and TMIN values for every day. These are the actually measured temperature data values in the GHCND. "TAVG" values are included in GHCNd, but are actually  derived products that are equal to (TMAX + TMIN )/2 for every day. *(Note. TMID in this dataset is the midpoint between TMAX and TMIN for each day. Most weather databases and almanacs instead call TMID the "average temperature" for the day --including TAVG in GHCHNd! -- but this is inaccurate. Most weather stations only record maximum and minimum temperature for a day, and not the average, and so in fact when they report "average temperature" they really report "midpoint temperature". Don't you think calling the midpoint the average, is misleading? I do.)* 

But TAVG values are provided only sporadically in the GHCNd dataset. This is why StationLoader just calculates the "TMID" value (equivalent to the TAVG value), based upon TMAX and TMIN values. Therefore every single entry in self.station_data will have a TMAX, TMIN, and TMID value, and TMID will be equal to (TMAX+TMIN)/2 for every row. 
