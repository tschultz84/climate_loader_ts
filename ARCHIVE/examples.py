# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 13:10:53 2021

@author: 14154
"""



#%%
#First test: Loading weather in Bozeman, Montana.
bzdata=LoadStation([45.647256643331126,-111.04060494981753],True) #Loading in Bozeman, MT coordinates

#%%
#Analyze the climate here for the 10 years from 2011 to 2021.
date1 = '2011-1-1'
date2 = '2021-12-31'
bzcalc=StationAnalyzer(bzdata.station_data,date1,date2,display=True)

#%%
obdata=load.LoadStation([37.755663644,-122.506497974],True) #Loading in Ocean Beach, SF coordinates
#obdata.calculate_tmid_new(obdata.all_data_np)


#%%
date1 = '2019-1-1'
date2 = '2020-12-1'
obcalc=analyze.StationAnalyzer(obdata.station_data,date1,date2,display=True)

#print(obcalc.kpi)
#obcalc.key_charts()

#%%
sddata=load.LoadStation([32.741947,-117.239571],True) #Loading in Ocean Beach, San Diego coordinates
#%%
date1 = '2021-1-1'
date2 = '2021-12-31'
sdcalc=analyze.StationAnalyzer(sddata.station_data,date1,date2,display=True)

#%%

dvtdata=load.LoadStation([32.659167, -116.099167],True) #Loading in Desert View Tower data
#%%
date1 = '2021-1-1'
date2 = '2021-12-31'
dvtcalc=analyze.StationAnalyzer(dvtdata.station_data,date1,date2,display=True)

#%%

atldata=load.LoadStation([33.6639, -84.428],True) #Loading in ATlanta GA data
#%%
date1 = '1984-1-1'
date2 = '1984-12-31'
atlcalc=analyze.StationAnalyzer(atldata.station_data,date1,date2,display=True)




