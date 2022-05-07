#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 07:11:53 2020

@author: scs-rd
"""

#%%


#%%

#First isolates stations which start their record prior to 1950
early = allstdata[allstdata['Firstyear']<1920]
#The further limits to only stations  which have records  which go through 2015
long = early[early['Lastyear']>2015]

short = allstdata[allstdata['Firstyear']>2000]

short = short[short['Element'].isin(["TAVG","TMAX"])]


#%%
#TEST ANALYSES -- MONTANA, CALIFORNIA
#The c
#%%
mtstations = tempdata[tempdata['STATE']=="MT"]

mtstationspart=mtstations[(mtstations['Firstyear']<1941) ]
mtstationspart=mtstationspart[(mtstationspart['Lastyear']>1998) ]

#%%
years1 =np.concatenate((np.arange(1900,1920),np.arange(1999,2019)))

#mtstationsd = StationData(mtstationspart,["TMEAN"],years1,True)
#mtstationsd.load_stats()

#%%
t1 = AnalyzeStationData(mtstationsd,20,"TMEAN",[1900,1920],[1999,2019],18,18,False)

statdata=t1.all_stats()
#t2 = AnalyzeStationData(mtstationsd,24,"TAVG",[1900,1920],[2000,2020],5,5)
#%%
t1.compdata('USW00024135')
t1.compdata('USC00242689')
#%%
t1.compdata('USC00244522')
t1.compdata('USC00248597')
t1.compdata('USC00240056')


#%%
t1.simpledatahist(['USC00240364'])
return5=t1.perioddata(['USC00240364'])
#%%

years2 =np.arange(1880,1950)
mtstationsd = StationData(mtstationspart.iloc[0:2],["TMEAN"],years2,True)
#mtstationsd2 = StationData(mtstationspart.iloc[0:2],["TMAX","TMIN"],years2,True)

data=mtstationsd.load_stats()
#data2=mtstationsd2.load_stats()
#%%
#Test the average
month=6
year=1900

x=data[(data['Year']==year)&(data['Month']==month)].iloc[:,4:]
y=data2[(data2['Year']==year)&(data2['Month']==month)].iloc[:,4:].mean()
pd.set_option('display.max_columns', 6)
test =x.divide(y)
print(test)


#%%
#ANALYSIS -- US
#Selects stations in USA
countries = allstdata['ID'].str[:2]
usa=allstdata[countries=="US"]
#%%
years10 = np.arange(1860,2020)

usadata = StationData(usa.sample(150),["TMAX"],years10,True)
r555=usadata.load_stats()

