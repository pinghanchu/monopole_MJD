#!/usr/bin/env python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os.path
from os import path

def fill(mass, energy):
    energyfile = "{}GeV_{}GeV".format(mass,energy)
    filename = "../data/data_{}".format(energyfile)

    df=pd.read_csv("{}.csv".format(filename))
    
    dfnew = pd.DataFrame(columns=['EventID','Channel','StartTime','EndTime','StartPx','EndPx','StartPy','EndPy','StartPz','EndPz','Edep','EventStartTime','EventEndTime','EventStartPx','EventEndPx','EventStartPy','EventEndPy','EventStartPz','EventEndPz','Esum',"M"]) 
    
    for ievent in df["EventID"].unique():
        df1=df[df["EventID"]==ievent]
        m = len(df1["Channel"].unique())
        eventstarttime=0
        eventendtime=0
        eventstartpx = df1["DetPX"].iloc[0] # m
        eventendpx = df1["DetPX"].iloc[len(df1)-1]
        eventstartpy = df1["DetPY"].iloc[0]
        eventendpy = df1["DetPY"].iloc[len(df1)-1]
        eventstartpz = df1["DetPZ"].iloc[0]
        eventendpz = df1["DetPZ"].iloc[len(df1)-1]
        Esum =0.
        for j in range(len(df1)):
            Esum = Esum + float(df1["Edep"].iloc[j]) # MeV
            count = 0
        for ichannel in df1["Channel"].unique():
            df2=df1[df1["Channel"]==ichannel]
            starttime = float(df2["Time"].iloc[0])*1e-9
            endtime = float(df2["Time"].iloc[len(df2)-1])*1e-9
            if(count==0):
                eventstarttime=endtime
            elif(count==m-1):
                eventendtime=endtime
            count=count+1
                    
        for ichannel in df1["Channel"].unique():
            df2=df1[df1["Channel"]==ichannel]
            Edep = 0.
            starttime = float(df2["Time"].iloc[0])*1e-9
            endtime = float(df2["Time"].iloc[len(df2)-1])*1e-9
            startpx = df2["DetPX"].iloc[0]
            endpx = df2["DetPX"].iloc[len(df2)-1]
            startpy = df2["DetPY"].iloc[0]
            endpy = df2["DetPY"].iloc[len(df2)-1]
            startpz = df2["DetPZ"].iloc[0]
            endpz = df2["DetPZ"].iloc[len(df2)-1]  
                        
            for i in range(len(df2)):
                Edep = Edep + float(df2["Edep"].iloc[i])
                            
            result = {'EventID':ievent,'Channel':ichannel,'StartTime':starttime,'EndTime':endtime,'StartPx':startpx,'EndPx':endpx,'StartPy':startpy,'EndPy':endpy,'StartPz':startpz,'EndPz':endpz,'Edep':Edep,'EventStartTime':eventstarttime,'EventEndTime':eventendtime,'EventStartPx':eventstartpx,'EventEndPx':eventendpx,'EventStartPy':eventstartpy,'EventEndPy':eventendpy,'EventStartPz':eventstartpz,'EventEndPz':eventendpz,'Esum':Esum,"M":m}
            
            dfnew=dfnew.append(result,ignore_index=True)
                            
    dfnew["DeltaT"]=dfnew["EndTime"]-dfnew["StartTime"]
    dfnew["EventDeltaT"]=dfnew["EventEndTime"]-dfnew["EventStartTime"]
    dfnew["EventDeltaL"]=np.sqrt((dfnew["EventEndPx"]-dfnew["EventStartPx"])*(dfnew["EventEndPx"]-dfnew["EventStartPx"])+(dfnew["EventEndPy"]-dfnew["EventStartPy"])*(dfnew["EventEndPy"]-dfnew["EventStartPy"])+(dfnew["EventEndPz"]-dfnew["EventStartPz"])*(dfnew["EventEndPz"]-dfnew["EventStartPz"]))
    dfnew["EventVelocity"]=dfnew["EventDeltaL"]/dfnew["EventDeltaT"]
    
    dfnew.to_hdf("../data/event_{}.h5".format(energyfile),key='df',mode='w')
    dfnew.to_csv("../data/event_{}.csv".format(energyfile))
             
    return dfnew

for j in range(10,31):
    mass="1e{}".format(j)
    for i in range(-5,11):
        energy = "1e{}".format(i)
        print(mass, energy)
        fill(mass,energy)
