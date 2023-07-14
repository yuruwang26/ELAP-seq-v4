import pandas as pd
import csv
import os
import sys
inFile = sys.argv[1]

reader = pd.read_csv(inFile, delimiter="\t", names=["chr","p", "position", "/", "strand", "/2", "base", "in_arrest", "in_rt","IP_arrest","IP_rt","in_sum","IP_sum","In_stop","IP_stop","peak","samp"], chunksize=1000000)
#reader.ncolumns = ["chr", "p","position","strand", "base", "arrest", "readthrough", "IP_stop", "IP_sum"]
for df in reader:
    #print (r + len(df))
    df["position"] = pd.to_numeric(df["position"])
    for i in range(len(df)):
        dis1 = df["position"].iloc[i] - 35
        dis2 = df["position"].iloc[i] + 35
        for t in range(1,35):
            if i-t > 0 and df["position"].iloc[i-t] > dis1 and df["position"].iloc[i-t] < df["position"].iloc[i]:
                read_difference_1 = (df["IP_sum"].iloc[i-t])/(df["IP_sum"].iloc[i])
                if read_difference_1 > 3.5:
                    print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["/"].iloc[i],df["strand"].iloc[i],df["/2"].iloc[i],df["base"].iloc[i],df["in_arrest"].iloc[i], df["in_rt"].iloc[i],df["IP_arrest"].iloc[i],df["IP_rt"].iloc[i], df["in_sum"].iloc[i], df["IP_sum"].iloc[i],df["In_stop"].iloc[i],df["IP_stop"].iloc[i],df["peak"].iloc[i],df["samp"].iloc[i]  )
                   
                    exit
            elif i+t < len(df) and df["position"].iloc[i+t] < dis2 and df["position"].iloc[i+t] > df["position"].iloc[i]:     
                read_difference_2 = (df["IP_sum"].iloc[i+t])/(df["IP_sum"].iloc[i])
                if read_difference_2 > 3.5:
                    print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["/"].iloc[i],df["strand"].iloc[i],df["/2"].iloc[i],df["base"].iloc[i],df["in_arrest"].iloc[i], df["in_rt"].iloc[i],df["IP_arrest"].iloc[i],df["IP_rt"].iloc[i], df["in_sum"].iloc[i], df["IP_sum"].iloc[i],df["In_stop"].iloc[i],df["IP_stop"].iloc[i],df["peak"].iloc[i],df["samp"].iloc[i]  )
                   
                    exit
       