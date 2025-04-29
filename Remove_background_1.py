import pandas as pd
import csv
import os
import sys
inFile = sys.argv[1]

reader = pd.read_csv(inFile, delimiter="\t", names=["chr","p", "position", "/", "strand", "/2", "base", "in_arrest", "in_rt","IP_arrest","IP_rt","in_sum","IP_sum","In_stop","IP_stop","peak","samp"], chunksize=1000000)
for df in reader:
    #print (r + len(df))
    df["position"] = pd.to_numeric(df["position"])
    for i in range(len(df)):
        if i-1 > 0 and df["position"].iloc[i] == df["position"].iloc[i-1]:
            continue
        elif i-1 > 0 and df["position"].iloc[i] != df["position"].iloc[i-1]:
            if i+1 < len(df) and df["position"].iloc[i] != df["position"].iloc[i+1]:
                print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["/"].iloc[i],df["strand"].iloc[i],df["/2"].iloc[i],df["base"].iloc[i],df["in_arrest"].iloc[i],df["in_rt"].iloc[i],df["IP_arrest"].iloc[i],df["IP_rt"].iloc[i],df["in_sum"].iloc[i],df["IP_sum"].iloc[i],df["In_stop"].iloc[i],df["IP_stop"].iloc[i],df["peak"].iloc[i],df["samp"].iloc[i])
            elif i+1 < len(df) and df["position"].iloc[i] == df["position"].iloc[i+1]:
                if i+2 < len(df) and df["position"].iloc[i] != df["position"].iloc[i+2]:
                    print(df["chr"].iloc[i+1],df["p"].iloc[i+1],df["position"].iloc[i+1],df["/"].iloc[i+1],df["strand"].iloc[i+1],df["/2"].iloc[i+1],df["base"].iloc[i+1],df["in_arrest"].iloc[i+1],df["in_rt"].iloc[i+1],df["IP_arrest"].iloc[i+1],df["IP_rt"].iloc[i+1],df["in_sum"].iloc[i+1],df["IP_sum"].iloc[i+1],df["IP_stop"].iloc[i+1],df["In_stop"].iloc[i+1],df["peak"].iloc[i+1],df["samp"].iloc[i+1]) 
                elif i+2 < len(df) and df["position"].iloc[i] == df["position"].iloc[i+2]:
                    if i+3 < len(df) and df["position"].iloc[i] != df["position"].iloc[i+3]:
                        print(df["chr"].iloc[i+2],df["p"].iloc[i+2],df["position"].iloc[i+2],df["/"].iloc[i+2],df["strand"].iloc[i+2],df["/2"].iloc[i+2],df["base"].iloc[i+2],df["in_arrest"].iloc[i+2],df["in_rt"].iloc[i+2],df["IP_arrest"].iloc[i+2],df["IP_rt"].iloc[i+2],df["in_sum"].iloc[i+2],df["IP_sum"].iloc[i+2],df["In_stop"].iloc[i+2],df["IP_stop"].iloc[i+2],df["peak"].iloc[i+2],df["samp"].iloc[i+2])

                    elif i+3 < len(df) and df["position"].iloc[i] == df["position"].iloc[i+3]:
                        print(df["chr"].iloc[i+3],df["p"].iloc[i+3],df["position"].iloc[i+3],df["/"].iloc[i+3],df["strand"].iloc[i+3],df["/2"].iloc[i+3],df["base"].iloc[i+3],df["in_arrest"].iloc[i+3],df["in_rt"].iloc[i+3],df["IP_arrest"].iloc[i+3],df["IP_rt"].iloc[i+3],df["in_sum"].iloc[i+3],df["IP_sum"].iloc[i+3],df["In_stop"].iloc[i+3],df["IP_stop"].iloc[i+3],df["peak"].iloc[i+3],df["samp"].iloc[i+3])
                        





        
