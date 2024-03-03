import pandas as pd
import csv
import sys
inFile = sys.argv[1]

reader = pd.read_csv(inFile, delimiter="\t", names=["chr","p", "position","strand", "base","rep1_III/IV_In_stop", "rep1_III/IV_in_sum", "rep1_III/IV_IP_sum","rep1_In_stop","rep1_IP_stop","rep1_in_sum","rep1_IP_sum","rep1_peak","rep1_samp","rep2_III/IV_In_stop", "rep2_III/IV_in_sum", "rep2_III/IV_IP_sum","rep2_In_stop","rep2_IP_stop","rep2_in_sum","rep2_IP_sum","rep2_peak","rep2_samp","In_avg","IP_avg"], chunksize=1000000)
#reader.ncolumns = ["chr", "p","position","strand", "base", "arrest", "readthrough", "IP_stop", "IP_sum"]
for df in reader:
    #print (r + len(df))
    for i in range(len(df)):
        if df["rep1_peak"].iloc[i] == "in":
            if df["rep2_peak"].iloc[i] == "in":
            
                if df["rep1_samp"].iloc[i] == "III":
                    if df["rep2_samp"].iloc[i] =="III":
                        if df["IP_avg"].iloc[i] > 0.2:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
                
                    elif df["rep2_samp"].iloc[i] == "III_IV":
                        if df["IP_avg"].iloc[i] > 0.18:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
                elif df["rep1_samp"].iloc[i] == "III_IV":
                    if df["rep2_samp"].iloc[i] =="III":
                        if df["IP_avg"].iloc[i] > 0.18:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
                
                    elif df["rep2_samp"].iloc[i] == "III_IV":
                        if df["IP_avg"].iloc[i] > 0.2:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
                
            elif df["rep2_peak"].iloc[i] == "out":
                if df["rep1_samp"].iloc[i] == "III":
                    if df["rep2_samp"].iloc[i] == "III":
                        if df["IP_avg"].iloc[i] > 0.3:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
                    elif df["rep2_samp"].iloc[i] == "III_IV":
                        if df["IP_avg"].iloc[i] > 0.18:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
                elif df["rep1_samp"].iloc[i] == "III_IV":
                    if df["rep2_samp"].iloc[i] == "III":
                        if df["IP_avg"].iloc[i] > 0.2:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
                    elif df["rep2_samp"].iloc[i] == "III_IV":
                        if df["IP_avg"].iloc[i] > 0.18:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
        elif df["rep1_peak"].iloc[i] == "out":            
            if df["rep2_peak"].iloc[i] == "in": 
                if df["IP_avg"].iloc[i] > 0.1:
                    if df["rep1_samp"].iloc[i] == "III":
                        if df["IP_avg"].iloc[i] > 0.2:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
                    elif df["rep1_samp"].iloc[i] == "III_IV":
                        if df["IP_avg"].iloc[i] > 0.18:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
                
            elif df["rep2_peak"].iloc[i] == "out":
                if df["rep1_samp"].iloc[i] == "III":
                    if df["rep2_samp"].iloc[i] == "III":
                        if df["IP_avg"].iloc[i] > 0.3:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
                    elif df["rep2_samp"].iloc[i] == "III_IV":
                        if df["IP_avg"].iloc[i] > 0.225:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
                elif df["rep1_samp"].iloc[i] == "III_IV":
                    if df["rep2_samp"].iloc[i] == "III":
                        if df["IP_avg"].iloc[i] > 0.225:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])
                    elif df["rep2_samp"].iloc[i] == "III_IV":
                        if df["IP_avg"].iloc[i] > 0.15:
                            print(df["chr"].iloc[i],df["p"].iloc[i],df["position"].iloc[i],df["strand"].iloc[i],df["base"].iloc[i],df["rep1_III/IV_stop"].iloc[i],df["rep1_III/IV_in_sum"].iloc[i],df["rep1_III/IV_IP_sum"].iloc[i],df["rep1_In_stop"].iloc[i],df["rep1_IP_stop"].iloc[i],df["rep1_in_sum"].iloc[i],df["rep1_IP_sum"].iloc[i],df["rep1_peak"].iloc[i],df["rep1_samp"].iloc[i],df["rep2_III/IV_stop"].iloc[i],df["rep2_III/IV_in_sum"].iloc[i],df["rep2_III/IV_IP_sum"].iloc[i],df["rep2_In_stop"].iloc[i],df["rep2_IP_stop"].iloc[i],df["rep2_in_sum"].iloc[i],df["rep2_IP_sum"].iloc[i],df["rep2_peak"].iloc[i],df["rep2_samp"].iloc[i],df["In_avg"].iloc[i],df["IP_avg"].iloc[i])