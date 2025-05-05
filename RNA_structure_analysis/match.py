import pandas as pd
import csv
import os
import sys
import numpy as np

file = open('5‘UTR-input.txt','r') 
Lines = file.readlines()
for line in Lines:
    m = 0
    a = [0]*200
    for i in range(1,200):
        if line[i] == ".":
            m = m +0
            if i == 100:
                print("unpair")
        elif line[i] != ".":
            if line[i] == "(":
                m = m+1
                a[m] = i
            elif line[i] == ")":
                if a[m] == 100 or i == 100:
                    print(a[m], i)
                m = m-1
     

 
