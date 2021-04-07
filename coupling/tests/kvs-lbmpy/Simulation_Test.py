#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import kvstest
from shutil import copyfile

fileList = ["file1","file2","file3"]
copyfile("kvstest.ini", "kvstest_bu.ini")
copyfile("kvs.xml", "kvs_bu.xml")    

for f in fileList:
    copyfile("simulation-test-templates/" + f + ".ini", "kvstest.ini")
    copyfile("simulation-test-templates/" + f + ".xml", "kvstest.xml")
    try:
        os.run("kvstest.py")
    except:
        print("Simulation test crashed.")
    
    df_test = pd.read_csv('unfiltered1.csv', delimiter=';')
    df_values = pd.read_csv('unfiltered2.csv', delimiter=';')
    df_comparison = pd.concat([df_values, df_test],axis='columns', keys=['Comparative values','Test results'])
    print (df_comparison.head(2))
    errorSumM = 0
    errorSumX = 0
    errorSumY = 0
    errorSumZ = 0
    i = 0

    for index,row in df_comparison.iterrows(): 
        cycle1 = row[0]
        x1 = row[1]
        y1 = row[2]
        z1 = row[3]
        mass1 = row[7]
        impulseX1 = row[8]
        impulseY1 = row[9]
        impulseZ1 = row[10]
        rowCsv1 = index
        for index, row in df_comparison.iterrows():
            cycle2 = row[12]
            x2 = row[13]
            y2 = row[14]
            z2 = row[15]
            mass2 = row[19]
            impulseX2 = row[20]
            impulseY2 = row[21]
            impulseZ2 = row[22]
            rowCsv2 = index
        
            if cycle1 == cycle2 and x1 == x2 and y1 == y2 and z1 == z2:
                errorM = float(mass1) - float(mass2)
                errorX = float(impulseX1) - float(impulseX2)
                errorY = float(impulseY1) - float(impulseY2)
                errorZ = float(impulseZ1) - float(impulseZ2)
            
                errorSqM = errorM * errorM
                errorSqX = errorX * errorX
                errorSqY = errorY * errorY
                errorSqZ = errorZ * errorZ
            
                errorSumM = errorSumM + errorSqM
                errorSumX = errorSumX + errorSqX
                errorSumY = errorSumY + errorSqY
                errorSumZ = errorSumZ + errorSqZ
                i = i+1
            
    mseM = errorSumM/i
    mseX = errorSumX/i
    mseY = errorSumY/i
    mseZ = errorSumZ/i
    print ("mseM:",mseM,"mseX:",mseX,"mseY:",mseY,"mseZ:",mseZ,". Es wurden",i,"Werte verglichen.")
    if abs(mseM) > 1:           
        print ("MSE der Masse zu hoch. mseM:",mseM,"mseX:",mseX,"mseY:",mseY,"mseZ:",mseZ)
    if abs(mseX) > 1:           
        print ("MSE des Impulses in X-Richtung zu hoch.","mseM:",mseM,"mseX:",mseX,"mseY:",mseY,"mseZ:",mseZ)
    if abs(mseY) > 1:           
        print ("MSE des Impulses in Y-Richtung zu hoch.","mseM:",mseM,"mseX:",mseX,"mseY:",mseY,"mseZ:",mseZ)
    if abs(mseZ) > 1:           
        print ("MSE des Impulses in Z Richtung zu hoch.","mseM:",mseM,"mseX:",mseX,"mseY:",mseY,"mseZ:",mseZ)
    if abs(mseM) <= 1 and abs(mseX) <= 1 and abs(mseY) <= 1 and abs(mseZ) <= 1:
        print ("Simulation erfolgreich.")

