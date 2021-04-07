#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
from shutil import copyfile

#TODO:
# rename i
# zeile für zeile
# kommentare
# auf englisch übersetzen

fileList = ["couette1", "couette2"]
copyfile("couette.xml", "couette_bu.xml")    

for f in fileList:
    copyfile("simulation-test-templates/" + f + ".xml", "couette.xml")
    copyfile("simulation-test-templates/" + f + ".csv", "md-comparison-data.csv")
    ret_code = os.system("./test")
    print("Simulation exit code: " + str(ret_code))
    
    df_test = pd.read_csv('.csv', delimiter=';')
    df_values = pd.read_csv('md-comparison-data.csv', delimiter=';')
    df_comparison = pd.concat([df_values, df_test],axis='columns', keys=['Comparative values','Test results'])
    errorSumM, errorSumX, errorSumY, errorSumZ, i = 0

    for index1,row1 in df_comparison.iterrows(): 
        #get iteration of sample
        cycle1 = row1[0]

        #get cellIndex of sample
        x1 = row1[1]
        y1 = row1[2]
        z1 = row1[3]

        #get mass and momentum of sample
        mass1 = row1[7]
        impulseX1 = row1[8]
        impulseY1 = row1[9]
        impulseZ1 = row1[10]

        for index2, row2 in df_comparison.iterrows():
            #get iteration of comparison
            cycle2 = row2[12]

            #get cellIndex of sample
            x2 = row2[13]
            y2 = row2[14]
            z2 = row2[15]
            mass2 = row2[19]
            
            #get mass and momentum of sample
            impulseX2 = row2[20]
            impulseY2 = row2[21]
            impulseZ2 = row2[22]
        
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

    #TODO: Handle all cases of failure with automated email notification to author of responsible commit
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

