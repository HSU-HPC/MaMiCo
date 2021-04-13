#!/usr/bin/env python
# coding: utf-8
#@author: Leonard Hannen, Felix Maurer

import os
import pandas as pd
import kvstest
import smtplib
from shutil import copyfile
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

#creates backup of the original .ini and .xml files
fileList = ["file1","file2","file3"]
copyfile("kvstest.ini", "kvstest_bu.ini")
copyfile("kvs.xml", "kvs_bu.xml")    

#sets up the email
msg = MIMEMultipart()
sender = 'sender@example.com'
recievers = ["person1@example.com", "person2@example.com", "person3@example.com"]
body = "While running the kvs test, following errors occurred:\n"

msg['Subject'] = 'Simulation Test error report'
msg['From'] = sender
msg['To'] = ','.join(recievers)

errorCheck = false

for f in fileList:
    #copys the templates into the right folder
    copyfile("simulation-test-templates/" + f + ".ini", "kvstest.ini")
    copyfile("simulation-test-templates/" + f + ".xml", "kvstest.xml")

    #starts the test
    try:
        os.run("kvstest.py")
    except:
        print("Simulation test crashed.")
    
    #reads in the .csv files
    df_test = pd.read_csv('unfiltered1.csv', delimiter=';')
    df_values = pd.read_csv('unfiltered2.csv', delimiter=';')

    #adds the data from both files into one dataframe
    df_comparison = pd.concat([df_values, df_test],axis='columns', keys=['Comparative values','Test results'])
    errorSumM = 0
    errorSumX = 0
    errorSumY = 0
    errorSumZ = 0
    #number of comparisons
    comp = 0

    for index,row in df_comparison.iterrows(): 

        #checks if a row is empty and skips it if true
        isNull = pd.isnull(df_comparison.iloc[index,1]) 
        if isNull:
            continue

        #get cellIndex of sample
        x1 = row[1]
        y1 = row[2]
        z1 = row[3]

        #get mass and impulse of sample
        mass1 = row[7]
        impulseX1 = row[8]
        impulseY1 = row[9]
        impulseZ1 = row[10]

        #get cellIndex of comparison
        x2 = row[13]
        y2 = row[14]
        z2 = row[15]

        #get mass and impulse of comparison
        mass2 = row[19]
        impulseX2 = row[20]
        impulseY2 = row[21]
        impulseZ2 = row[22]
        
        #gets the differences between sample and comparison
        errorM = float(mass1) - float(mass2)
        errorX = float(impulseX1) - float(impulseX2)
        errorY = float(impulseY1) - float(impulseY2)
        errorZ = float(impulseZ1) - float(impulseZ2)
    
        #squares the differences
        errorSqM = errorM * errorM
        errorSqX = errorX * errorX
        errorSqY = errorY * errorY
        errorSqZ = errorZ * errorZ
    
        #sums up the squared errors
        errorSumM = errorSumM + errorSqM
        errorSumX = errorSumX + errorSqX
        errorSumY = errorSumY + errorSqY
        errorSumZ = errorSumZ + errorSqZ
        
        #counts the comparisons
        comp = comp+1
    
    #get the MSE
    mseM = errorSumM/comp
    mseX = errorSumX/comp
    mseY = errorSumY/comp
    mseZ = errorSumZ/comp

    #turns values into strings
    fStr = str(f)
    compStr = str(comp)
    mseMStr = str(mseM)
    mseMStr = str(mseM)
    mseMStr = str(mseM)
    mseMStr = str(mseM)

    #the value that determines the success or failure of a simulation can be changed here
    if abs(mseM) > 1:        
        body += "Using configuration "+fStr+ ", the MSE of the mass was to high.\nmseM:"+mseMStr+ " mseX:"+mseXStr+" mseY:"+mseYStr+" mseZ:"+mseZStr+ ".There were "+compStr+ " comparisons made.\n\n"
        errorCheck = true
    if abs(mseX) > 1:           
        body += "Using configuration "+fStr+", the MSE of the impulse in x-direction was to high.\nmseM:" +mseMStr+ " mseX:"+mseXStr+" mseY:"+mseYStr+" mseZ:"+mseZStr+ ".There were "+compStr+ " comparisons made.\n\n"
        errorCheck = true
    if abs(mseY) > 1:           
        body += "Using configuration "+fStr+ ", the MSE of the impulse in y-direction was to high.\nmseM:" +mseMStr+ " mseX:"+mseXStr+" mseY:"+mseYStr+" mseZ:"+mseZStr+ ".There were "+compStr+ " comparisons made.\n\n"
        errorCheck = true
    if abs(mseZ) > 1:           
        body += "Using configuration "+fStr+ ", the MSE of the impulse in z-direction was to high.\nmseM:" +mseMStr+ " mseX:"+mseXStr+" mseY:"+mseYStr+" mseZ:"+mseZStr+ ".There were "+compStr+ " comparisons made.\n\n"
        errorCheck = true
    if abs(mseM) <= 1 and abs(mseX) <= 1 and abs(mseY) <= 1 and abs(mseZ) <= 1:
        print ("Simulation successful.")

msg.attach(MIMEText(body, 'plain‘))
msgText = msg.as_string()

#sends the email
mail=smtplib.SMTP('host')
if errorCheck:
    mail.sendmail(sender, recievers, msgText)
    mail.quit