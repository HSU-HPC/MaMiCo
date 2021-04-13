#!/usr/bin/env python
# coding: utf-8
#@author: Leonard Hannen, Felix Maurer

import os
import pandas as pd
import smtplib
from shutil import copyfile
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

#######################################

#Accepted mean squared error for mass
ACCEPTED_MSE_M = 1

#Accepted mean squared error for momentum in x, y and z direction
ACCEPTED_MSE_X = 0
ACCEPTED_MSE_Y = 1
ACCEPTED_MSE_Z = 1

#sets up the email
sender = 'sender@example.com'
recievers = ["person1@example.com", "person2@example.com", "person3@example.com"]

#######################################

msg = MIMEMultipart()
body = "While running the couette test, following errors occurred:\n"
msg['Subject'] = 'Simulation Test error report'
msg['From'] = sender
msg['To'] = ','.join(recievers)

#creates backup of the original .xml file
fileList = ["kvs1", "kvs2"]
copyfile("kvs.xml", "kvs_bu.xml")
copyfile("kvstest.ini", "kvstest_bu.ini")

errorCheck = False

for f in fileList:
    #copys the templates into the right folder
    copyfile("simulation-test-templates/" + f + ".xml", "kvs.xml")
    copyfile("simulation-test-templates/" + f + ".ini", "kvs.ini")
    copyfile("simulation-test-templates/" + f + ".csv", "md-comparison-data.csv")

    #starts the test
    ret_code = os.system("python3 kvstest.py")
    print("Simulation exit code: " + str(ret_code))

    if(ret_code != 0):
        body += "Simulation with config " + f + " exited with nonzero exit code: " + str(ret_code) + "\n"
        errorCheck = True
        continue
    
    #reads in the .csv files
    df_test = pd.read_csv('md-sample.csv', delimiter=';')
    df_values = pd.read_csv('md-comparison-data.csv', delimiter=';')

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

        #get mass and momentum of sample
        mass1 = row[7]
        momentumX1 = row[8]
        momentumY1 = row[9]
        momentumZ1 = row[10]

        #get cellIndex of comparison
        x2 = row[13]
        y2 = row[14]
        z2 = row[15]

        #get mass and momentum of comparison
        mass2 = row[19]
        momentumX2 = row[20]
        momentumY2 = row[21]
        momentumZ2 = row[22]
        
        #gets the differences between sample and comparison
        errorM = float(mass1) - float(mass2)
        errorX = float(momentumX1) - float(momentumX2)
        errorY = float(momentumY1) - float(momentumY2)
        errorZ = float(momentumZ1) - float(momentumZ2)
    
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
    mseXStr = str(mseX)
    mseYStr = str(mseY)
    mseZStr = str(mseZ)

    #the value that determines the success or failure of a simulation can be changed here
    if abs(mseM) > ACCEPTED_MSE_M:        
        body += "Using configuration "+fStr+ ", the MSE of the mass was too high.\nmseM:"+mseMStr+ " mseX:"+mseXStr+" mseY:"+mseYStr+" mseZ:"+mseZStr+ ".There were "+compStr+ " comparisons made.\n\n"
        errorCheck = True
    if abs(mseX) > ACCEPTED_MSE_X:           
        body += "Using configuration "+fStr+", the MSE of the momentum in x-direction was too high.\nmseM:" +mseMStr+ " mseX:"+mseXStr+" mseY:"+mseYStr+" mseZ:"+mseZStr+ ".There were "+compStr+ " comparisons made.\n\n"
        errorCheck = True
    if abs(mseY) > ACCEPTED_MSE_Y:           
        body += "Using configuration "+fStr+ ", the MSE of the momentum in y-direction was too high.\nmseM:" +mseMStr+ " mseX:"+mseXStr+" mseY:"+mseYStr+" mseZ:"+mseZStr+ ".There were "+compStr+ " comparisons made.\n\n"
        errorCheck = True
    if abs(mseZ) > ACCEPTED_MSE_Z:           
        body += "Using configuration "+fStr+ ", the MSE of the momentum in z-direction was too high.\nmseM:" +mseMStr+ " mseX:"+mseXStr+" mseY:"+mseYStr+" mseZ:"+mseZStr+ ".There were "+compStr+ " comparisons made.\n\n"
        errorCheck = True
    if abs(mseM) <= ACCEPTED_MSE_M and abs(mseX) <= ACCEPTED_MSE_X and abs(mseY) <= ACCEPTED_MSE_Y and abs(mseZ) <= ACCEPTED_MSE_Z:
        print ("Simulation successful.")

print(body)

msg.attach(MIMEText(body, 'plain'))
msgText = msg.as_string()

#sends the email
mail=smtplib.SMTP('host')
if errorCheck:
    mail.sendmail(sender, recievers, msgText)
    mail.quit
