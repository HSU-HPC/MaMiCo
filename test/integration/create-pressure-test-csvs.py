#!/usr/bin/env python3

import csv

from numpy import random

with open('original_signal.csv', 'w', newline='') as f1, \
open('noisy_signal.csv', 'w', newline='') as f2, \
open('noisy_read.csv', 'w', newline='') as f3:

	csv_signal = csv.writer(f1, delimiter=';', quotechar=' ', quoting=csv.QUOTE_ALL)
	csv_noisy = csv.writer(f2, delimiter=';', quotechar=' ', quoting=csv.QUOTE_ALL)
	noisy_read = csv.writer(f3, delimiter=';', quotechar=' ', quoting=csv.QUOTE_ALL)

	for t in range(105):
		for z in range(6):
			for y in range(6):
				for x in range(6):
					signal = 1 + ((x+t)%2 - 0.5) / 10
					noisy = signal + random.normal(scale=0.2)
					csv_signal.writerow([t+1,x+4,y+4,z+4, "%.5f" % signal])
				#	    "%.5f" % random.normal(),"%.5f" % random.normal(),"%.5f" % random.normal()])
					csv_noisy.writerow([t+1,x+4,y+4,z+4, "%.5f" % noisy ])
#						"%.5f" % random.normal(),"%.5f" % random.normal(),"%.5f" % random.normal()])
					noisy_read.writerow([t+1,x,y,z,x+4,y+4,z+4, "%.5f" % noisy, \
						"%.5f" % random.normal(),"%.5f" % random.normal(),"%.5f" % random.normal()])
		csv_signal.writerow([])
		csv_noisy.writerow([])
		noisy_read.writerow([])
