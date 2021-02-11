#!/usr/bin/env python3

import csv
from numpy import random

with open('original_signal.csv', 'w', newline='') as f1, \
open('noisy_signal.csv', 'w', newline='') as f2:

	csv_signal = csv.writer(f1, delimiter=';', quotechar=' ', quoting=csv.QUOTE_ALL)
	csv_noisy = csv.writer(f2, delimiter=';', quotechar=' ', quoting=csv.QUOTE_ALL)

	for t in range(60):
		for z in range(6):
			for y in range(6):
				for x in range(6):
					signal = (x+t)%2
					noisy = signal + random.normal(scale=2)
					csv_signal.writerow([t+1,x,y,z,x+3,y+3,z+3, "%.5f" % signal,\
					    "%.5f" % random.normal(),"%.5f" % random.normal(),"%.5f" % random.normal()])
					csv_noisy.writerow([t+1,x,y,z,x+3,y+3,z+3, "%.5f" % noisy, \
						"%.5f" % random.normal(),"%.5f" % random.normal(),"%.5f" % random.normal()])
		csv_signal.writerow([])
		csv_noisy.writerow([])