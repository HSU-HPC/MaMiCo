#!/usr/bin/env python3

import math
import matplotlib.pyplot as mplt

#YOU HAVE TO CUSTOMIZE THIS
number_procs = [1, 2, 3, 4, 5, 6, 7, 8]
runtimes = [29.616, 14.883, 9.975, 7.526, 6.170, 5.194, 4.449, 3.923]

mplt.style.use("seaborn")
#produce the plot
mplt.plot(number_procs, runtimes, color="red", label = "Total testing runtime")

#get labeling right
mplt.xlabel('Number of Processes')
mplt.ylabel('s')
mplt.legend()

mplt.tight_layout()
mplt.savefig("ut_runtime_plot.png")


