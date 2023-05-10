from os import sep
import numpy as np
import argparse
from enum import Enum
import sys
import matplotlib.pyplot as plt
import csv

id = "46658455"
RPT = 1
rpt = 10 #repetition in simulation
# array = np.zeros([50,RPT*rpt*100])
# for t in range(50):
#     print(".",flush=True,end='')
#     for n in np.arange(100)+1:
#         for i in range(rpt):
#             try: #ignore faulty simulations
#                 with open(''.join(["/scratch/users/ladmon/"+ id + "/" + id + "_",str(n),",i=",str(i),",t=",str(t),",energy.out"]), 'r') as f:
#                     for r in range(1,RPT+1):
#                         firstline = f.readline().rstrip()
#                         array[t,(r-1)*rpt*100 + i*100 + n-1] = int(firstline)
#             except: 
#                 for r in range(1,RPT+1):
#                     array[t,(r-1)*rpt*100 + i*100 + n-1] = int(0)
# np.savetxt(id + '_energy_1000.csv', array, delimiter=",", fmt='%i')

array = np.zeros([50,RPT*rpt*100])
for t in range(50):
    print(".",flush=True,end='')
    for n in np.arange(100)+1:
        for i in range(rpt):
            try:
                with open(''.join(["/scratch/users/ladmon/" + id + "/" + id + "_",str(n),",i=",str(i),",t=",str(t),",magnetization.out"]), 'r') as f:
                    for r in range(1,RPT+1):
                        firstline = f.readline().rstrip()
                        array[t,(r-1)*rpt*100 + i*100 + n-1] = int(firstline)
            except: 
                for r in range(1,RPT+1):
                    array[t,(r-1)*rpt*100 + i*100 + n-1] = int(0)
np.savetxt(id + '_magnetization_1000.csv', array, delimiter=",", fmt='%i')

#"/scratch/users/ladmon/44901342_1,i=1,energy.out"