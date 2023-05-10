from os import sep
import numpy as np
import argparse
from enum import Enum
import sys
import csv
import pandas as pd
from numpy import genfromtxt

id = "46711217"

m = []
for n in np.arange(100)+1:
	df = pd.read_csv(''.join(["/scratch/users/ladmon/" + id + "/" + id + "_",str(n),"_m_gather.out"]), sep=' ',header=None)
	m.append(df.values)
array = np.hstack(np.array(m))
np.savetxt(id + '_magnetization_1000.csv', array, delimiter=",", fmt='%i')

e = []
for n in np.arange(100)+1:
	df = pd.read_csv(''.join(["/scratch/users/ladmon/" + id + "/" + id + "_",str(n),"_e_gather.out"]), sep=' ',header=None)
	e.append(df.values)
array = np.hstack(np.array(e))
np.savetxt(id + '_energy_1000.csv', array, delimiter=",", fmt='%i')