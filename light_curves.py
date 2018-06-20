import pandas as pd
import numpy as np
import csv
import math
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib
import sys

def date_time(x):
    return datetime.strptime(sources['DATETIME'][x],'%Y-%m-%d %H:%M:%S.%f').timestamp()


with open('sources.csv') as csvFile:
    reader = csv.reader(csvFile)
    keys = next(reader)

sources = dict()
for i in keys:
    df = pd.read_csv('sources.csv')
    # sources[i]=np.array(df[i])
    test = np.append(np.array(df[i]),np.array(df[i]))
    sources[i] = test

# 726-053334
# identifier = input("Enter UCAC4 ID or 'C' to input coordinates instead: ")
identifier = '726-053334'

if identifier=='C':
    RA = float(input("RA coordinates in decimal decrees: "))
    DEC = float(input("DEC coordinates in decimal decrees: "))

if not identifier=='C':
    indices = np.nonzero(sources['id']==identifier)[0]
    # data = []
    # for x in indices:
    #     data.append([sources[i][x] for i in keys])
    print('%s data points found' % len(indices))


filt = input("Filter: ")
print("Enter dates as yyyy/mm/dd/hh/mm in GMT and 24 hr time or enter 'all' to see all data")
start = input("Start time : ")
if not start=='all':
    end = input("End time: ")
    start = datetime.strptime(start,'%Y/%m/%d/%H/%M').timestamp()
    end = datetime.strptime(end,'%Y/%m/%d/%H/%M').timestamp()

mags,time = [],[]
for x in indices:
    mags.append(sources['MAG_'+filt][x])
    time.append(datetime.strptime(sources['DATETIME'][x], '%Y-%m-%d %H:%M:%S.%f').timestamp())

mags = [mags[x] for x in range(len(mags)) if isinstance(mags[x], float)]
try:
    time = [time[x] for x in range(len(time)) if isinstance(mags[x], float)]
except IndexError:
    print('No %s magnitude data found' % filt)
    sys.exit()

plt.figure()
plt.scatter(time, mags)
if not start=='all':
    plt.xlim(start,end)
plt.show()
