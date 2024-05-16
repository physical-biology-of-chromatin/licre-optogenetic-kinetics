#! /usr/bin/python3
#
# USAGE: nohup ./progressLogging.py &
#
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import least_squares
import time
import sys
from datetime import datetime

Niter=20
eps=10

logfile = open('task.log', 'a', encoding="utf-8")
startime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
l = "\n" + 'STARTING JOB AT ' + startime + "\n"
logfile.write(l)
logfile.close()

for it in range(Niter):
    logfile = open('task.log', 'a', encoding="utf-8")
    l = 'I try a string with ' + repr(it) + ' and ' + "kk \n"
    time.sleep(1.5)
    logfile.write(l)
    logfile.close()

logfile = open('task.log', 'a', encoding="utf-8")
endtime =  datetime.now().strftime('%Y-%m-%d %H:%M:%S')
l = 'ENDING JOB AT ' + endtime + "\n"
logfile.write(l)
logfile.close()

