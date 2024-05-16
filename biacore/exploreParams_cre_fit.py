#! /usr/bin/python3
#
# USAGE: nohup ./thisscript.py &
#
# Progress can be monitored by typing:
#
#      tail exploreParams_progress.log
#
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import least_squares
import time
import sys
from datetime import datetime

##########################
#
# PROGRESS MONITORING
#
##########################
# Log file to monitor progress
# (it is important to open and close it
#  at every step so that content is flushed in it
#  before the end of the run)
logfilename = 'exploreParams_cre_fit_progress.log'
logfile = open(logfilename, 'a', encoding="utf-8")
startime = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
msg = "\n" + 'STARTING JOB AT ' + startime + "\n"
logfile.write(msg)
logfile.close()

##########################
#
#  MODEL IMPLEMENTATION 
#
##########################

# differential equation of the half and full loxP system
def halfb(y,t,param): #y=(C,B_1)
    alp,C0,ta,RUmax,kt,k1,km1=param
    dy=[0,0] #(dC/dt,dB_1/dt)
    B0=RUmax-y[1]
    if t<ta:
        dy[0]=kt*(C0-y[0])-k1/alp*y[0]*B0+km1/alp*y[1]
        dy[1]=k1*y[0]*B0-km1*y[1]
    else:
        dy[0]=kt*(-y[0])-k1/alp*y[0]*B0+km1/alp*y[1]
        dy[1]=k1*y[0]*B0-km1*y[1]
    return dy

def full1b(y,t,param): #y=(C,B1,B2)
    alp,C0,ta,RUmax,kt,k1,km1,k2,km2=param
    dy=[0,0,0] #dC/dt,dB_1/dt,dB_2/dt
    B0=RUmax/2-y[1]-y[2]
    if t<ta:
        dy[0]=kt*(C0-y[0])-2*(k1/alp)*y[0]*B0+(km1/alp)*y[1]-(k2/alp)*y[0]*y[1]+2*(km2/alp)*y[2]
        dy[1]= 2*k1*y[0]*B0-km1*y[1]-k2*y[0]*y[1]+2*km2*y[2]  
        dy[2]= k2*y[0]*y[1]-2*km2*y[2] 
    else:
        dy[0]=kt*(-y[0])-2*k1/alp*y[0]*B0+km1/alp*y[1]-k2/alp*y[0]*y[1]+2*km2/alp*y[2]
        dy[1]= 2*k1*y[0]*B0-km1*y[1]-k2*y[0]*y[1]+2*km2*y[2]  
        dy[2]= k2*y[0]*y[1]-2*km2*y[2] 
    return dy

##########################
#
# FITTING FUNCTION 
#
#########################

#infer based on relative
def fitallrel(param,xdatah,xdataf,ydatah,ydataf):
    Nconch,Lconch,Ntimeh,Ltimeh,tah,epsh=xdatah
    Nconcf,Lconcf,Ntimef,Ltimef,taf,epsf=xdataf
    alp,RUmaxh,RUmaxf,kt,k1,km1,k2,km2=param
    yh=np.zeros((Ntimeh,Nconch))
    for c in range(Nconch):
        yh[:,c]=odeint(halfb,[0,0],Ltimeh,args=([alp,Lconch[c],tah,RUmaxh,kt,k1,km1],))[:,1]  
    yf=np.zeros((Ntimef,Nconcf))
    for c in range(Nconcf):
        s=odeint(full1b,[0,0,0],Ltimef,args=([alp,Lconcf[c],taf,RUmaxf,kt,k1,km1,k2,km2],))
        yf[:,c]=s[:,1]+2*s[:,2]
    return np.concatenate(((yh-ydatah)/(ydatah+epsh),(yf-ydataf)/(ydataf+epsf)),axis=1).flatten()


###########################
#
# EXPERIMENTAL DATA
#
###########################

logfile = open(logfilename, 'a', encoding="utf-8")
msg = '...Loading experimental data\n'
logfile.write(msg)
logfile.close()

#load biacore data for Cre and extract data
itimemin=601 #index to start the fit (corresponding to t=0)
itimemax=4001 #index to end the fit
idt=10 #index of the step between two fitted points

tah=90 #in sec
dtsh=np.loadtxt('../../data/2023-06-01/biacore_cre_half_loxP.txt')
Nconch=len(np.unique(dtsh[:,1]))
Lconch=np.zeros(Nconch)
for c in range(1,Nconch+1):
    a=np.argwhere(dtsh[:,1]==c)
    Lconch[c-1]=dtsh[a[1],0] #in nM 
Ltimeh=dtsh[a,2]
Ltimeh=Ltimeh[np.arange(itimemin,itimemax,idt)]
Ntimeh=len(Ltimeh)

taf=90 #in sec
dtsf=np.loadtxt('../../data/2023-06-01/biacore_cre_full_loxP.txt')
Nconcf=len(np.unique(dtsf[:,1]))
Lconcf=np.zeros(Nconcf)
for c in range(1,Nconcf+1):
    a=np.argwhere(dtsf[:,1]==c)
    Lconcf[c-1]=dtsf[a[1],0] #in nM 
Ltimef=dtsf[a,2]
Ltimef=Ltimef[np.arange(itimemin,itimemax,idt)]
Ntimef=len(Ltimef)
               
Datah=np.zeros((Ntimeh,Nconch))
for c in range(Nconch):
    a=np.argwhere(dtsh[:,1]==(c+1))
    Datah[:,c]=dtsh[a[np.arange(itimemin,itimemax,idt)],3].ravel()
Dataf=np.zeros((Ntimef,Nconcf))
for c in range(Nconcf):
    a=np.argwhere(dtsf[:,1]==(c+1))
    Dataf[:,c]=dtsf[a[np.arange(itimemin,itimemax,idt)],3].ravel()


#############################################
#
# Fit using various initial parameter values
#
#############################################

Niter=2
eps=10
## For Cre:
# RUmaxh = 0.033*32*38*1 ~ 40
# RUmaxf = 0.033*11*38*2 ~ 27
# alph ~ 1.5
p_inf=[1,  36,24,0.08,-3,-3,-4,-4] #alp,RUmaxh,RUmaxf,kt,k1,km1,k2,km2
p_sup=[2.5,44,30,1.1 , 2, 2, 2, 0]
islog=[0,0,0,0,1,1,1,1] #0 if linear scale, 1 if log10 scale
p_init=np.zeros((Niter,8))+np.nan
score =np.zeros(Niter)+np.nan
p_optim=np.zeros((Niter,8))+np.nan

# Perform the fit
for it in range(Niter):

    # Progress message
    logfile = open(logfilename, 'a', encoding="utf-8")
    msg = '...Performing fit iteration ' + repr(it) + ' of ' + repr(Niter) + "\n"
    logfile.write(msg)
    logfile.close()

    # Generate random initial parameters within boundaries
    for i in range(8):
        if islog[i]==0:
            p_init[it,i]=np.random.uniform(p_inf[i], p_sup[i])
        else:
            p_init[it,i]=10.**(np.random.uniform(p_inf[i], p_sup[i]))

    # Fit model starting from these initial params
    pfit=least_squares(fitallrel,p_init[it,:],jac='3-point',method='trf',bounds=(0,np.inf),args=([Nconch,Lconch,Ntimeh,Ltimeh.ravel(),tah,eps],[Nconcf,Lconcf,Ntimef,Ltimef.ravel(),taf,eps],Datah,Dataf))
    score[it]=pfit.cost/(Ntimeh*Nconch+Ntimef*Nconcf)
    p_optim[it,:]=pfit.x


##########################################
#
# Save estimated parameters and scores
#
##########################################
outfileroot = "exploreParams_cre_fit_results_"
np.savetxt(outfileroot + 'p_optim.txt', p_optim)
np.savetxt(outfileroot + 'score.txt', score)
np.savetxt(outfileroot + 'p_init.txt', p_init)
np.savetxt(outfileroot + 'eps.txt', [eps])


############################################
# Write completion time in progress logfile
logfile = open(logfilename, 'a', encoding="utf-8")
endtime =  datetime.now().strftime('%Y-%m-%d %H:%M:%S')
msg = 'ENDING JOB AT ' + endtime + "\n"
logfile.write(msg)
logfile.close()
            

