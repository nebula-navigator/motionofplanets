#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  5 14:29:09 2025

@author: sohaib
"""

import numpy as np
from tabulate import tabulate
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd


def leapfrog(dt,steps,x0,y0,vx0,vy0):


    r0=np.sqrt(x0**2+y0**2) #radial position of planet at x0,y0
    print("r0 :",r0)
    
    ax0=-x0/r0**3 #acceleration component at x0
    ay0=-y0/r0**3 #acceleration component at y0
    print("ax(0) :",ax0)
    print("ay(0) :",ay0)
    
    vxhalf=vx0+ax0*(dt/2) #half integer time step velocity in x direction
    vyhalf=vy0+ay0*(dt/2) #half integer time step velocity in y direction
    print("vx(e/2) :",vxhalf)
    print("vy(e/2) :",vyhalf)
    
    time = 0.0 #started at time 0
    
    #bunch of initializations for arrays to store data later
    
    results=[]
    
    results.append([time, x0, y0, ax0, ay0, vx0, vy0,r0])
    
    
    #just to match the table in Feynman's book, I included an additional row in the data table to show half step velocities
    results.append([None, None, None, None, None, vxhalf, vyhalf,None]) 
    
    xcurr = x0
    ycurr = y0
    vxhalf_curr = vxhalf
    vyhalf_curr = vyhalf
    
    #looping and calculating parameters until time interval reaches the limit
    
    for n in range(steps + 1):
        xnew= xcurr + vxhalf_curr*dt
        ynew= ycurr + vyhalf_curr*dt
        
        time+=dt
        
        rnew= np.sqrt(xnew**2+ynew**2)
        
        
        axnew=-xnew/rnew**3
        aynew=-ynew/rnew**3
        
        vxnew=vxhalf_curr + axnew*(dt/2)
        vynew=vyhalf_curr + aynew*(dt/2)
        
        results.append([time,xnew, ynew, axnew, aynew, vxnew, vynew,rnew])
        
        vxhalf_curr= vxnew + (dt/2)*axnew
        vyhalf_curr = vynew + (dt/2)*aynew
        
        results.append([None, None, None, None, None, vxhalf_curr, vyhalf_curr,None])
        
        xcurr = xnew
        ycurr = ynew
        
        
        
 
    #to print the table
    
    headers = ["time","x", "y", "ax", "ay", "vx", "vy", "r"]
    print(tabulate(
        results,
        headers=headers,
        tablefmt="plain",
        floatfmt=".3f",
        missingval=""
    ))
    
    #since i used seaborn to plot, converting the stored values in arrays to a pandas df was useful for plotting
    
    df = pd.DataFrame(results, columns=["time", "x", "y", "ax", "ay", "vx", "vy","r"])
    
    ax=sns.lineplot(data=df, x="x", y="y", marker=False, legend=0,sort=False)
    sns.scatterplot(x=[0], y=[0], color="yellow", s=150, edgecolor="black", ax=ax)
    
    
    
    
    ax.set(xlabel='x position ', ylabel= 'y position')
    ax.set(title=f'(dt: {dt}, steps: {steps}, x0: {x0}, y0: {y0}, vx0: {vx0}, vy0: {vy0})')
    plt.suptitle('Motion of a planet around Sun')
    plt.figure(dpi=300)
    plt.show()
    
    return 


#1st iteration - a coarse timestep
dt=0.001
steps=525
x0=0.300
y0=0.000
vx0=0.0
vy0=0.5

leapfrog(dt,steps,x0,y0,vx0,vy0)


#2nd iteration
dt=0.01
steps=2000
x0=0.500
y0=0.000
vx0=0
vy0=1.630
leapfrog(dt,steps,x0,y0,vx0,vy0)

#3rd iteration
dt=0.01
steps=2000
x0=0.500
y0=0.00
vx0=0
vy0=0.630
leapfrog(dt,steps,x0,y0,vx0,vy0)






























