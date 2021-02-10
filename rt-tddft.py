#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 15:06:55 2021

@author: yavar001
"""
import re
import numpy as np

from scipy.fft import fft, fftfreq

import math
import matplotlib.pyplot as plt





rotated=np.empty((0,2))
#======================
def rotate_spectrum (data):
    rot_fft=rotated
    rownum=data.size
    for i in range(rownum):
        rowdata=data[i]
        re = rowdata.real
        im = rowdata.imag
        ab= abs(complex(re,im))
        r = math.sqrt (re**2 + im**2)
        if abs (r - ab) > 1e-6:
            raise Exception ("abs not equal to sqrt(re^2 + im^2)")

        angle = abs (math.atan2 (im, re))

        if angle > math.pi:
            raise Exception ("atan2 out of range")

        re_out = ab * math.cos (angle)
        im_out = ab * math.sin (angle)
        rot_fft=np.append(rot_fft, np.array([[re_out,im_out]]), axis=0)

    return rot_fft 





def rttddft_data(fname):
    
    
    list1=[]
    
    ff=open(fname, "r")
    for line in ff:
        if re.search("# Dipole moment", line):
            list1.append(line)
    ff.close()        
    
    

    
    xt=np.array([float(x.split()[1]) for x in list1 if "kick_x" in x.split()])
    xf=np.array([float(x.split()[2]) for x in list1 if "kick_x" in x.split()])
    yt=np.array([float(x.split()[1]) for x in list1 if "kick_y" in x.split()])
    yf=np.array([float(x.split()[3]) for x in list1 if "kick_y" in x.split()])
    zt=np.array([float(x.split()[1]) for x in list1 if "kick_z" in x.split()])
    zf=np.array([float(x.split()[4]) for x in list1 if "kick_z" in x.split()])
    
    
    
    
    #########    SHIFT FIRST POINT TO ZERO


    xf=xf-xf[0]
    yf=yf-yf[0]
    zf=zf-zf[0]
    



    ##############  damping


    x=xt
    xend=x[-1]



    ##  to control damp

    tau=-xend/(np.log(0.02))


    damp=np.exp(-x/tau)


    xdatadamp=xf*damp
    ydatadamp=yf*damp
    zdatadamp=zf*damp


    #######  PADDING




    #######################################################

    dx=x[1]-x[0]
    npadd=50000





    ##########  PADDING  #####################################
    #    paddend=npadd*dx


    paddata=np.zeros((npadd,))


    padxdata=np.hstack((xdatadamp, paddata))
    padydata=np.hstack((ydatadamp, paddata))
    padzdata=np.hstack((zdatadamp, paddata))

    #########################################################



    xfft=fft(padxdata)
    yfft=fft(padydata)
    zfft=fft(padzdata)


    xrotfft=rotate_spectrum(xfft)
    yrotfft=rotate_spectrum(yfft)
    zrotfft=rotate_spectrum(zfft)



    xwdata=np.loadtxt("xw.dat", skiprows=1)


    n=xfft.size
    period = (n-1)*dx


    dw=2*np.pi/period
    m=n//2

    wmin=0.0
    wmax=m*dw

    xw=np.linspace(wmin, wmax, m)

    ss=xw*(xrotfft[:m,1]+yrotfft[:m,1]+zrotfft[:m,1])
    
    rtdata=np.vstack((xw,ss)).T
    basename=fname.split(".")[0]
    
    np.savetxt(basename+".rtspdata", rtdata, fmt='%8.4f')
    return rtdata




def rttddft_plotspectrum(filename, win0, win1):
    data1=rttddft_data(filename)
    basename=filename.split(".")[0]
    xdat=27.21*data1[:,0]
    ydat=data1[:,1]
    yrange=ydat[ (xdat > win0) & (xdat < win1) ]
    ymax=np.max(yrange)*1.1
    plt.ylim(0, ymax)
    plt.xlim(win0, win1)
    plt.xlabel("Energy (eV)", fontsize=16)
    plt.ylabel("Absorption (arb. unit.)", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(basename)
    plt.plot(27.21*data1[:,0], data1[:,1])
    plt.tight_layout()
    plt.savefig(basename+"_rttd.png",dpi=300, bbox_inches = "tight") 
    plt.show()    
    return


import argparse
parser = argparse.ArgumentParser()
parser.add_argument("nwoutput", help="name of your nwchem output file")
parser.add_argument('-wi', action='store', dest='windstart', type=float, default=0.0,
                    help='starting point of energy window')
parser.add_argument('-wf', action='store', dest='windend', type=float, default=25.0,
                    help='end point for energy window')

args = parser.parse_args()
rttddft_plotspectrum(args.nwoutput, args.windstart, args.windend)
    
