import ROOT
import numpy as np
import sys
import sys, os
from utils import *
from glob import glob

fs = 5000 #MHz (to match 5 gsps scope data)

hnll = ROOT.TH1F("","",750,0,150)
hdata = ROOT.TH1F("","",750,0,150)

def nll(meant,sigma,t0,pole0,pole1,zero,scale1,scale2,time1,time2,scale):
  y_close = np.zeros(hnll.GetNbinsX())
  for i in range(hnll.GetNbinsX()):
    t_gaus = hnll.GetBinCenter(i+1)
    val_gaus_close = 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(t_gaus-meant)**2/(2*sigma**2))
    for j in range(hnll.GetNbinsX()-i):
      t_tail = hnll.GetBinCenter(j+1)
      val_close = val_gaus_close / (t_tail + t0) + 1e-5
      y_close[i+j] += val_close
   
  dc1 = digital_coefficients([pole0,pole1,500],[zero],fs)
  output1 = apply_filter_from_dc(dc1,y_close)
  reflection1 = np.interp(np.linspace(0-time1,150-time1,len(output1)),np.linspace(0,150,len(output1)),output1,left=0)
  reflection2 = np.interp(np.linspace(0-time2,150-time2,len(output1)),np.linspace(0,150,len(output1)),output1,left=0)
  output1 += reflection1*scale1
  output1 += reflection2*scale2

  output1 *= 10**8
  output1 *= scale

  for i in range(hnll.GetNbinsX()):
    hnll.SetBinContent(i+1,output1[i])

  tot = 0
  for i in range(hnll.GetNbinsX()):
    tot += ((hnll.GetBinContent(i+1)-hdata.GetBinContent(i+1))/hdata.GetBinError(i+1))**2
  #print meant,sigma,t0,pole0,pole1,zero,scale,":",tot
  return tot

if __name__ == "__main__":
  from iminuit import Minuit

  straws = [-3,3,36,90]
  num_positions = {-3: 5, 3: 5, 36: 5, 90: 3}

  fout = open("FitCurrentPulse.dat","w")

  for straw in straws:
    for position in range(1,num_positions[straw]+1):
      for cal_or_hv in range(2): # cal = 0, hv = 1

        print "Straw %d, position %d" % (straw,position),
        if cal_or_hv == 0:
          print "cal"
        else:
          print "hv"
    
        data = np.genfromtxt("/data/HD2/rbonventre/prototypeData/analogDRAC/wvf_ch%d_pos%d.dat" % (straw,position))
        data[:,1] -= np.mean(data[:250,1])
        data[:,2] -= np.mean(data[:250,2])
        for i in range(hdata.GetNbinsX()):
          if cal_or_hv == 0:
            hdata.SetBinContent(i+1,data[250+i,1]) # cal side
            hdata.SetBinError(i+1,data[250+i,4])
          else:
            hdata.SetBinContent(i+1,data[250+i,2]) # hv side
            hdata.SetBinError(i+1,data[250+i,5])
    
        m0 = Minuit(nll,print_level=0,pedantic=False,
            meant=32,           limit_meant=(20,50),
            sigma=2.5,          fix_sigma=False,
            t0=3,               fix_t0=False,
            pole0=160,          fix_pole0=True,
            pole1=6,            fix_pole1=True,
            zero=0.2,           fix_zero=True,
            scale1=0.2,         limit_scale1=(0,1),
            scale2=0.1,         limit_scale2=(0,1),
            time1=8.5,          limit_time1=(2,12),
            #time1=5,           limit_time1=(2,8),
            time2=17,           limit_time2=(12,20),
            #time2=10,          limit_time2=(8,20),
            scale=1)
        m0.migrad()
        #nll(**m0.values)
  
        print m0.values
        print m0.errors
        keys = ["meant","sigma","t0","pole0","pole1","zero","scale","time1","scale1","time2","scale2"]
        fout.write("%d %d %d " % (straw,position,cal_or_hv))
        for key in keys:
          fout.write("%.8f %.8f " % (m0.values[key],m0.errors[key]))
        fout.write("\n")
  
        #hdata.Draw()
        #hdata.SetLineColor(ROOT.kRed)
        #hnll.Draw("same")
        #raw_input()
  fout.close() 
