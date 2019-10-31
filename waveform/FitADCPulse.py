import ROOT
import numpy as np
import sys
import sys, os
from utils import *
from glob import glob

straw_lengths = {-3: 120., 3: 120.,  36: 100., 90: 48}
distances = {-3: np.array([10,35,60,85,110]), 3: np.array([10,35,60,85,110]), 36: np.array([9,29.5,50,70.5,91]), 90: np.array([9,24,39])}

straws = straw_lengths.keys()
if len(sys.argv) > 1:
  straws = [int(sys.argv[1])]

attenuation_length = 140.

results = {straw: {'cal': [], 'hv': []} for straw in straw_lengths}
errors = {straw: {'cal': [], 'hv': []} for straw in straw_lengths}
fin = open("FitCurrentPulse.dat")
for line in fin:
  data = line.split()
  straw = int(data[0])
  position = int(data[1])
  cal_or_hv = ['cal','hv'][int(data[2])]
  keys = ["meant","sigma","t0","pole0","pole1","zero","scale","time1","scale1","time2","scale2"]
  vals = {}
  errs = {}
  for j in range(len(keys)):
    vals[keys[j]] = float(data[3+j*2])
    errs[keys[j]] = float(data[3+j*2+1])
  results[straw][cal_or_hv].append(vals)
  errors[straw][cal_or_hv].append(errs)

for straw in straws:
  results[straw]['hv'].reverse()
  errors[straw]['hv'].reverse()



fs = 5000 #MHz (to match 5 gsps scope data)

hnllcal = ROOT.TH1F("","",750,0,150)
hnllhv = ROOT.TH1F("","",750,0,150)
hnll = ROOT.TH1F("","",750,0,150)
hdata = ROOT.TH1F("","",750,0,150)

def fillhisto(h,meant_override,meant,sigma,t0,pole0,pole1,zero,scale1,scale2,time1,time2,scale):
  y_close = np.zeros(hnll.GetNbinsX())
  for i in range(hnll.GetNbinsX()):
    t_gaus = hnll.GetBinCenter(i+1)
    val_gaus_close = 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(t_gaus-meant_override)**2/(2*sigma**2))
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

  for i in range(h.GetNbinsX()):
    h.SetBinContent(i+1,output1[i])

whichdist = 0

def nll(meant_cal, meant_hv, scale, pole):
  caldist = whichdist
  hvdist = 4-whichdist
  fillhisto(hnllcal,meant_cal,**results[straw]['cal'][caldist])
  fillhisto(hnllhv,meant_hv,**results[straw]['hv'][hvdist])

  output2 = [(hnllcal.GetBinContent(i+1)+hnllhv.GetBinContent(i+1))*scale for i in range(hnll.GetNbinsX())]
  dc2 = digital_coefficients([pole],[],fs)
  output3 = apply_filter_from_dc(dc2,output2)
  output3 /= np.max(output3)

  tot = 0
  for i in range(hnll.GetNbinsX()):
    hnll.SetBinContent(i+1,output3[i])
    tot += ((hnll.GetBinContent(i+1)-hdata.GetBinContent(i+1))/hdata.GetBinError(i+1))**2
  #print meant,sigma,t0,pole0,pole1,zero,scale,":",tot
  return tot

if __name__ == "__main__":
  from iminuit import Minuit

  straws = [-3]
  num_positions = {-3: 5, 3: 5, 36: 5, 90: 3}

  fout = open("FitADCPulse.dat","w")

  for straw in straws:
    for position in range(1,num_positions[straw]+1):
        whichdist = position-1
        print "Straw %d, position %d" % (straw,position),
    
        data = np.genfromtxt("/data/HD2/rbonventre/prototypeData/analogDRAC/wvf_ch%d_pos%d.dat" % (straw,position))
        data[:,3] -= np.mean(data[:250,3])
        data[:,6] /= np.max(data[:,3])
        data[:,3] /= np.max(data[:,3])
        for i in range(hdata.GetNbinsX()):
          hdata.SetBinContent(i+1,data[250+i,3])
          hdata.SetBinError(i+1,data[250+i,6])
    
        m0 = Minuit(nll,print_level=0,pedantic=False,
            meant_cal=35, meant_hv=35,
            scale=1, fix_scale=True,
            pole=1.6)
        m0.migrad()
        #nll(**m0.values)
  
        print m0.values
        print m0.errors
        keys = ["meant_cal","meant_hv","pole"]
        fout.write("%d %d " % (straw,position))
        for key in keys:
          fout.write("%.8f %.8f " % (m0.values[key],m0.errors[key]))
        fout.write("\n")
  
        #hdata.Draw()
        #hdata.SetLineColor(ROOT.kRed)
        #hnll.Draw("same")
        #raw_input()
  fout.close() 
