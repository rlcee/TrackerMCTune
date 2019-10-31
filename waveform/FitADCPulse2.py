import ROOT
import numpy as np
import sys
import sys, os
from utils import *
from glob import glob

fs = 2500 #MHz (to match 5 gsps scope data)

# set this to the model parameters for a 48 cm straw in the center
reflection_time = 5.08 + (2*48.5-2*24)/51.4
reflection_scale = 0.58 * np.exp(-(2*48.5-2*24)/256.8)
t0 = 0.8 + 0.015*24
parameters = {"meant": 30, "sigma": 2.42, "t0": t0, "pole0": 160, "pole1": 6, "zero": 0.2, "scale1": reflection_scale, "scale2": 0, "time1": reflection_time, "time2": 0, "scale": 1}


hnll = ROOT.TH1F("","",800,0,320)
hnll2 = ROOT.TH1F("","",1600,0,640)
hdata = ROOT.TH1F("","",16,0,320)

def fillhisto(meant_override,meant,sigma,t0,pole0,pole1,zero,scale1,scale2,time1,time2,scale):
  y_close = np.zeros(hnll2.GetNbinsX())
  for i in range(hnll2.GetNbinsX()):
    t_gaus = hnll2.GetBinCenter(i+1)
    val_gaus_close = 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(t_gaus-meant_override)**2/(2*sigma**2))
    for j in range(hnll2.GetNbinsX()-i):
      t_tail = hnll2.GetBinCenter(j+1)
      val_close = val_gaus_close / (t_tail + t0) + 1e-5
      y_close[i+j] += val_close
   
  dc1 = digital_coefficients([pole0,pole1],[zero],fs)
  output1 = apply_filter_from_dc(dc1,y_close)
  reflection1 = np.interp(np.linspace(0-time1,640-time1,len(output1)),np.linspace(0,640,len(output1)),output1,left=0)
  reflection2 = np.interp(np.linspace(0-time2,640-time2,len(output1)),np.linspace(0,640,len(output1)),output1,left=0)
  output1 += reflection1*scale1
  output1 += reflection2*scale2

  output1 *= 10**8
  output1 *= scale
  reflection1 *= scale1*scale*10**8

  return output1

output2_x = np.array([hnll2.GetBinCenter(i+1) for i in range(hnll2.GetNbinsX())])
output2 = fillhisto(300,**parameters)

def nll(meant, pole, scale):
  dc2 = digital_coefficients([pole],[],fs)
  output3 = apply_filter_from_dc(dc2,output2)
  output3 /= np.max(output3)

  output3 *= scale

  offset = (300 - meant)
  

  tot = 0
  for i in range(hdata.GetNbinsX()):
    sim_val = np.interp(hdata.GetBinLowEdge(i+1) + offset,output2_x,output3)
    tot += ((sim_val-hdata.GetBinContent(i+1)))**2*100
  return tot

def nll_draw(meant, pole, scale):
  dc2 = digital_coefficients([pole],[],fs)
  output3 = apply_filter_from_dc(dc2,output2)
  output3 /= np.max(output3)

  output3 *= scale

  offset = int((300 - meant)/hnll.GetBinWidth(1))

  for i in range(hnll.GetNbinsX()):
    hnll.SetBinContent(i+1,output3[i+offset])

if __name__ == "__main__":
  from iminuit import Minuit

  fout = open("FitADCPulse2.dat","w")
  ROOT.gSystem.Load("../event/Dict.so")
  f = ROOT.TFile("/data/HD2/rbonventre/prototypeData/lblPrototype/run_2683.root")
  t = f.Get("T")
  idone = 0
  for i in range(t.GetEntries()):
    t.GetEntry(i)
    for j in range(t.events.straws.size()):
      if t.events.straws[j].channel != 5:
        continue
      data = [t.events.straws[j].samples[k] for k in range(16)]
      data -= np.mean(data[:3])
      peak = np.max(data)
      if np.max(data) > 800 or np.max(data) < 500:
        continue
      data /= np.max(data)
      for k in range(hdata.GetNbinsX()):
        hdata.SetBinContent(k+1,data[k])
    
      m0 = Minuit(nll,print_level=0,pedantic=False,
            meant=80, pole=3, scale=1)
#      nll(**m0.values)
#      nll_draw(**m0.values)
#      x = np.array([hdata.GetBinLowEdge(i+1) for i in range(hdata.GetNbinsX())])
#      y = np.array([hdata.GetBinContent(i+1) for i in range(hdata.GetNbinsX())])
#      g = ROOT.TGraph(hdata.GetNbinsX(),x,y)
#      g.SetMarkerStyle(21)
#      g.Draw("ALP")
      #hdata.Draw()
      #hdata.SetLineColor(ROOT.kRed)
#      hnll.Draw("same")
#      raw_input()
      m0.migrad()
        #nll(**m0.values)
  
#      print m0.values
#      print m0.errors
      if not m0.migrad_ok() or m0.fval > 50:
        print i,j,"BAD"
        continue
      idone += 1
      print i,j,"SUCCESS",idone
      keys = ["meant","pole","scale"]
      fout.write("%d %d %.2f %.2f " % (i,j,m0.fval,peak))
      for key in keys:
        fout.write("%.8f %.8f " % (m0.values[key],m0.errors[key]))
      fout.write("\n")
  
#      if m0.values['pole'] > 3.5 or m0.values['pole'] < 3:
#      if True:
#        print m0.values['pole'],m0.values['scale'],peak
#        nll_draw(**m0.values)
#        x = np.array([hdata.GetBinLowEdge(i+1) for i in range(hdata.GetNbinsX())])
#        y = np.array([hdata.GetBinContent(i+1) for i in range(hdata.GetNbinsX())])
#        g = ROOT.TGraph(hdata.GetNbinsX(),x,y)
#        g.SetMarkerStyle(21)
#        g.Draw("ALP")
##        hdata.Draw()
##        hdata.SetLineColor(ROOT.kRed)
#        hnll.Draw("same")
#        raw_input()
