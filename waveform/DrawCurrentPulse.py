import numpy as np

straw_lengths = {-3: 120., 3: 120.,  36: 100., 90: 48}
distances = {-3: np.array([10,35,60,85,110]), 3: np.array([10,35,60,85,110]), 36: np.array([9,29.5,50,70.5,91]), 90: np.array([9,24,39])}

if __name__ == "__main__":
  import sys
  import os
  import numpy as np
  import ROOT
  from array import array
  sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'src'))
  from utils import *

  straw = int(sys.argv[1])
  cal_or_hv = int(sys.argv[2]) # hv = 0, cal = 1
  cal_hv_string = ['cal','hv']
  pos = int(sys.argv[3])

  fs = 5000 #MHz (to match 5 gsps scope data)
  
  hbf = ROOT.TH1F("hbf_%d" % pos,"hbf",750,0,150)
  hdata = ROOT.TH1F("hdata_%d" % pos,"Straw %d %s: %d cm" % (straw,cal_hv_string[cal_or_hv],distances[straw][pos-1]),750,0,150)
  hparam = ROOT.TH1F("hparam_%d" % pos,"hparam",750,0,150)
  
  def nll(h,meant,sigma,t0,pole0,pole1,zero,scale1,scale2,time1,time2,scale):
    y_close = np.zeros(hdata.GetNbinsX())
    for i in range(hdata.GetNbinsX()):
      t_gaus = hdata.GetBinCenter(i+1)
      val_gaus_close = 1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(t_gaus-meant)**2/(2*sigma**2))
      for j in range(hdata.GetNbinsX()-i):
        t_tail = hdata.GetBinCenter(j+1)
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
  
    for i in range(hdata.GetNbinsX()):
      h.SetBinContent(i+1,output1[i])
  
   
  data = np.genfromtxt("/data/HD2/rbonventre/prototypeData/analogDRAC/wvf_ch%d_pos%d.dat" % (straw,pos))
  data[:,1] -= np.mean(data[:250,1])
  data[:,2] -= np.mean(data[:250,2])
  for i in range(hdata.GetNbinsX()):
    if cal_or_hv == 0:
      hdata.SetBinContent(i+1,data[250+i,1]) # cal side
    else:
      hdata.SetBinContent(i+1,data[250+i,2]) # hv side
  
  results = np.genfromtxt("FitCurrentPulse.dat")
  results_dict = {}
  keys = ["meant","sigma","t0","pole0","pole1","zero","scale","time1","scale1","time2","scale2"]
  for line in results:
    if int(line[0]) == straw and int(line[1]) == pos and int(line[2]) == cal_or_hv:
      for i in range(len(keys)):
        results_dict[keys[i]] = line[3+i*2]
  nll(hbf,**results_dict)
  params = results_dict.copy()
#  params['sigma0'] = 2.31
#   params['sigma0'] = 2.162879 + 0.000158*(straw_lengths[straw]/2.-distances[straw][results_position])**2
#  params['tau0'] = 1.06 + 0.032194*distances[straw][results_position]
#  params['time1'] = 4.62 + 0.022324*(straw_lengths[straw] + (straw_lengths[straw]-distances[straw][results_position])-distances[straw][results_position])
#  params['time2'] = 2*4.62 + 2*straw_lengths[straw]/35.
#  params['scale1'] = 0.61*np.exp(-(straw_lengths[straw] + (straw_lengths[straw]-distances[straw][results_position])-distances[straw][results_position])/140.)
#  params['scale2'] = 0.61*0.61*np.exp(-2*straw_lengths[straw]/160.)
#  params['sigma0'] = 2.75
#  params['sigma0'] = 2.03
#  params['tau0'] = 0.83
#  params['time1'] = 8.41
#  params['time2'] = 16.86
#  params['scale1'] = 0.31
#  params['scale2'] = 0.092

  if cal_or_hv == 0:
    results_position = pos-1
  else:
    results_position = 5-pos
  params['sigma'] = 2.42
  params['t0'] = 0.8 + 0.01525*distances[straw][results_position]
  params['time1'] = 5.08 + 0.019463*(straw_lengths[straw] + (straw_lengths[straw]-distances[straw][results_position])-distances[straw][results_position])
  params['time2'] = 2*5.08 + 0.019463*2*straw_lengths[straw]
  params['scale1'] = 0.58*np.exp(-(straw_lengths[straw] + (straw_lengths[straw]-distances[straw][results_position])-distances[straw][results_position])/256.8)
  params['scale2'] = 0.58**2*np.exp(-2*straw_lengths[straw]/256.8)
#  params['sigma0'] = 2.75
#  params['sigma0'] = 2.03
#  params['tau0'] = 0.83
#  params['time1'] = 8.41
#  params['time2'] = 16.86
#  params['scale1'] = 0.31
#  params['scale2'] = 0.092

  nll(hparam,**params)

  deltat = hbf.GetBinCenter(hbf.GetMaximumBin())-hparam.GetBinCenter(hparam.GetMaximumBin())
  params['meant'] += deltat
  nll(hparam,**params)

  hparam.Scale(hbf.GetMaximum()/hparam.GetMaximum())
 

  print "sigma: \t%7.2f\t%7.2f" % (results_dict["sigma"],params["sigma"])
  print "t0:    \t%7.2f\t%7.2f" % (results_dict["t0"],params["t0"])
  print "time1: \t%7.2f\t%7.2f" % (results_dict["time1"],params["time1"])
  print "scale1:\t%7.2f\t%7.2f" % (results_dict["scale1"],params["scale1"])
  print "time2: \t%7.2f\t%7.2f" % (results_dict["time2"],params["time2"])
  print "scale2:\t%7.3f\t%7.2f" % (results_dict["scale2"],params["scale2"])

  hdata.SetStats(0)
  hdata.GetXaxis().SetTitle("Time (ns)")
  hdata.GetYaxis().SetTitle("Pulse size (mV)")
  hdata.SetLineColor(ROOT.kBlack)
  hdata.SetLineWidth(2)
  hdata.SetLineStyle(2)
  hdata.Draw("HIST")
  hbf.SetLineColor(ROOT.kRed)
  hbf.Draw("same HIST")
  hparam.SetLineColor(ROOT.kBlue)
  hparam.Draw("same HIST")
  l = ROOT.TLegend(0.45,0.65,0.85,0.85)
  ###hdata.GetXaxis().SetRangeUser(20,50)
  l.AddEntry(hdata,"Data","l")
  l.AddEntry(hbf,"Best fit","l")
  l.AddEntry(hparam,"Param","l")
  l.Draw()
  raw_input()
