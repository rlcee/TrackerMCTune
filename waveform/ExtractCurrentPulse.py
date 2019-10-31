import ROOT
from array import array
from utils import *
import sys
import os

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
hnll = ROOT.TH1F("","",750,0,150)
def histo(hnll,meant,sigma,t0,pole0,pole1,zero,scale1,scale2,time1,time2,scale):
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


# FIT PROP VELOCITY
def fit_prop_velocity():
  mg = ROOT.TMultiGraph("mg","mg")
  l = ROOT.TLegend(0.15,0.65,0.35,0.85)
  color_index = 1
  for straw in straws:
    y = array('d',2*distances[straw]-straw_lengths[straw])
    x = array('d',[results[straw]['cal'][i]['meant']-results[straw]['hv'][len(distances[straw])-1-i]['meant'] for i in range(len(distances[straw]))])
    mean_val = x[len(distances[straw])/2]
    for i in range(len(x)):
      x[i] -= mean_val
    g = ROOT.TGraph(len(distances[straw]),x,y)

    g.SetMarkerStyle(21)
    g.SetMarkerColor(color_index)
    g.SetLineColor(color_index)
    color_index += 1
    l.AddEntry(g,"Straw %d" % straw,"p")
    
    mg.Add(g)
  
  fr = mg.Fit("pol1","SQ")
  mg.SetTitle("Effective Propagation Velocity")
  mg.Draw("ALP")
  l.Draw()
  mg.GetXaxis().SetTitle("#Delta t (ns)")
  mg.GetYaxis().SetTitle("#Delta d (cm)")
  params = fr.GetParams()
  paramerrs = fr.GetErrors()
  print "Prop velocity: ",params[1],"+-",paramerrs[1]
  p = ROOT.TPaveText(0.5,0.2,0.7,0.3,"NDC NB")
  p.AddText("v = %.1f #pm %.1f (cm/ns)" % (params[1],paramerrs[1]))
  p.SetFillStyle(0)
  p.SetBorderSize(0)
  p.Draw()
  raw_input()



def calc_sigma():
  for straw in straws:
    # CALCULATE PULSE SHAPE SIGMA
    print straw,"Mean sigma: ",np.mean(np.append([results[straw]['hv'][i]['sigma'] for i in range(len(distances[straw]))],[results[straw]['hv'][i]['sigma'] for i in range(len(distances[straw]))])),"+-",np.std(np.append([results[straw]['hv'][i]['sigma'] for i in range(len(distances[straw]))],[results[straw]['hv'][i]['sigma'] for i in range(len(distances[straw]))]))

  mg = ROOT.TMultiGraph("mg","mg")
  l = ROOT.TLegend(0.15,0.65,0.35,0.85)
  color_index = 1
  sigmas = []
  for straw in straws:
    x = array('d',straw_lengths[straw]/2.-distances[straw])
    y = array('d',[results[straw]['cal'][position]['sigma'] for position in range(len(distances[straw]))])
    xerr = array('d',[0 for i in range(len(distances[straw]))])
    yerr = array('d',[errors[straw]['cal'][position]['sigma'] for position in range(len(distances[straw]))])
    gc = ROOT.TGraphErrors(len(distances[straw]),x,y,xerr,yerr)
    x = array('d',straw_lengths[straw]/2.-distances[straw])
    y = array('d',[results[straw]['hv'][position]['sigma'] for position in range(len(distances[straw]))])
    xerr = array('d',[0 for i in range(len(distances[straw]))])
    yerr = array('d',[errors[straw]['hv'][position]['sigma'] for position in range(len(distances[straw]))])
    gh = ROOT.TGraphErrors(len(distances[straw]),x,y,xerr,yerr)
    gc.SetMarkerStyle(21)
    gh.SetMarkerStyle(22)
    gc.SetMarkerColor(color_index)
    gc.SetLineColor(color_index)
    gh.SetMarkerColor(color_index)
    gh.SetLineColor(color_index)
    color_index += 1
    l.AddEntry(gc,"Cal %d" % straw,"p")
    l.AddEntry(gh,"HV %d" % straw,"p")
    
    mg.Add(gc)
    mg.Add(gh)

    for hvcal in ['hv','cal']:
      for position in range(len(distances[straw])):
        sigmas.append(results[straw][hvcal][position]['sigma'])
  print "Overall mean sigma: ",np.mean(sigmas),"+-",np.std(sigmas)

  f1 = ROOT.TF1("f1","[0] + [1]*x*x",-120,120)
  fr = mg.Fit("f1","SQ")
  mg.SetTitle("Pulse shape sigma")
  mg.Draw("ALP")
  l.Draw()
  mg.GetXaxis().SetTitle("Distance from straw end (cm)")
  mg.GetYaxis().SetTitle("sigma (ns)")
  
  params = fr.GetParams()
  paramerrs = fr.GetErrors()
  print "sigma [ns] = (%f +- %f) + (%f +- %f) * dist + (%f +- %f) * dist^2" % (params[0],paramerrs[0],0,0,params[1],paramerrs[1])
  p = ROOT.TPaveText(0.13,0.55,0.63,0.61,"NDC NB")
  p.AddText("t_{0} (ns) = (%.2f #pm %.2f) + (%.3f #pm %.3f) *  dist[cm]" % (params[0],paramerrs[0],params[1],paramerrs[1]))
  p.SetFillStyle(0)
  p.SetBorderSize(0)
  p.Draw()

  raw_input()



def fit_t0():
  mg = ROOT.TMultiGraph("mg","mg")
  l = ROOT.TLegend(0.15,0.65,0.35,0.85)
  color_index = 1
  for straw in straws:
    # FIT PULSE SHAPE T0
    x = array('d',distances[straw])
    y = array('d',[results[straw]['cal'][i]['t0'] for i in range(len(distances[straw]))])
    xerr = array('d',[0 for i in range(len(distances[straw]))])
    yerr = array('d',[errors[straw]['cal'][i]['t0'] for i in range(len(distances[straw]))])
    gc = ROOT.TGraphErrors(len(distances[straw]),x,y,xerr,yerr)
    x = array('d',distances[straw])
    y = array('d',[results[straw]['hv'][i]['t0'] for i in range(len(distances[straw]))])
    xerr = array('d',[0 for i in range(len(distances[straw]))])
    yerr = array('d',[errors[straw]['hv'][i]['t0'] for i in range(len(distances[straw]))])
    gh = ROOT.TGraphErrors(len(distances[straw]),x,y,xerr,yerr)
    gc.SetMarkerStyle(21)
    gh.SetMarkerStyle(22)
    gc.SetMarkerColor(color_index)
    gc.SetLineColor(color_index)
    gh.SetMarkerColor(color_index)
    gh.SetLineColor(color_index)
    color_index += 1
    l.AddEntry(gc,"Cal %d" % straw,"p")
    l.AddEntry(gh,"HV %d" % straw,"p")
    
    mg.Add(gc)
    mg.Add(gh)
  fr = mg.Fit("pol1","SQ")
  mg.SetTitle("Pulse shape t0")
  mg.Draw("ALP")
  l.Draw()
  mg.GetXaxis().SetTitle("Distance from straw end (cm)")
  mg.GetYaxis().SetTitle("t0 (ns)")
  
  params = fr.GetParams()
  paramerrs = fr.GetErrors()
  print "t0 [ns] = (%f +- %f) + (%f +- %f) * dist" % (params[0],paramerrs[0],params[1],paramerrs[1])
  p = ROOT.TPaveText(0.13,0.55,0.63,0.61,"NDC NB")
  p.AddText("t_{0} (ns) = (%.2f #pm %.2f) + (%.3f #pm %.3f) *  dist[cm]" % (params[0],paramerrs[0],params[1],paramerrs[1]))
  p.SetFillStyle(0)
  p.SetBorderSize(0)
  p.Draw()

  raw_input()



def fit_scaling():
  # FIT SCALING
  x = array('d',distances[straw])
  y = array('d',[results[straw]['cal'][i]['scale'] for i in range(len(distances[straw]))])
  gc = ROOT.TGraph(len(distances[straw]),x,y)
  x = array('d',distances[straw])
  y = array('d',[results[straw]['hv'][i]['scale'] for i in range(len(distances[straw]))])
  gh = ROOT.TGraph(len(distances[straw]),x,y)
  
  mg = ROOT.TMultiGraph("mg","mg")
  mg.Add(gc)
  mg.Add(gh)
  fr = mg.Fit("pol1","SQ")
  mg.Draw("ALP")
  mg.GetXaxis().SetTitle("Distance from straw end (cm)")
  mg.GetYaxis().SetTitle("Scaling factor")
  
  params = fr.GetParams()
  paramerrs = fr.GetErrors()
  print "scaling factor = (%f +- %f) + (%f +- %f) * dist" % (params[0],paramerrs[0],params[1],paramerrs[1])
  raw_input()







def fit_attenuation_length():
  # FIT ATTENUATION LENGTH
  mg = ROOT.TMultiGraph("mg","mg")
  l = ROOT.TLegend(0.65,0.65,0.85,0.85)
  color_index = 1
  import copy
  for straw in straws:
    x = array('d',distances[straw])
    y = array('d') 
    for i in range(len(distances[straw])):
      temp = copy.copy(results[straw]['cal'][i])
      temp['scale1'] = 0
      temp['scale2'] = 0
      histo(hnll,**temp)
      y.append(hnll.GetMaximum())
    gc = ROOT.TGraph(len(distances[straw]),x,y)
    x = array('d',distances[straw])
    y = array('d') 
    for i in range(len(distances[straw])):
      temp = copy.copy(results[straw]['hv'][i])
      temp['scale1'] = 0
      temp['scale2'] = 0
      histo(hnll,**temp)
      y.append(hnll.GetMaximum())
    gh = ROOT.TGraph(len(distances[straw]),x,y)
    gc.SetMarkerStyle(21)
    gh.SetMarkerStyle(22)
    gc.SetMarkerColor(color_index)
    gc.SetLineColor(color_index)
    gh.SetMarkerColor(color_index)
    gh.SetLineColor(color_index)
    color_index += 1
    l.AddEntry(gc,"Cal %d" % straw,"p")
    l.AddEntry(gh,"HV %d" % straw,"p")
    
    mg.Add(gc)
    mg.Add(gh)
  f0 = ROOT.TF1("f0","[0]*exp(-x/[1])")
  f0.SetParameter(0,y[0])
  f0.SetParameter(1,240.)
  fr = mg.Fit("f0","SQ")
  mg.SetTitle("Attenuation Length")
  mg.Draw("ALP")
  l.Draw()
  mg.GetXaxis().SetTitle("Distance from straw end (cm)")
  mg.GetYaxis().SetTitle("Peak minus pedestal without reflections (arb.)")
  
  params = fr.GetParams()
  paramerrs = fr.GetErrors()
  print "Atten. Length [cm] = (%f +- %f)" % (params[1],paramerrs[1])
  p = ROOT.TPaveText(0.17,0.24,0.5,0.3,"NDC NB")
  p.AddText("L = %.1f #pm %.1f (cm)" % (params[1],paramerrs[1]))
  p.SetFillStyle(0)
  p.SetBorderSize(0)
  p.Draw()
  raw_input()


def fit_reflection_time():
  # REFLECTION TIME SHIFT AND VELOCITY
  mg = ROOT.TMultiGraph("mg","mg")
  l = ROOT.TLegend(0.65,0.65,0.85,0.85)
  color_index = 1
  for straw in straws:
    x = array('d',[(straw_lengths[straw] + (straw_lengths[straw]-distances[straw][i])) - distances[straw][i] for i in range(len(distances[straw]))])
    y = array('d',[results[straw]['hv'][i]['time1'] for i in range(len(distances[straw]))])
#    if straw == 3:
#      y[0] += 8
    gh = ROOT.TGraph(len(distances[straw]),x,y)
    x = array('d',[(straw_lengths[straw] + (straw_lengths[straw]-distances[straw][i])) - distances[straw][i] for i in range(len(distances[straw]))])
    y = array('d',[results[straw]['cal'][i]['time1'] for i in range(len(distances[straw]))])
    gc = ROOT.TGraph(len(distances[straw]),x,y)
    gc.SetMarkerStyle(21)
    gh.SetMarkerStyle(22)
    gc.SetMarkerColor(color_index)
    gc.SetLineColor(color_index)
    gh.SetMarkerColor(color_index)
    gh.SetLineColor(color_index)
    color_index += 1
    l.AddEntry(gc,"Cal %d" % straw,"p")
    l.AddEntry(gh,"HV %d" % straw,"p")
  
    mg.Add(gc)
    mg.Add(gh)
  fr = mg.Fit("pol1","S")
  mg.SetTitle("Time delay of first reflection")
  mg.Draw("ALP")
  l.Draw()
  mg.GetXaxis().SetTitle("Distance travelled by reflection (cm)")
  mg.GetYaxis().SetTitle("#delta t_{reflection} (ns)")
  
  params = fr.GetParams()
  paramerrs = fr.GetErrors()
  print "reflection t [ns] = (%f +- %f) + (%f +- %f) * dist" % (params[0],paramerrs[0],params[1],paramerrs[1])
  p = ROOT.TPaveText(0.4,0.18,0.85,0.25,"NDC NB")
  p.AddText("t (ns) = (%.2f #pm %.2f) + (%.3f #pm %.3f) * dist [cm]" % (params[0],paramerrs[0],params[1],paramerrs[1]))
  p.SetFillStyle(0)
  p.SetBorderSize(0)
  p.Draw()

  raw_input()



def fit_reflection_fraction():
  # REFLECTED FRACTION
  mg = ROOT.TMultiGraph("mg","mg")
  l = ROOT.TLegend(0.65,0.65,0.85,0.85)
  color_index = 1
  rfs = []
  for straw in straws:
    x = array('d',distances[straw])
    y = array('d',[results[straw]['cal'][i]['scale1']/(np.exp(-(2*straw_lengths[straw]-2*distances[straw][i])/attenuation_length)) for i in range(len(distances[straw]))])
    gc = ROOT.TGraph(len(distances[straw]),x,y)
    for i in range(len(y)):
      rfs.append(y[i])
    x = array('d',distances[straw])
    y = array('d',[results[straw]['hv'][i]['scale1']/(np.exp(-(2*straw_lengths[straw]-2*distances[straw][i])/attenuation_length)) for i in range(len(distances[straw]))])
    gh = ROOT.TGraph(len(distances[straw]),x,y)
    for i in range(len(y)):
      rfs.append(y[i])
    gc.SetMarkerStyle(21)
    gh.SetMarkerStyle(22)
    gc.SetMarkerColor(color_index)
    gc.SetLineColor(color_index)
    gh.SetMarkerColor(color_index)
    gh.SetLineColor(color_index)
    color_index += 1
    l.AddEntry(gc,"Cal %d" % straw,"p")
    l.AddEntry(gh,"HV %d" % straw,"p")
  
    mg.Add(gc)
    mg.Add(gh)

    print straw,"Mean reflection fraction:",np.mean(np.append([results[straw]['cal'][i]['scale1']/(np.exp(-(2*straw_lengths[straw]-2*distances[straw][i])/attenuation_length)) for i in range(len(distances[straw]))],[results[straw]['hv'][i]['scale1']/(np.exp(-(2*straw_lengths[straw]-2*distances[straw][i])/attenuation_length)) for i in range(len(distances[straw]))])),"+-",np.std(np.append([results[straw]['cal'][i]['scale1']/(np.exp(-(2*straw_lengths[straw]-2*distances[straw][i])/attenuation_length)) for i in range(len(distances[straw]))],[results[straw]['hv'][i]['scale1']/(np.exp(-(2*straw_lengths[straw]-2*distances[straw][i])/attenuation_length)) for i in range(len(distances[straw]))]))
  fr = mg.Fit("pol1","S")
  mg.SetTitle("Calculated fraction of signal reflected in first reflection (A.L. = %.1f cm)" % attenuation_length)
  mg.Draw("ALP")
  l.Draw()
  mg.GetXaxis().SetTitle("Distance from straw end (cm)")
  mg.GetYaxis().SetTitle("Reflected fraction")
  p = ROOT.TPaveText(0.4,0.18,0.85,0.25,"NDC NB")
  p.AddText("fraction = (%.2f #pm %.2f)" % (np.mean(rfs),np.std(rfs)))
  p.SetFillStyle(0)
  p.SetBorderSize(0)
  p.Draw()
  raw_input()


def fit_reflection_fraction2():
  # REFLECTED FRACTION
  mg = ROOT.TMultiGraph("mg","mg")
  l = ROOT.TLegend(0.65,0.65,0.85,0.85)
  color_index = 1
  frs = {}
  f1s = {}
  for straw in straws:
    x = array('d',[2*straw_lengths[straw]-2*distances[straw][i] for i in range(len(distances[straw]))])
    y = array('d',[results[straw]['cal'][i]['scale1'] for i in range(len(distances[straw]))])
    gc = ROOT.TGraph(len(distances[straw]),x,y)
    f1 = ROOT.TF1("f1_%d" % straw,"[0]*exp(-x/[1])",x[0],x[-1])
    #f1 = ROOT.TF1("f1_%d" % straw,"[0]-x/[1]",x[0],x[-1])
    f1.SetParameter(0,y[0])
    f1.SetParameter(1,200)
    fr = gc.Fit("f1_%s" % straw,"SQ0")
    f1s[straw] = [f1]
    frs[straw] = [fr]

    print straw,"cal: rf = ",fr.GetParams()[0]," al = ",fr.GetParams()[1]

    x = array('d',[2*straw_lengths[straw]-2*distances[straw][i] for i in range(len(distances[straw]))])
    y = array('d',[results[straw]['hv'][i]['scale1'] for i in range(len(distances[straw]))])
    gh = ROOT.TGraph(len(distances[straw]),x,y)
    f2 = ROOT.TF1("f2_%d" % straw,"[0]*exp(-x/[1])",x[0],x[-1])
    #f2 = ROOT.TF1("f2_%d" % straw,"[0]-x/[1]",x[0],x[-1])
    f2.SetParameter(0,y[0])
    f2.SetParameter(1,200)
    fr2 = gh.Fit("f2_%d" % straw,"SQ0")
    f1s[straw].append(f2)
    frs[straw].append(fr2)

    print straw,"hv: rf = ",fr2.GetParams()[0]," al = ",fr2.GetParams()[1]

    gc.SetMarkerStyle(21)
    gh.SetMarkerStyle(22)
    gc.SetMarkerColor(color_index)
    gc.SetLineColor(color_index)
    f1.SetLineColor(color_index)
    f1.SetLineStyle(2)
    gh.SetMarkerColor(color_index)
    gh.SetLineColor(color_index)
    f2.SetLineColor(color_index)
    f2.SetLineStyle(3)
    color_index += 1
    l.AddEntry(gc,"Cal %d" % straw,"p")
    l.AddEntry(gh,"HV %d" % straw,"p")
  
    mg.Add(gc)
    mg.Add(gh)

  mg.SetTitle("Calculated fraction of signal reflected in first reflection")
  mg.Draw("ALP")
  l.Draw()
  rfs = []
  als = []
  for straw in straws:
    f1s[straw][0].Draw("same")
    f1s[straw][1].Draw("same")
    rfs.append(frs[straw][0].GetParams()[0])
    rfs.append(frs[straw][1].GetParams()[0])
    als.append(frs[straw][0].GetParams()[1])
    als.append(frs[straw][1].GetParams()[1])
  mg.GetXaxis().SetTitle("Distance travelled by reflection (cm)")
  mg.GetYaxis().SetTitle("Reflected fraction")
  p = ROOT.TPaveText(0.4,0.18,0.85,0.25,"NDC NB")
  p.AddText("fraction = (%.2f #pm %.2f)" % (np.mean(rfs),np.std(rfs)))
  p.AddText("attenuation length = (%.1f #pm %.1f)" % (np.mean(als),np.std(als)))
  p.SetFillStyle(0)
  p.SetBorderSize(0)
  p.Draw()
  raw_input()




def fit_reflection2_time():
  # REFLECTION TIME SHIFT AND VELOCITY for SECOND REFLECTION
  mg = ROOT.TMultiGraph("mg","mg")
  l = ROOT.TLegend(0.65,0.65,0.85,0.85)
  color_index = 1

  # have to assume some phase shift time
  phase_shift = 4.62 #from fit_reflection_time()
  velocities = {}
  paramerrs = {}
  colors = {}
  for straw in straws:
    rfs = []
    x = array('d',distances[straw])
    y = array('d',[results[straw]['hv'][i]['time2'] for i in range(len(distances[straw]))])
    gh = ROOT.TGraph(len(distances[straw]),x,y)
    for i in range(len(y)):
      rfs.append(y[i])
    x = array('d',distances[straw])
    y = array('d',[results[straw]['cal'][i]['time2'] for i in range(len(distances[straw]))])
    gc = ROOT.TGraph(len(distances[straw]),x,y)
    for i in range(len(y)):
      rfs.append(y[i])
    colors[straw] = color_index

    gc.SetMarkerStyle(21)
    gh.SetMarkerStyle(22)
    gc.SetMarkerColor(color_index)
    gc.SetLineColor(color_index)
    gh.SetMarkerColor(color_index)
    gh.SetLineColor(color_index)
    color_index += 1
    l.AddEntry(gc,"Cal %d" % straw,"p")
    l.AddEntry(gh,"HV %d" % straw,"p")
  
    mg.Add(gc)
    mg.Add(gh)

    velocities[straw] = 2*straw_lengths[straw]/(np.mean(rfs)-2*phase_shift)
    paramerrs[straw] = velocities[straw] - 2*straw_lengths[straw]/(np.mean(rfs)+np.std(rfs)-2*phase_shift)
  #fr = mg.Fit("pol1","S")
  mg.SetTitle("Time delay of second reflection")
  mg.Draw("ALP")
  l.Draw()
  mg.GetXaxis().SetTitle("Distance from straw end (cm)")
  mg.GetYaxis().SetTitle("#delta t_{reflection} (ns)")
  mg.GetYaxis().SetRangeUser(4,22)

  p = ROOT.TPaveText(0.4,0.18,0.85,0.25,"NDC NB")
  for straw in straws:
    tt = p.AddText("Straw %d: Eff. vel. = %.2f #pm %.2f (cm/ns)" % (straw,velocities[straw],paramerrs[straw]))
    tt.SetTextColor(colors[straw])
    tt.Draw()
#  p.AddText("fraction = (%.2f \pm %.2f)" % (np.mean(rfs),np.std(rfs)))
#  p.AddText("fraction = (%.2f \pm %.2f)" % (np.mean(rfs),np.std(rfs)))
#  p.AddText("fraction = (%.2f \pm %.2f)" % (np.mean(rfs),np.std(rfs)))
  p.SetFillStyle(0)
  p.SetBorderSize(0)
  p.Draw()
 
  
  #params = fr.GetParams()
  #paramerrs = fr.GetErrors()
  #print "2nd reflection t [ns] = (%f +- %f) + (%f +- %f) * dist" % (params[0],paramerrs[0],params[1],paramerrs[1])
  raw_input()


def fit_reflection2_fraction():
  # REFLECTED FRACTION for SECOND REFLECTION
  mg = ROOT.TMultiGraph("mg","mg")
  l = ROOT.TLegend(0.65,0.65,0.85,0.85)
  color_index = 1
  rfs = []
  for straw in straws:
    x = array('d',distances[straw])
    y = array('d',[np.sqrt(results[straw]['cal'][i]['scale2']/(np.exp(-(2*straw_lengths[straw])/attenuation_length))) for i in range(len(distances[straw]))])
    gc = ROOT.TGraph(len(distances[straw]),x,y)
    for i in range(len(y)):
      rfs.append(y[i])
    x = array('d',distances[straw])
    y = array('d',[np.sqrt(results[straw]['hv'][i]['scale1']/(np.exp(-(2*straw_lengths[straw])/attenuation_length))) for i in range(len(distances[straw]))])
    gh = ROOT.TGraph(len(distances[straw]),x,y)
    for i in range(len(y)):
      rfs.append(y[i])

    gc.SetMarkerStyle(21)
    gh.SetMarkerStyle(22)
    gc.SetMarkerColor(color_index)
    gc.SetLineColor(color_index)
    gh.SetMarkerColor(color_index)
    gh.SetLineColor(color_index)
    color_index += 1
    l.AddEntry(gc,"Cal %d" % straw,"p")
    l.AddEntry(gh,"HV %d" % straw,"p")
  
    mg.Add(gc)
    mg.Add(gh)
#  fr = mg.Fit("pol1","S")
  mg.SetTitle("Calculated fraction of signal reflected in second reflection (A.L. = %.1f cm)" % attenuation_length)
  mg.Draw("ALP")
  l.Draw()
  mg.GetXaxis().SetTitle("Distance from straw end (cm)")
  mg.GetYaxis().SetTitle("Reflected fraction")
  p = ROOT.TPaveText(0.4,0.18,0.85,0.25,"NDC NB")
  p.AddText("fraction = (%.2f #pm %.2f)" % (np.mean(rfs),np.std(rfs)))
  p.SetFillStyle(0)
  p.SetBorderSize(0)
  p.Draw()

  raw_input()




def fit_scaling_integral():
  # FIT SCALING INTEGRAL


#  hnlls = []
#  c = ROOT.TCanvas("c","c",800,600)
#  l = ROOT.TLegend(0.6,0.6,0.8,0.8)
#  for i in range(len(distances[straw])):
#    hnlls.append(ROOT.TH1F("hnllcal_%d" % i,"hnllcal_%d" % i,750,0,150))
#    temp = copy.copy(results[straw]['cal'][i])
#    temp['scale1'] = 0
#    temp['scale2'] = 0
#    histo(hnlls[-1],**temp)
#    peaks.append(hnlls[-1].GetMaximum())
#    ints.append(hnlls[-1].Integral(0,650))
#    hnlls[-1].SetLineColor(i+1)
#    l.AddEntry(hnlls[-1],"%f" % distances[straw][i],"l")
#    if i == 0:
#      hnlls[-1].Draw()
#    else:
#      hnlls[-1].Draw("same")
#  l.Draw()
#  hnlls2 = []
#  c2 = ROOT.TCanvas("c2","c2",800,600)
#  l2 = ROOT.TLegend(0.6,0.6,0.8,0.8)
#  for i in range(len(distances[straw])):
#    hnlls2.append(ROOT.TH1F("hnll2cal_%d" % i,"hnll2cal_%d" % i,750,0,150))
#    temp = copy.copy(results[straw]['cal'][i])
#    temp['scale1'] = 0
#    temp['scale2'] = 0
#    temp['sigma'] = sigma
#    temp['t0'] = t0[i]
#    temp['scale'] = 1
#    histo(hnlls2[-1],**temp)
#    hnlls2[-1].SetLineColor(i+1)
#    l2.AddEntry(hnlls2[-1],"%f" % distances[straw][i],"l")
#    if i == 0:
#      hnlls2[-1].Draw()
#    else:
#      hnlls2[-1].Draw("same")
#  l2.Draw()
#  raw_input()

  import copy
  mg = ROOT.TMultiGraph("mg","mg")
  l = ROOT.TLegend(0.65,0.15,0.85,0.35)
  color_index = 1
  for straw in straws:
    print "STRAW",straw
    # HAVE TO ASSUME FUNCTIONAL FORM WE ARE USING FOR SIGMA AND TAU FIXME
    #sigma = 2.384
    #t0 = 2.105488 + 0.047196*distances[straw]
    sigma = 2.31
    t0 = 1.06 + 0.032*distances[straw]

    x = array('d',distances[straw])
    y = array('d') 
    for i in range(len(distances[straw])):
      temp = copy.copy(results[straw]['cal'][i])
      temp['scale1'] = 0
      temp['scale2'] = 0
      histo(hnll,**temp)
      pmp_orig = hnll.GetMaximum()
      temp['sigma'] = sigma
      temp['t0'] = t0[i]
      temp['scale'] = 1
      histo(hnll,**temp)
      pmp_fit = hnll.GetMaximum()
      y.append(pmp_orig/pmp_fit)
    print "y",y
    gc = ROOT.TGraph(len(distances[straw]),x,y)
    x2 = array('d',distances[straw])
    y2 = array('d') 
    for i in range(len(distances[straw])):
      temp = copy.copy(results[straw]['hv'][i])
      temp['scale1'] = 0
      temp['scale2'] = 0
      histo(hnll,**temp)
      pmp_orig = hnll.GetMaximum()
      temp['sigma'] = sigma
      temp['t0'] = t0[i]
      temp['scale'] = 1
      histo(hnll,**temp)
      pmp_fit = hnll.GetMaximum()
      y2.append(pmp_orig/pmp_fit)
    print "y2",y2
    gh = ROOT.TGraph(len(distances[straw]),x2,y2)

    gc.SetMarkerStyle(21)
    gh.SetMarkerStyle(22)
    gc.SetMarkerColor(color_index)
    gc.SetLineColor(color_index)
    gh.SetMarkerColor(color_index)
    gh.SetLineColor(color_index)
    color_index += 1
    l.AddEntry(gc,"Cal %d" % straw,"p")
    l.AddEntry(gh,"HV %d" % straw,"p")

#    print [np.mean([y[i],y2[i]]) for i in range(len(y))]
#    print [np.std([y[i],y2[i]]) for i in range(len(y))]
   
    mg.Add(gc)
    mg.Add(gh)
#  fr = mg.Fit("pol1","SQ")
  mg.Draw("ALP")
  l.Draw()
  mg.GetXaxis().SetTitle("Distance from straw end (cm)")
  mg.GetYaxis().SetTitle("Scaling factor")
  
#  params = fr.GetParams()
#  paramerrs = fr.GetErrors()
#  print "scaling factor = (%f +- %f) + (%f +- %f) * dist" % (params[0],paramerrs[0],params[1],paramerrs[1])
  raw_input()
  
def fit_scaling_integral2():
  import copy
  # use fit attenuation length 
  fit_alength = 139.
  mg = ROOT.TMultiGraph("mg","mg")
  l = ROOT.TLegend(0.65,0.15,0.85,0.35)
  color_index = 1
  for straw in straws:
    x = array('d',distances[straw])
    y = array('d') 
    for i in range(len(distances[straw])):
      temp = copy.copy(results[straw]['cal'][i])
      temp['scale1'] = 0
      temp['scale2'] = 0
      histo(hnll,**temp)
      pmp_orig = hnll.GetMaximum()
      y.append(pmp_orig/np.exp(-distances[straw][i]/fit_alength))
    print "y",y
    gc = ROOT.TGraph(len(distances[straw]),x,y)
    x2 = array('d',distances[straw])
    y2 = array('d') 
    for i in range(len(distances[straw])):
      temp = copy.copy(results[straw]['hv'][i])
      temp['scale1'] = 0
      temp['scale2'] = 0
      histo(hnll,**temp)
      pmp_orig = hnll.GetMaximum()
      y2.append(pmp_orig/np.exp(-distances[straw][i]/fit_alength))
    print "y2",y2
    gh = ROOT.TGraph(len(distances[straw]),x2,y2)

    gc.SetMarkerStyle(21)
    gh.SetMarkerStyle(22)
    gc.SetMarkerColor(color_index)
    gc.SetLineColor(color_index)
    gh.SetMarkerColor(color_index)
    gh.SetLineColor(color_index)
    color_index += 1
    l.AddEntry(gc,"Cal %d" % straw,"p")
    l.AddEntry(gh,"HV %d" % straw,"p")

#    print [np.mean([y[i],y2[i]]) for i in range(len(y))]
#    print [np.std([y[i],y2[i]]) for i in range(len(y))]
   
    mg.Add(gc)
    mg.Add(gh)
#  fr = mg.Fit("pol1","SQ")
  mg.Draw("ALP")
  l.Draw()
  mg.GetXaxis().SetTitle("Distance from straw end (cm)")
  mg.GetYaxis().SetTitle("Scaling factor")
  
#  params = fr.GetParams()
#  paramerrs = fr.GetErrors()
#  print "scaling factor = (%f +- %f) + (%f +- %f) * dist" % (params[0],paramerrs[0],params[1],paramerrs[1])
  raw_input()
 
 
def fit_scaling_integral3():
  import copy
  # use fit attenuation length 
  fit_alength = 137.7
  x = np.linspace(0,120,5)
  print x
  y = []
  temp_results = {"meant": 30, "sigma": 0, "t0": 0, "pole0": 160, "pole1": 6, "zero": 1.2, "scale1": 0, "scale2": 0, "time1": 0, "time2": 0, "scale": 1}
  for i in range(len(x)):
    print i
    sigma = 2.42
    t0 = 0.8 + 0.015*x[i]
    temp_results["sigma"] = sigma
    temp_results["t0"] = t0
    histo(hnll,**temp_results)
    pmp = hnll.GetMaximum()
    if i == 0:
      original_pmp = pmp
      y.append(temp_results['scale'])
    else:
      goal_pmp = original_pmp*np.exp(-x[i]/fit_alength)
      y.append(temp_results['scale']*goal_pmp/pmp)
  y = np.array(y)
  print y
  g = ROOT.TGraph(len(x),x,y)
  g.SetMarkerStyle(21)
  g.Draw("ALP")
  g.GetXaxis().SetTitle("Distance from straw end (cm)")
  g.GetYaxis().SetTitle("Scaling factor")
  raw_input()
 

if __name__ == "__main__":
  calc_sigma()
  fit_t0()
  fit_attenuation_length()
  fit_reflection_time()
  fit_reflection_fraction()
  fit_reflection_fraction2()
  fit_reflection2_time()
  fit_reflection2_fraction()
  fit_prop_velocity()
####  fit_scaling_integral()
  fit_scaling_integral3()
####  fit_scaling_integral2()
####  fit_scaling()
