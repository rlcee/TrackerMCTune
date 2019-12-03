#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>

#include <Minuit2/FCNBase.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnSimplex.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>
#include <Minuit2/FunctionMinimum.h>
#include <TTree.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TError.h>
#include <TLine.h>
#include <TStyle.h>
#include <TFile.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TCanvas.h>

#include "prototype/Utils.hh"
#include "fit/StrawFit.hh"

int main(int argc, char** argv)
{
  std::string filename = "";
  std::string outname = "";
  std::string help = "./tot_sim -i <input filename> -o <output filename>";
  std::string inputcommand = std::string(argv[0]);
  for (int i=1;i<argc;i++)
    inputcommand += " " + std::string(argv[i]);
  
  int argc2 = 0;char **argv2;TApplication theApp("tapp", &argc2, argv2);

  int c;
  while ((c = getopt (argc, argv, "hi:o:")) != -1){
    switch (c){
      case 'h':
        std::cout << help << std::endl;
        return 0;
      case 'i':
        filename = std::string(optarg);
        break;
      case 'o':
        outname = std::string(optarg);
        break;
      case '?':
        if (optopt == 'i' || optopt == 'o')
          std::cout << "Option -" << optopt << " requires an argument." << std::endl;
        else
          std::cout << "Unknown option `-" << (char)optopt << "'." << std::endl;
        return 1;
    }
  }

  if (filename.size() == 0 || outname.size() == 0){
    std::cout << help << std::endl;
    return 1;
  }
 
  gStyle->SetOptFit(1);

  std::vector<double> times;
  std::vector<double> tots;
  std::vector<double> docas;
  TFile *f = new TFile(filename.c_str());
  /*
  TDirectory *d = (TDirectory*) f->Get("SHD");
  TTree *t = (TTree*) d->Get("shdiag");
  Int_t straw, panel;
  Int_t mcproc;
  Float_t time[2], mcshd;
  Float_t ewm;
  Double_t mcsptime;
  Float_t tot[2];

  t->SetBranchAddress("straw",&straw);
  t->SetBranchAddress("panel",&panel);
  t->SetBranchAddress("time",&time);
  t->SetBranchAddress("tot",&tot);
  t->SetBranchAddress("mcsptime",&mcsptime);
  t->SetBranchAddress("mcshd",&mcshd);
  t->SetBranchAddress("ewmoffset",&ewm);
  t->SetBranchAddress("mcproc",&mcproc);

  for (int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    if (mcproc != 56)
      continue;
    // the prototype is only the top 8 straws
    // and our generator is pointed at panel 1
    if (straw < 88 || panel != 1)
      continue;
    times.push_back(std::min(time[0],time[1]) - mcsptime + ewm);
    tots.push_back((tot[0]+tot[1])/2.);
  }
  */

  TDirectory *d = (TDirectory*) f->Get("makeSD");
  TTree *t = (TTree*) d->Get("sddiag");
  Int_t straw, panel;
  Int_t tdc[2], tot[2], mcproc;
  Double_t mctime;
  Float_t mcdca, ecptime, wdist[2], ewmoffset;
  t->SetBranchAddress("mcproc",&mcproc);
  t->SetBranchAddress("straw",&straw);
  t->SetBranchAddress("panel",&panel);
  t->SetBranchAddress("tdc",&tdc);
  t->SetBranchAddress("tot",&tot);
  t->SetBranchAddress("mctime",&mctime);
  t->SetBranchAddress("mcdca",&mcdca);
  t->SetBranchAddress("wdist",&wdist);
  t->SetBranchAddress("ecptime",&ecptime);
  t->SetBranchAddress("ewmoffset",&ewmoffset);


  for (int i=0;i<t->GetEntries();i++){
    t->GetEntry(i);
    //FIXME should i cut on mcproc?
//    if (mcproc != 56)
//      continue;
    if (straw < 88 || panel != 1)
      continue;
    times.push_back((tdc[0]+tdc[1])/2.0*2*0.015625 - (fmod(mctime,1695) + ecptime - ewmoffset));
    docas.push_back(mcdca);
    double totcal = tot[0] * 4 + (0x7F - tdc[0] & 0x7F)*0.015625*2;
    double tothv  = tot[1] * 4 + (0x7F - tdc[1] & 0x7F)*0.015625*2;
//    std::cout << "Cal: " << tot[0]*4 << " " << totcal << std::endl;
//    std::cout << "HV: " << tot[1]*4 << " " << tothv << std::endl;
    tots.push_back((totcal+tothv)/2.); 
  }


  TFile *fout = new TFile(outname.c_str(),"RECREATE");
  TTree *tout = new TTree("tot","tot");
  double t_time, t_tot, t_doca;
  tout->Branch("time",&t_time,"time/D");
  tout->Branch("tot",&t_tot,"tot/D");
  tout->Branch("doca",&t_doca,"doca/D");
  for (int i=0;i<times.size();i++){
    t_tot = tots[i];
    t_time = times[i];
    t_doca = docas[i];
    tout->Fill();
  }
  tout->Write();

  TObjString *tl = new TObjString(inputcommand.c_str());
  fout->WriteObject(tl,"command");
  fout->Close();
















  TH2F *h1 = new TH2F("h1","h1",32,0,64,120,-20,100);
  TH1F *h1d = new TH1F("h1d","h1d",32,0,64);
//  timeoffset = -20;

  for (int i=0;i<times.size();i++){
    h1->Fill(tots[i],times[i]);
  }

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  TProfile *hp = (TProfile*) h1->ProfileX("hp",1,-1,"S");
  hp->SetLineColor(kBlue);
  hp->SetTitle("");
  hp->SetStats(0);
  hp->GetXaxis()->SetTitle("Time over threshold (ns)");
  hp->GetYaxis()->SetTitle("Drift time (ns)");
  hp->SetMarkerStyle(22);
  hp->GetXaxis()->SetRangeUser(0,50);
//  hp->Draw();
  h1->Draw();
//  h1d->Draw();

  theApp.Run();

  return 0;
}
