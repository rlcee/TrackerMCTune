#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TApplication.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TMath.h>
#include <TProfile.h>


int main(int argc, char** argv)
{
  std::string filename[2];
  filename[0] = "";
  filename[1] = "";
  std::string help = "./tot_compare -d <data filename> -s <sim filename>";
  std::string inputcommand = std::string(argv[0]);
  for (int i=1;i<argc;i++)
    inputcommand += " " + std::string(argv[i]);
  

  int c;
  while ((c = getopt (argc, argv, "hd:s:")) != -1){
    switch (c){
      case 'h':
        std::cout << help << std::endl;
        return 0;
      case 'd':
        filename[0] = std::string(optarg);      
        break;
      case 's':
        filename[1] = std::string(optarg);
        break;
      case '?':
        if (optopt == 'd' || optopt == 's')
          std::cout << "Option -" << optopt << " requires an argument." << std::endl;
        else
          std::cout << "Unknown option `-" << optopt << "'." << std::endl;
        return 1;
    }
  }

  if (filename[0].size() == 0 && filename[1].size() == 0 && optind == argc-2){
    for (int index = optind; index < argc; index++){
      filename[index-optind] = std::string(argv[index]);
    }
  }

  if (filename[0].size() == 0 || filename[1].size() == 0){
    std::cout << help << std::endl;
    return 1;
  }
 
  int argc2 = 0;char **argv2;TApplication theApp("tapp", &argc2, argv2);

  std::vector<double> times[2];
  std::vector<double> docas[2];
  std::vector<double> tots[2];

  TH2F* h[2];

  for (int ifile=0;ifile<2;ifile++){
    TFile *f = new TFile(filename[ifile].c_str());
    TTree *t = (TTree*) f->Get("tot");
    double t_time, t_tot, t_doca;
    t->SetBranchAddress("tot",&t_tot);
    t->SetBranchAddress("time",&t_time);
    t->SetBranchAddress("doca",&t_doca);
    for (int i=0;i<t->GetEntries();i++){
      t->GetEntry(i);
      tots[ifile].push_back(t_tot);
      times[ifile].push_back(t_time);
      docas[ifile].push_back(t_doca);
    }

    h[ifile] = new TH2F(TString::Format("tot_%d",ifile),"tot",32,0,64,120,-20,100);

    for (int j=0;j<tots[ifile].size();j++){
      h[ifile]->Fill(tots[ifile][j],times[ifile][j]);
    }
  }

  TCanvas *cscatter = new TCanvas("cscatter","scatter",600,600);
  h[0]->SetStats(0);
  h[0]->SetTitle("");
  h[0]->GetXaxis()->SetTitle("Time over threshold (ns)");
  h[0]->GetYaxis()->SetTitle("Drift time (ns)");
  h[0]->Draw();
  h[1]->SetMarkerColor(kRed);
  h[1]->Draw("same");

  TLegend *l = new TLegend(0.55,0.65,0.85,0.85);
  l->AddEntry(h[0],"8-Straw Prototype Data","P");
  l->AddEntry(h[1],"G4 + Straw Simulation","P");
  l->SetBorderSize(0);
  l->Draw();


  TCanvas *chist = new TCanvas("chist","chist",600,600);
  TProfile *hdata = (TProfile*) h[0]->ProfileX("hdata",1,-1,"S");
  TProfile *hsim = (TProfile*) h[1]->ProfileX("hsim",1,-1,"S");
  hdata->SetLineColor(kBlue);
  hdata->SetTitle("");
  hdata->SetStats(0);
  hdata->GetXaxis()->SetTitle("Time over threshold (ns)");
  hdata->GetYaxis()->SetTitle("Drift time (ns)");
  hdata->SetMarkerStyle(22);
//  hdata->GetXaxis()->SetRangeUser(4,50);
  hdata->Draw();
  hsim->SetLineColor(kRed);
  hsim->SetFillColor(kRed);
  hsim->SetMarkerColor(kRed);
  hsim->SetFillStyle(3001);
  hsim->Draw("same E2");
  TProfile *hsim2 = (TProfile*) hsim->Clone("hsim2");
  hsim2->SetLineColor(kRed);
  hsim2->SetFillStyle(0);
  hsim2->Draw("hist same");


  TLegend *l2 = new TLegend(0.55,0.65,0.85,0.85);
  l2->AddEntry(hdata,"8-Straw Prototype Data","PL");
  l2->AddEntry(hsim,"G4 + Straw Simulation","FL");
  l2->SetBorderSize(0);
  l2->Draw();

  theApp.Run();

  return 0;

}

