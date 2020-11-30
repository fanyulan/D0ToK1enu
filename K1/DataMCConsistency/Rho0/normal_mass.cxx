#include<cstdio>
#include<string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <iomanip>

#include "RooAbsReal.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "TMath.h"


using namespace RooFit;
void normal_mass(){

TStyle *bes3Style= new TStyle("bes3","bes3 style");
Int_t icol=0;
bes3Style->SetFrameBorderMode(icol);
bes3Style->SetCanvasBorderMode(icol);
bes3Style->SetPadBorderMode(icol);
bes3Style->SetPadColor(icol);
bes3Style->SetCanvasColor(icol);
bes3Style->SetStatColor(icol);
bes3Style->SetTitleFillColor(icol);
bes3Style->SetPalette(1); 
bes3Style->SetPaperSize(TStyle::kUSLetter);
bes3Style->SetPadTopMargin(.12);
bes3Style->SetPadLeftMargin(.15);
bes3Style->SetPadRightMargin(.08);
bes3Style->SetPadBottomMargin(.15);
Int_t font=22; 
Double_t tsize=0.05;
bes3Style->SetTextFont(font);
bes3Style->SetTextSize(tsize);
bes3Style->SetLabelFont(font,"xyz");
bes3Style->SetLabelSize(tsize,"xyz");
bes3Style->SetLabelOffset(0.01,"xyz");
bes3Style->SetTitleFont(font,"xyz");
bes3Style->SetTitleSize(tsize,"xyz");
bes3Style->SetTitleXOffset(0.9);
bes3Style->SetTitleYOffset(1.1);
bes3Style->SetTitleBorderSize(2.);
bes3Style->SetMarkerStyle(0);
bes3Style->SetMarkerSize(0.8);
bes3Style->SetFrameBorderMode (0.);
bes3Style->SetFrameLineWidth  (2.);
bes3Style->SetLineWidth(1.0);
bes3Style->SetHistLineWidth(2.);
bes3Style->SetLineStyleString(2,"[12 12]");
bes3Style->SetErrorX(0.001);
bes3Style->SetOptTitle(0);
bes3Style->SetOptStat(0);
bes3Style->SetLineStyleString(2,"[30 10]");
bes3Style->SetLineStyleString(3,"[4 8]");
bes3Style->SetLineStyleString(4,"[15 12 4 12]");
bes3Style->SetLineStyleString(5,"[15 15]");
bes3Style->SetLineStyleString(6,"[15 12 4 12 4 12]");
bes3Style->SetOptDate(0);
bes3Style->SetDateY(.98);
bes3Style->SetStripDecimals(kFALSE);
bes3Style->SetEndErrorSize(0.0);
gROOT->SetStyle("bes3");
gROOT->ForceStyle();


using namespace RooFit;
gSystem->Load("libRooFit");
using namespace RooFit;
gStyle->SetLabelSize(0.07,"xyz");
gStyle->SetNdivisions(507,"xyz");
gROOT->SetStyle("Plain");
gStyle->SetLabelSize(0.06,"xyz");
gStyle->SetNdivisions(405,"xyz");
gStyle->SetPadTopMargin(.10);
gStyle->SetPadLeftMargin(.15);
gStyle->SetPadRightMargin(.05);
gStyle->SetPadBottomMargin(.15);
gStyle->SetTitleSize(0.06,"xyz");
gStyle->SetOptTitle(0);
gStyle->SetMarkerSize(0.5);
gStyle->SetTitle("");
gStyle->SetOptStat(0);
gStyle->SetEndErrorSize(0);

      TChain tdata("tagD");  
      tdata.Add("/scratchfs/bes/fangyl/data1/data/all/*.root");//data
      TChain tsig("tagD");
      tsig.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/PDG2020/DDecay_root/red_truth/truth*.root");
      //tsig.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/BELLE/fit2/red_truth/truth*.root");
      TChain tbkg("tagD");
      tbkg.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/new/inclusiveMC/umiss2fit/reduce_truth/root/truth*.root");

      TCut cut_signal="isTrueTag&&isSignal&&abs(KpTrueID)==321&&abs(pi1TrueID)==211&&abs(pi2TrueID)==211&&abs(elecTrueID)==11&&abs(elecMomID)==421&&(((mcmodea==666||mcmodea==222)&&mcmode1==103)||((mcmodeb==-666||mcmodeb==-222)&&mcmode2==-103))";
      TCut cuta = "isbest==1&&type==1&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)&&m_pipi>0.31&&fabs(m_pipi-0.497611)>0.01&&fabs(deltaE_pipzswp)>0.012";
      //TCut cutb = "fabs(Umissfit)<0.05&&mhad>1.1&&mhad<1.4";
      TCut cutb = "fabs(Umissfit)<0.025&&mhad>1.172&&mhad<1.372";
      TCut cutc = "elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge";
      TCut cut0 = "mode==0&&clveto==1&&deltaE>-0.029&&deltaE<0.027";
      TCut cut1 = "mode==1&&deltaE>-0.069&&deltaE<0.038";
      TCut cut3 = "mode==3&&deltaE>-0.031&&deltaE<0.028";
      TCut K_2  = "mcmode1==104||mcmode2==-104";
     
      TCut Krho     = "abs(pi1MomID)==113&&abs(pi2MomID)==113";//K rho
      TCut K1430pi  = "abs(KpMomID)==10311&&abs(pi1MomID)==10311";//K0*0(1430) pi
      TCut K892pi   = "abs(KpMomID)==313&&abs(pi1MomID)==313";
      TCut Komega   = "abs(pi1MomID)==223&&abs(pi2MomID)==223"; 

      double bins = 25, x1min=0.6, x1max=1.2;
      hsig_1         = new TH1F("hsig_1","",        bins,x1min,x1max);
      hsig_1Krho     = new TH1F("hsig_1Krho","",    bins,x1min,x1max);
      hsig_1K1430pi  = new TH1F("hsig_1K1430pi","", bins,x1min,x1max);
      hsig_1K892pi   = new TH1F("hsig_1K892pi","",  bins,x1min,x1max);
      hsig_1Komega   = new TH1F("hsig_1Komega","",  bins,x1min,x1max);
      hsig_1other    = new TH1F("hsig_1other","",   bins,x1min,x1max);
      hbkg_1         = new TH1F("hbkg_1","",        bins,x1min,x1max);
      hdata_1        = new TH1F("hdata_1","",       bins,x1min,x1max);
      double x2min=0.2, x2max=0.9;
      hsig_2         = new TH1F("hsig_2","",        bins,x2min,x2max);
      hsig_2Krho     = new TH1F("hsig_2Krho","",    bins,x2min,x2max);
      hsig_2K1430pi  = new TH1F("hsig_2K1430pi","", bins,x2min,x2max);
      hsig_2K892pi   = new TH1F("hsig_2K892pi","",  bins,x2min,x2max);
      hsig_2Komega   = new TH1F("hsig_2Komega","",  bins,x2min,x2max);
      hsig_2other    = new TH1F("hsig_2other","",   bins,x2min,x2max);
      hbkg_2         = new TH1F("hbkg_2","",        bins,x2min,x2max);
      hdata_2        = new TH1F("hdata_2","",       bins,x2min,x2max);

      tsig.Draw("m_Kpi>>hsig_1",       (cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal);  
      tsig.Draw("m_Kpi>>hsig_1Krho",   (cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal&&Krho);
      tsig.Draw("m_Kpi>>hsig_1K1430pi",(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal&&K1430pi);
      tsig.Draw("m_Kpi>>hsig_1K892pi", (cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal&&K892pi);
      tsig.Draw("m_Kpi>>hsig_1Komega", (cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal&&Komega);

      tdata.Draw("m_Kpi>>hdata_1",     (cut0||cut1||cut3)&&cuta&&cutb&&cutc);
      tbkg.Draw("m_Kpi>>hbkg_1",       (cut0||cut1||cut3)&&cuta&&cutb&&cutc&&!cut_signal&&!K_2);
 
      tsig.Draw("m_pipi>>hsig_2"  ,    (cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal);
      tsig.Draw("m_pipi>>hsig_2Krho",   (cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal&&Krho);
      tsig.Draw("m_pipi>>hsig_2K1430pi",(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal&&K1430pi);
      tsig.Draw("m_pipi>>hsig_2K892pi", (cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal&&K892pi);
      tsig.Draw("m_pipi>>hsig_2Komega", (cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal&&Komega);
      
      tdata.Draw("m_pipi>>hdata_2",    (cut0||cut1||cut3)&&cuta&&cutb&&cutc);
      tbkg.Draw("m_pipi>>hbkg_2",      (cut0||cut1||cut3)&&cuta&&cutb&&cutc&&!cut_signal&&!K_2);
      
      double nbkg1 = hbkg_1->GetEntries(), nbkg2 = hbkg_2->GetEntries();
      double ndata1= hdata_1->GetEntries(),ndata2= hdata_2->GetEntries();
      double Nsig1 = ndata1 - nbkg1/10.8, Nsig2 = ndata2 - nbkg2/10.8;
      double nsig1 = hsig_1->GetEntries(), nsig2 = hsig_2->GetEntries();     
      double nsig1_Krho    = hsig_1Krho->GetEntries(),    nsig2_Krho = hsig_2Krho->GetEntries(),
             nsig1_K1430pi = hsig_1K1430pi->GetEntries(), nsig2_K1430pi = hsig_2K1430pi->GetEntries(),
             nsig1_K892pi  = hsig_1K892pi->GetEntries(),  nsig2_K892pi  = hsig_2K892pi->GetEntries(),
             nsig1_Komega  = hsig_1Komega->GetEntries(),  nsig2_Komega  = hsig_2Komega->GetEntries();

     cout<<"nsig1: "<<nsig1<<" K rho: "<<nsig1_Krho<<" K(1430)pi: "<<nsig1_K1430pi<<" K(892)pi: "<<nsig1_K892pi<<" K omega: "<<nsig1_Komega<<endl;
     double aa=0.01, bb=0.50, cc=0.50, dd=0.99, ee=0.66, ff=0.98;

     double bottom_m = 0.14, left_m = 0.12;
     TCanvas *c = new TCanvas("c","", 0,0,800,250);
     //c->Range(0,0,1,1);
     c->SetBorderSize(1);
     c->SetLeftMargin(left_m);
     c->SetRightMargin(0.05);
     c->SetTopMargin(0.05);
     c->SetBottomMargin(bottom_m);
     c->SetFrameFillColor(0);

     TPad *c1_4 = new TPad("c1_4", "c1_4",aa,0.02,bb,0.98);
     c1_4->Draw();
     c1_4->cd();
     //c1_4->Range(0,0,1,1);
     c1_4->SetFillStyle(4000);
     c1_4->SetBorderMode(0);
     c1_4->SetBorderSize(1);
     c1_4->SetTickx();
     c1_4->SetTicky();
     c1_4->SetLeftMargin(left_m);
     c1_4->SetRightMargin(0.05);
     c1_4->SetTopMargin(0.05);
     c1_4->SetBottomMargin(bottom_m);
     c1_4->SetFrameFillColor(0);

     hbkg_1->Scale(1.0/10.8);
     //hsig_1Krho->Scale(Nsig1*nsig1_Krho/(nsig1*nsig1));
     //hsig_1K1430pi->Scale(Nsig1*nsig1_K1430pi/(nsig1*nsig1));                          
     //hsig_1K892pi->Scale(Nsig1*nsig1_K892pi/(nsig1*nsig1));
     //hsig_1Komega->Scale(Nsig1*nsig1_Komega/(nsig1*nsig1));
    
     hsig_1Krho->Scale(Nsig1/nsig1);
     hsig_1K1430pi->Scale(Nsig1/nsig1);                          
     hsig_1K892pi->Scale(Nsig1/nsig1);
     hsig_1Komega->Scale(Nsig1/nsig1);

     hsig_1K892pi->Add(hbkg_1);
     hsig_1Krho->Add(hsig_1K892pi);
     hsig_1K1430pi->Add(hsig_1Krho);
     hsig_1Komega->Add(hsig_1K1430pi);
     
      hsig_1Komega->Draw();
      hsig_1Komega->SetFillColor(4);
      hsig_1Komega->SetLineWidth(2);
      hsig_1Komega->SetMarkerStyle(20);
      hsig_1Komega->SetMarkerSize(0.4);
      hsig_1Komega->GetXaxis()->SetTitle("M_{K^{-}#pi^{+}}(GeV/#font[22]{c}^{2})");
      hsig_1Komega->GetXaxis()->SetTitleSize(0.07);
      hsig_1Komega->GetXaxis()->CenterTitle();
      hsig_1Komega->GetYaxis()->SetTitle("Events(/0.024 GeV/#font[22]{c}^{2})");
      hsig_1Komega->GetYaxis()->SetTitleSize(0.06);
      hsig_1Komega->GetYaxis()->CenterTitle();
      hsig_1Komega->GetYaxis()->SetRangeUser(0,1.4*hdata_1->GetMaximum());
      
      hsig_1K1430pi->Draw("same");
      hsig_1K1430pi->SetFillColor(8);
      hsig_1Krho->Draw("same");
      hsig_1Krho->SetFillColor(2);
      hsig_1K892pi->Draw("same");
      hsig_1K892pi->SetFillColor(30);

      hbkg_1->Draw("same");
      //hbkg_1->SetLineWidth(2);
      hbkg_1->SetFillColor(5);
     
      hdata_1->Draw("E same");
      hdata_1->SetLineWidth(2);
      hdata_1->SetMarkerStyle(20);
      hdata_1->SetMarkerSize(0.8);

      leg = new TLegend(0.6,0.68,0.9,0.87);
      leg->SetBorderSize(0);
      leg->SetTextFont(22);
      leg->SetLineColor(0);
      leg->SetFillColor(0);
      leg->AddEntry(hdata_1,"data");
      leg->AddEntry(hbkg_1,"inclusive MC","f");
      leg->AddEntry(hsig_1Krho,"K #rho","f");
      leg->AddEntry(hsig_1K1430pi,"K(1430) #pi","f");
      leg->AddEntry(hsig_1K892pi,"K(892) pi ","f");
      leg->AddEntry(hsig_1Komega,"K omega","f");

      leg->Draw();

      
      c1_4->Modified();
      c->cd(); 
   
      c1_5 = new TPad("c1_5", "c1_5",cc,0.02,dd,0.98);
      c1_5->Draw();
      c1_5->cd();
      c1_5->Range(0,0,1,1);
      c1_5->SetFillStyle(4000);
      c1_5->SetBorderMode(0);
      c1_5->SetBorderSize(2);
      c1_5->SetTickx();
      c1_5->SetTicky();
      c1_5->SetLeftMargin(left_m);
      c1_5->SetRightMargin(0.05);
      c1_5->SetTopMargin(0.05);
      c1_5->SetBottomMargin(bottom_m);
      c1_5->SetFrameFillColor(0);
  
     hbkg_2->Scale(1.0/10.8);
     hsig_2Krho->Scale(Nsig2/nsig2);
     hsig_2K1430pi->Scale(Nsig2/nsig2);                          
     hsig_2K892pi->Scale(Nsig2/nsig2);
     hsig_2Komega->Scale(Nsig2/nsig2);
     
     hsig_2Krho->Add(hbkg_2);
     hsig_2K892pi->Add(hsig_2Krho);
     hsig_2K1430pi->Add(hsig_2K892pi);
     hsig_2Komega->Add(hsig_2K1430pi);

      hsig_2Komega->Draw();
      hsig_2Komega->SetFillColor(4);
      hsig_2Komega->SetMarkerStyle(20);
      hsig_2Komega->SetMarkerSize(0.8);
      hsig_2Komega->GetXaxis()->SetTitle("M_{#pi^{+}#pi^{-}}(GeV/#font[22]{c}^{2})");
      hsig_2Komega->GetYaxis()->SetTitle("Events(/0.028 GeV/#font[22]{c}^{2})");
      hsig_2Komega->GetXaxis()->SetTitleSize(0.07);
      hsig_2Komega->GetYaxis()->SetTitleSize(0.06);
      hsig_2Komega->GetXaxis()->SetNdivisions(4);
      hsig_2Komega->GetYaxis()->SetNdivisions(6);
      hsig_2Komega->GetXaxis()->CenterTitle();
      hsig_2Komega->GetYaxis()->CenterTitle();
      hsig_2Komega->GetYaxis()->SetRangeUser(0,1.4*hdata_2->GetMaximum());

      hsig_2K1430pi->Draw("same");
      hsig_2K1430pi->SetFillColor(8);
      hsig_2K892pi->Draw("same");
      hsig_2K892pi->SetFillColor(30);
      hsig_2Krho->Draw("same");
      hsig_2Krho->SetFillColor(2);
      
      hbkg_2->Draw("same");
      hbkg_2->SetFillColor(5);
      hdata_2->Draw("E same");
      hdata_2->SetLineWidth(2);
      hdata_2->SetMarkerStyle(20);
      hdata_2->SetMarkerSize(0.8);
  
     leg = new TLegend(0.2,0.68,0.55,0.87);
     leg->SetBorderSize(0);
     leg->SetTextFont(22);
     leg->SetLineColor(0);
     leg->SetFillColor(0);
     leg->AddEntry(hdata_2,"data");
     leg->AddEntry(hbkg_2,"inclusive MC","f");
     leg->AddEntry(hsig_2Krho,"K #rho","f");
     leg->AddEntry(hsig_2K1430pi,"K(1430) #pi","f");
     leg->AddEntry(hsig_2K892pi,"K(892) pi ","f");
     leg->AddEntry(hsig_2Komega,"K omega","f");
     leg->Draw();


    c1_5->Modified();
    c->cd();
    c->SetSelected(c);
    c->Print("normal_mass.eps");
}
