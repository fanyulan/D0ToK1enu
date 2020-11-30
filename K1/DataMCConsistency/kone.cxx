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
void kone(){

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
bes3Style->SetPadBottomMargin(.12);
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
      tdata.Add("/scratchfs/bes/fangyl/data1/data/all/*.root");

      TChain tsig("tagD");
      tsig.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/BELLE/fit2/red_truth/truth*.root");
      TChain tbkg("tagD");
      tbkg.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/new/inclusiveMC/umiss2fit/reduce_truth/root/truth*.root");//inclusive mc
      TCut cut_signal="isTrueTag&&isSignal&&abs(KpTrueID)==321&&abs(pi1TrueID)==211&&abs(pi2TrueID)==211&&abs(elecTrueID)==11&&abs(elecMomID)==421&&(((mcmodea==666||mcmodea==222)&&mcmode1==103)||((mcmodeb==-666||mcmodeb==-222)&&mcmode2==-103))";
      TCut cuta = "isbest==1&&type==1&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)&&m_pipi>0.31&&fabs(m_pipi-0.497611)>0.01&&fabs(deltaE_pipzswp)>0.012";
      TCut cutb = "fabs(Umissfit)<0.1";
      TCut cutc = "elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge";
      TCut cut0 = "mode==0&&clveto==1&&deltaE>-0.029&&deltaE<0.027";
      TCut cut1 = "mode==1&&deltaE>-0.069&&deltaE<0.038";
      TCut cut3 = "mode==3&&deltaE>-0.031&&deltaE<0.028";
      TCut K_2  = "mcmode1==104||mcmode2==-104";
     
      
      int bins = 60;
      hsig_1        = new TH1F("hsig_1","", bins,0.9,1.68);
      hbkg_1        = new TH1F("hbkg_1","", bins,0.9,1.68);
      hdata_1       = new TH1F("hdata_1","",bins,0.9,1.68);
     
      tdata.Draw("mhad>>hdata_1", (cut0||cut1||cut3)&&cuta&&cutb&&cutc);
      tsig.Draw( "mhad>>hsig_1",  (cut0||cut1||cut3)&&cuta&&cutb&&cutc&&(cut_signal));  
      tbkg.Draw( "mhad>>hbkg_1",  (cut0||cut1||cut3)&&cuta&&cutb&&cutc&&!(cut_signal)&&!K_2);
     
      double ndata1= hdata_1->GetEntries();
      double nbkg = hbkg_1->GetEntries();
      double nsig = hsig_1->GetEntries();
      double Nbkg = nbkg/10.8;
      double Nsig = ndata1-Nbkg;

      TCanvas *c = new TCanvas("c","", 0,0, 800,600);
      c->Range(0,0,1,1);
      c->SetBorderSize(2);
      c->SetLeftMargin(0.05);
      c->SetRightMargin(0.05);
      c->SetTopMargin(0.05);
      c->SetBottomMargin(0.05);
      c->SetFrameFillColor(0);

     TPad *c1_1 = new TPad("c1_1", "c1_1",0.05,0.05,0.95,0.95);
     c1_1->Draw();
     c1_1->cd();
     c1_1->Range(0,0,1,1);
     c1_1->SetFillStyle(4000);
     c1_1->SetBorderMode(0);
     c1_1->SetBorderSize(0.3);
     c1_1->SetTickx();
     c1_1->SetTicky();
     c1_1->SetLeftMargin(0.12);
     c1_1->SetRightMargin(0.05);
     c1_1->SetTopMargin(0.1);
     c1_1->SetBottomMargin(0.12);
     c1_1->SetFrameFillColor(0);

     hbkg_1->Scale(1.0/10.8);     
     hbkg_1->Draw();
     hbkg_1->SetLineWidth(2);
     hbkg_1->SetFillColor(kYellow);
     hbkg_1->SetFillColor(5);
     hbkg_1->SetMarkerStyle(20);
     hbkg_1->SetMarkerSize(0.8);
     hbkg_1->GetXaxis()->SetTitle("M_{K^{-}#pi^{+}#pi^{-}} (GeV/c^{2})");
     hbkg_1->GetYaxis()->SetTitle("Events(/0.013 GeV/#font[22]{c})");
     hbkg_1->GetXaxis()->SetTitleSize(0.05);
     hbkg_1->GetYaxis()->SetTitleSize(0.05);
     hbkg_1->GetXaxis()->SetNdivisions(4);
     hbkg_1->GetYaxis()->SetNdivisions(6);
     hbkg_1->GetXaxis()->CenterTitle();
     hbkg_1->GetYaxis()->CenterTitle();
     hbkg_1->GetYaxis()->SetRangeUser(0,1.55*hdata_1->GetMaximum());
      hdata_1->Draw("E same");
      hdata_1->SetLineWidth(2);
      hdata_1->SetMarkerStyle(20);
      hdata_1->SetMarkerSize(0.8);

      hsig_1->Scale(Nsig/nsig);
      hsig_1->Add(hbkg_1,1);
      hsig_1->Draw("same");
      hsig_1->SetLineWidth(2);
      hsig_1->SetLineColor(2);
      hsig_1->SetMarkerStyle(20);
      hsig_1->SetMarkerSize(0.8);

      leg = new TLegend(0.6,0.7,0.95,0.89);
      leg->SetBorderSize(0);
      leg->SetTextFont(22);
      leg->SetLineColor(0);
      leg->SetFillColor(0);
      leg->AddEntry(hdata_1,"data");
      leg->AddEntry(hbkg_1,"inclusive MC","f");
      leg->AddEntry(hsig_1,"signal MC","f");
      leg->Draw();
      
      c1_1->Modified();
      c->cd(); 
      c->Print("./figure/kone.eps");
}
