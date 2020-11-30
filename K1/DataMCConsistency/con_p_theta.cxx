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
void con_p_theta(){

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
Double_t tsize=0.06;
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


using namespace RooFit;

      TChain tdata("tagD");  
      tdata.Add("/scratchfs/bes/fangyl/data1/data/all/*.root");
      TChain tsig("tagD");
      tsig.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/BELLE/fit2/red_truth/truth*.root");
      TChain tbkg("tagD");
      tbkg.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/new/inclusiveMC/umiss2fit/reduce_truth/root/truth*.root");

      TCut cut_signal="isTrueTag&&isSignal&&abs(KpTrueID)==321&&abs(pi1TrueID)==211&&abs(pi2TrueID)==211&&abs(elecTrueID)==11&&abs(elecMomID)==421&&(((mcmodea==666||mcmodea==222)&&mcmode1==103)||((mcmodeb==-666||mcmodeb==-222)&&mcmode2==-103))";
      TCut cuta = "isbest&&type&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)&&m_pipi>0.31&&fabs(m_pipi-0.497611)>0.01&&fabs(deltaE_pipzswp)>0.012";
      //TCut cutb = "fabs(Umissfit)<0.05&&mhad>1.1&&mhad<1.4";
      TCut cutb = "fabs(Umissfit)<0.025&&mhad>1.172&&mhad<1.372";
      TCut cutc = "elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge";
      TCut cut0 = "mode==0&&clveto&&deltaE>-0.029&&deltaE<0.027";
      TCut cut1 = "mode==1&&deltaE>-0.069&&deltaE<0.038";
      TCut cut3 = "mode==3&&deltaE>-0.031&&deltaE<0.028";
      TCut cut_K2 = "mcmode1==104||mcmode2==-104";

      int bins = 25;
      hsig_1        = new TH1F("hsig_1","", bins,0.,0.8);
      hbkg_1        = new TH1F("hbkg_1","", bins,0.,0.8);
      hdata_1       = new TH1F("hdata_1","",bins,0.,0.8);
      
      hsig_2        = new TH1F("hsig_2","", bins,0.,0.8);
      hbkg_2        = new TH1F("hbkg_2","", bins,0.,0.8);
      hdata_2       = new TH1F("hdata_2","",bins,0.,0.8);

      hsig_3        = new TH1F("hsig_3","", bins,0.,0.8);
      hbkg_3        = new TH1F("hbkg_3","", bins,0.,0.8);
      hdata_3       = new TH1F("hdata_3","",bins,0.,0.8);
 
      hsig_4        = new TH1F("hsig_4","", bins,0.,0.8);
      hbkg_4        = new TH1F("hbkg_4","", bins,0.,0.8);
      hdata_4       = new TH1F("hdata_4","",bins,0.,0.8);    
      
      hsig_5        = new TH1F("hsig_5","", bins,-1.0,1.0);
      hbkg_5        = new TH1F("hbkg_5","", bins,-1.0,1.0);
      hdata_5       = new TH1F("hdata_5","",bins,-1.0,1.0);

      hsig_6        = new TH1F("hsig_6","", bins,-1.0,1.0);
      hbkg_6        = new TH1F("hbkg_6","", bins,-1.0,1.0);
      hdata_6       = new TH1F("hdata_6","",bins,-1.0,1.0);  

      hsig_7        = new TH1F("hsig_7","", bins,-1.0,1.0);
      hbkg_7        = new TH1F("hbkg_7","", bins,-1.0,1.0);
      hdata_7       = new TH1F("hdata_7","",bins,-1.0,1.0);

      hsig_8        = new TH1F("hsig_8","", bins,-1.0,1.0);
      hbkg_8        = new TH1F("hbkg_8","", bins,-1.0,1.0);
      hdata_8       = new TH1F("hdata_8","",bins,-1.0,1.0);
 
      tdata.Draw("p_kaonp3>>hdata_1",(cut0||cut1||cut3)&&cuta&&cutb&&cutc);
      tsig.Draw( "p_kaonp3>>hsig_1"  ,(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal);  
      tbkg.Draw( "p_kaonp3>>hbkg_1",(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&!cut_signal&&!cut_K2);
 
      tdata.Draw("p_pionp3>>hdata_2",(cut0||cut1||cut3)&&cuta&&cutb&&cutc);
      tsig.Draw( "p_pionp3>>hsig_2"  ,(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal);
      tbkg.Draw( "p_pionp3>>hbkg_2",(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&!cut_signal&&!cut_K2);
     
      tdata.Draw("p_pion2p3>>hdata_3",(cut0||cut1||cut3)&&cuta&&cutb&&cutc);
      tsig.Draw( "p_pion2p3>>hsig_3"  ,(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal);
      tbkg.Draw( "p_pion2p3>>hbkg_3",(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&!cut_signal&&!cut_K2);
      
      tdata.Draw("elecP>>hdata_4",(cut0||cut1||cut3)&&cuta&&cutb&&cutc);
      tsig.Draw( "elecP>>hsig_4"  ,(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal);
      tbkg.Draw( "elecP>>hbkg_4",(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&!cut_signal&&!cut_K2);
      
      tdata.Draw("cosK>>hdata_5",(cut0||cut1||cut3)&&cuta&&cutb&&cutc);
      tsig.Draw( "cosK>>hsig_5"  ,(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal);
      tbkg.Draw( "cosK>>hbkg_5",(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&!cut_signal&&!cut_K2);

      tdata.Draw("cosPi1>>hdata_6",(cut0||cut1||cut3)&&cuta&&cutb&&cutc);
      tsig.Draw( "cosPi1>>hsig_6"  ,(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal);
      tbkg.Draw( "cosPi1>>hbkg_6",(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&!cut_signal&&!cut_K2);

      tdata.Draw("cosPi2>>hdata_7",(cut0||cut1||cut3)&&cuta&&cutb&&cutc);
      tsig.Draw( "cosPi2>>hsig_7"  ,(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal);
      tbkg.Draw( "cosPi2>>hbkg_7",(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&!cut_signal&&!cut_K2);

      tdata.Draw("cosE>>hdata_8",(cut0||cut1||cut3)&&cuta&&cutb&&cutc);
      tsig.Draw( "cosE>>hsig_8"  ,(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&cut_signal);
      tbkg.Draw( "cosE>>hbkg_8",(cut0||cut1||cut3)&&cuta&&cutb&&cutc&&!cut_signal&&!cut_K2);

      double s_bkg = 10.8;
      double nbkg1 = hbkg_1->GetEntries(),nbkg2 = hbkg_2->GetEntries(),nbkg3 = hbkg_3->GetEntries(),nbkg4 = hbkg_4->GetEntries();
      double nbkg5 = hbkg_5->GetEntries(),nbkg6 = hbkg_6->GetEntries(),nbkg7 = hbkg_7->GetEntries(),nbkg8 = hbkg_8->GetEntries();
      double ndata1= hdata_1->GetEntries(),ndata2= hdata_2->GetEntries(),ndata3= hdata_3->GetEntries(),ndata4= hdata_4->GetEntries();
      double ndata5= hdata_5->GetEntries(),ndata6= hdata_6->GetEntries(),ndata7= hdata_7->GetEntries(),ndata8= hdata_8->GetEntries();
      double nsig1 = hsig_1->GetEntries(),nsig2 = hsig_2->GetEntries(),nsig3 = hsig_3->GetEntries(),nsig4 = hsig_4->GetEntries();
      double nsig5 = hsig_5->GetEntries(),nsig6 = hsig_6->GetEntries(),nsig7 = hsig_7->GetEntries(),nsig8 = hsig_8->GetEntries();
      double Nsig1 = ndata1-nbkg1/10.8, Nsig2 = ndata2-nbkg2/10.8,Nsig3 = ndata3-nbkg3/10.8,Nsig4 = ndata4-nbkg4/10.8;
      double Nsig5 = ndata5-nbkg5/10.8, Nsig6 = ndata6-nbkg6/10.8,Nsig7 = ndata7-nbkg7/10.8,Nsig8 = ndata8-nbkg8/10.8;
      double aa=0.01, bb=0.25, cc=0.25, dd=0.49, ee=0.49, ff=0.74, gg=0.74, hh=0.98;
      double bottom_m = 0.15,  left_m = 0.1,  title_s = 0.07;
      TCanvas *c = new TCanvas("c","", 0,0,1200,600);
      c->Range(0,0,1,1);
      c->SetBorderSize(2);
      c->SetLeftMargin(left_m);
      c->SetRightMargin(0.05);
      c->SetTopMargin(0.05);
      c->SetBottomMargin(bottom_m);
      c->SetFrameFillColor(0);
      
      TPad *c1_1 = new TPad("c1_1", "c1_1",aa,0.51,bb,0.99);
      c1_1->Draw();
      c1_1->cd();
      c1_1->Range(0,0,1,1);
      c1_1->SetFillStyle(4000);
      c1_1->SetBorderMode(0);
      c1_1->SetBorderSize(2);
      c1_1->SetTickx();
      c1_1->SetTicky();
      c1_1->SetLeftMargin(0.15);
      c1_1->SetRightMargin(0.05);
      c1_1->SetTopMargin(0.02);
      c1_1->SetBottomMargin(bottom_m);
      c1_1->SetFrameFillColor(0);
 
     hbkg_1->Scale(1.0/10.8);
     hbkg_1->Draw();
     hbkg_1->SetLineWidth(2);
     hbkg_1->SetLineColor(1);
     hbkg_1->SetFillColor(5);
     hbkg_1->SetMarkerStyle(20);
     hbkg_1->SetMarkerSize(0.8);
     hbkg_1->GetXaxis()->SetTitle("p_{K^{-}}(GeV/#font[22]{c})");
     hbkg_1->GetYaxis()->SetTitle("Events(/0.032 GeV/#font[22]{c})");
     hbkg_1->GetXaxis()->SetTitleSize(title_s);
     hbkg_1->GetYaxis()->SetTitleSize(title_s);
     hbkg_1->GetXaxis()->SetNdivisions(4);
     hbkg_1->GetYaxis()->SetNdivisions(6);
     hbkg_1->GetXaxis()->CenterTitle();
     hbkg_1->GetYaxis()->CenterTitle();
     hbkg_1->GetYaxis()->SetRangeUser(0,1.55*hdata_1->GetMaximum());
      hdata_1->Draw("E same");
      hdata_1->SetLineWidth(2);
      hdata_1->SetMarkerStyle(20);
      hdata_1->SetMarkerSize(0.8);

      hsig_1->Scale(Nsig1/nsig1);
      hsig_1->Add(hbkg_1,1);
      hsig_1->Draw("same");
      hsig_1->SetLineWidth(2);
      hsig_1->SetLineColor(2);      
      hsig_1->SetMarkerStyle(20);
      hsig_1->SetMarkerSize(0.8);
      double legx0=0.6, legy0=0.7;
      leg = new TLegend(legx0,legy0,0.91,0.89);
      leg->SetBorderSize(0);
      leg->SetTextFont(22);
      leg->SetLineColor(0);
      leg->SetFillColor(0);
      //leg->SetHeader("inclusive MC");
      leg->AddEntry(hdata_1,"data");
      leg->AddEntry(hbkg_1,"inclusive MC","f");
      leg->AddEntry(hsig_1,"signal MC","f");
      leg->Draw();
      
      c1_1->Modified();
      c->cd(); 
   
      c1_2 = new TPad("c1_2", "c1_2",cc,0.51,dd,0.99);
      c1_2->Draw();
      c1_2->cd();
      c1_2->Range(0,0,1,1);
      c1_2->SetFillStyle(4000);
      c1_2->SetBorderMode(0);
      c1_2->SetBorderSize(2);
      c1_2->SetTickx();
      c1_2->SetTicky();
      c1_2->SetLeftMargin(left_m);
      c1_2->SetRightMargin(0.05);
      c1_2->SetTopMargin(0.05);
      c1_2->SetBottomMargin(bottom_m);
      c1_2->SetFrameFillColor(0);

    hbkg_2->Scale(1.0/10.8);
    hbkg_2->Draw();
    hbkg_2->SetLineWidth(2);
    hbkg_2->SetLineColor(1);
    hbkg_2->SetFillColor(5);
    hbkg_2->SetMarkerStyle(20);
    hbkg_2->SetMarkerSize(0.8);
    hbkg_2->GetXaxis()->SetTitle("p_{#pi^{+}}(GeV/#font[22]{c})");
    //hbkg_2->GetYaxis()->SetTitle("Events(/0.032 GeV/#font[22]{c})");
    hbkg_2->GetXaxis()->SetTitleSize(title_s);
    hbkg_2->GetYaxis()->SetTitleSize(title_s);
    hbkg_2->GetXaxis()->SetNdivisions(4);
    hbkg_2->GetYaxis()->SetNdivisions(6);
    hbkg_2->GetXaxis()->CenterTitle();
    hbkg_2->GetYaxis()->CenterTitle();
    hbkg_2->GetYaxis()->SetRangeUser(0,1.55*hdata_2->GetMaximum());
    hdata_2->Draw("E same");
    hdata_2->SetLineWidth(2);
    hdata_2->SetMarkerStyle(20);
    hdata_2->SetMarkerSize(0.8);
 
    hsig_2->Scale(Nsig2/nsig2);     
    hsig_2->Add(hbkg_2,1);
    hsig_2->Draw("same");
    hsig_2->SetLineWidth(2);
    hsig_2->SetLineColor(2);
    hsig_2->SetMarkerStyle(20);
    hsig_2->SetMarkerSize(0.8);
    leg = new TLegend(legx0,legy0,0.91,0.89);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    //leg->SetHeader("inclusive MC");
    leg->AddEntry(hdata_2,"data");
    leg->AddEntry(hbkg_2,"inclusive MC","f");
    leg->AddEntry(hsig_2,"signal MC","f");
    leg->Draw();

    c1_2->Modified();
    c->cd();

    c1_3 = new TPad("c1_3", "c1_3",ee,0.51,ff,0.99);
    c1_3->Draw();
    c1_3->cd();
    c1_3->Range(0,0,1,1);
    c1_3->SetFillStyle(4000);
    c1_3->SetBorderMode(0);
    c1_3->SetBorderSize(2);
    c1_3->SetTickx();
    c1_3->SetTicky();
    c1_3->SetLeftMargin(left_m);
    c1_3->SetRightMargin(0.05);
    c1_3->SetTopMargin(0.05);
    c1_3->SetBottomMargin(bottom_m);
    c1_3->SetFrameFillColor(0);

   hbkg_3->Scale(1.0/10.8);
   hbkg_3->Draw();
   hbkg_3->SetLineWidth(2);
   hbkg_3->SetLineColor(1);
   hbkg_3->SetFillColor(5);
   hbkg_3->SetMarkerStyle(20);
   hbkg_3->SetMarkerSize(0.8);
   hbkg_3->GetXaxis()->SetTitle("p_{#pi^{-}}(GeV/#font[22]{c})");
   //hbkg_3->GetYaxis()->SetTitle("Events(/0.032 GeV/#font[22]{c})");
   hbkg_3->GetXaxis()->SetTitleSize(title_s);
   hbkg_3->GetYaxis()->SetTitleSize(title_s);
   hbkg_3->GetXaxis()->SetNdivisions(4);
   hbkg_3->GetYaxis()->SetNdivisions(6);
   hbkg_3->GetXaxis()->CenterTitle();
   hbkg_3->GetYaxis()->CenterTitle();
   hbkg_3->GetYaxis()->SetRangeUser(0,1.55*hdata_3->GetMaximum());
   hdata_3->Draw("E same");
   hdata_3->SetLineWidth(2);
   hdata_3->SetMarkerStyle(20);
   hdata_3->SetMarkerSize(0.8);
   
    hsig_3->Scale(Nsig3/nsig3);
    hsig_3->Add(hbkg_3,1);
    hsig_3->Draw("same");
    hsig_3->SetLineWidth(2);
    hsig_3->SetLineColor(2);
    hsig_3->SetMarkerStyle(20);
    hsig_3->SetMarkerSize(0.8);
    leg = new TLegend(legx0,legy0,0.91,0.89);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    //leg->SetHeader("inclusive MC");
    leg->AddEntry(hdata_3,"data");
    leg->AddEntry(hbkg_3,"inclusive MC","f");
    leg->AddEntry(hsig_3,"signal MC","f");
    leg->Draw();


    c1_3->Modified();
    c->cd();

     c1_4 = new TPad("c1_4", "c1_4",gg,0.51,hh,0.99);
     c1_4->Draw();
     c1_4->cd();
     c1_4->Range(0,0,1,1);
     c1_4->SetFillStyle(4000);
     c1_4->SetBorderMode(0);
     c1_4->SetBorderSize(2);
     c1_4->SetTickx();
     c1_4->SetTicky();
     c1_4->SetLeftMargin(left_m);
     c1_4->SetRightMargin(0.05);
     c1_4->SetTopMargin(0.05);
     c1_4->SetBottomMargin(bottom_m);
     c1_4->SetFrameFillColor(0);

   hbkg_4->Scale(1.0/10.8);
   hbkg_4->Draw();
   hbkg_4->SetLineWidth(2);
   hbkg_4->SetLineColor(1);
   hbkg_4->SetFillColor(5);
   hbkg_4->SetMarkerStyle(20);
   hbkg_4->SetMarkerSize(0.8);
   hbkg_4->GetXaxis()->SetTitle("p_{e^{+}}(GeV/#font[22]{c})");
   //hbkg_4->GetYaxis()->SetTitle("Events(/0.032 GeV/#font[22]{c})");
   hbkg_4->GetXaxis()->SetTitleSize(title_s);
   hbkg_4->GetYaxis()->SetTitleSize(title_s);
   hbkg_4->GetXaxis()->SetNdivisions(4);
   hbkg_4->GetYaxis()->SetNdivisions(6);
   hbkg_4->GetXaxis()->CenterTitle();
   hbkg_4->GetYaxis()->CenterTitle();
   hbkg_4->GetYaxis()->SetRangeUser(0,1.55*hdata_4->GetMaximum());
   hdata_4->Draw("E same");
   hdata_4->SetLineWidth(2);
   hdata_4->SetMarkerStyle(20);
   hdata_4->SetMarkerSize(0.8);

    hsig_4->Scale(Nsig4/nsig4);
    hsig_4->Add(hbkg_4,1);
    hsig_4->Draw("same");
    hsig_4->SetLineWidth(2);
    hsig_4->SetLineColor(2);
    hsig_4->SetMarkerStyle(20);
    hsig_4->SetMarkerSize(0.8);
    leg = new TLegend(legx0,legy0,0.91,0.89);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    //leg->SetHeader("inclusive MC");
    leg->AddEntry(hdata_4,"data");
    leg->AddEntry(hbkg_4,"inclusive MC","f");
    leg->AddEntry(hsig_4,"signal MC","f");
    leg->Draw();

    
    c1_4->Modified();
    c->cd();


    c1_5 = new TPad("c1_5", "c1_5",aa,0.01,bb,0.50);
    c1_5->Draw();
    c1_5->cd();
    c1_5->Range(0,0,1,1);
    c1_5->SetFillStyle(4000);
    c1_5->SetBorderMode(0);
    c1_5->SetBorderSize(2);
    c1_5->SetTickx();
    c1_5->SetTicky();
    c1_5->SetLeftMargin(0.15);
    c1_5->SetRightMargin(0.05);
    c1_5->SetTopMargin(0.05);
    c1_5->SetBottomMargin(bottom_m);
    c1_5->SetFrameFillColor(0);

    hbkg_5->Scale(1.0/10.8);
    hbkg_5->Draw();
    hbkg_5->SetLineWidth(2);
    hbkg_5->SetLineColor(1);
    hbkg_5->SetFillColor(5);
    hbkg_5->SetMarkerStyle(20);
    hbkg_5->SetMarkerSize(0.8);
    hbkg_5->GetXaxis()->SetTitle("cos(K^{-})");
    hbkg_5->GetYaxis()->SetTitle("Events(/0.08)");
    hbkg_5->GetXaxis()->SetTitleSize(title_s);
    hbkg_5->GetYaxis()->SetTitleSize(title_s);
    hbkg_5->GetXaxis()->SetNdivisions(4);
    hbkg_5->GetYaxis()->SetNdivisions(6);
    hbkg_5->GetXaxis()->CenterTitle();
    hbkg_5->GetYaxis()->CenterTitle();
    hbkg_5->GetYaxis()->SetRangeUser(0,1.55*hdata_5->GetMaximum());

   hdata_5->Draw("E same");
   hdata_5->SetLineWidth(2);
   hdata_5->SetMarkerStyle(20);
   hdata_5->SetMarkerSize(0.8);

    hsig_5->Scale(Nsig5/nsig5);
    hsig_5->Add(hbkg_5,1);
    hsig_5->Draw("same");
    hsig_5->SetLineWidth(2);
    hsig_5->SetLineColor(2);
    hsig_5->SetMarkerStyle(20);
    hsig_5->SetMarkerSize(0.8);

    leg = new TLegend(legx0,legy0,0.91,0.89);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    //leg->SetHeader("inclusive MC");
    leg->AddEntry(hdata_5,"data");
    leg->AddEntry(hbkg_5,"inclusive MC","f");
    leg->AddEntry(hsig_5,"signal MC","f");
    leg->Draw();


    c1_5->Modified();
    c->cd();

   c1_6 = new TPad("c1_6", "c1_6",cc,0.01,dd,0.50);
   c1_6->Draw();
   c1_6->cd();
   c1_6->Range(0,0,1,1);
   c1_6->SetFillStyle(4000);
   c1_6->SetBorderMode(0);
   c1_6->SetBorderSize(2);
   c1_6->SetTickx();
   c1_6->SetTicky();
   c1_6->SetLeftMargin(left_m);
   c1_6->SetRightMargin(0.05);
   c1_6->SetTopMargin(0.05);
   c1_6->SetBottomMargin(bottom_m);
   c1_6->SetFrameFillColor(0);
   hbkg_6->Scale(1.0/10.8);
   hbkg_6->Draw();
   hbkg_6->SetLineWidth(2);
   hbkg_6->SetLineColor(1);
   hbkg_6->SetFillColor(5);
   hbkg_6->SetMarkerStyle(20);
   hbkg_6->SetMarkerSize(0.8);
   hbkg_6->GetXaxis()->SetTitle("cos(#pi^{+})");
   //hbkg_6->GetYaxis()->SetTitle("Events/(0.08)");
   hbkg_6->GetXaxis()->SetTitleSize(title_s);
   hbkg_6->GetYaxis()->SetTitleSize(title_s);
   hbkg_6->GetXaxis()->SetNdivisions(4);
   hbkg_6->GetYaxis()->SetNdivisions(6);
   hbkg_6->GetXaxis()->CenterTitle();
   hbkg_6->GetYaxis()->CenterTitle();
   hbkg_6->GetYaxis()->SetRangeUser(0,1.55*hdata_6->GetMaximum());

   hdata_6->Draw("E same");
   hdata_6->SetLineWidth(2);
   hdata_6->SetMarkerStyle(20);
   hdata_6->SetMarkerSize(0.8);

    hsig_6->Scale(Nsig6/nsig6);
    hsig_6->Add(hbkg_6,1);
    hsig_6->Draw("same");
    hsig_6->SetLineWidth(2);
    hsig_6->SetLineColor(2);
    hsig_6->SetMarkerStyle(20);
    hsig_6->SetMarkerSize(0.8);
    leg = new TLegend(legx0,legy0,0.91,0.89);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    //leg->SetHeader("inclusive MC");
    leg->AddEntry(hdata_6,"data");
    leg->AddEntry(hbkg_6,"inclusive MC","f");
    leg->AddEntry(hsig_6,"signal MC","f");
    leg->Draw();


    c1_6->Modified();
    c->cd();


     c1_7 = new TPad("c1_7", "c1_7",ee,0.01,ff,0.50);
     c1_7->Draw();
     c1_7->cd();
     c1_7->Range(0,0,1,1);
     c1_7->SetFillStyle(4000);
     c1_7->SetBorderMode(0);
     c1_7->SetBorderSize(2);
     c1_7->SetTickx();
     c1_7->SetTicky();
     c1_7->SetLeftMargin(left_m);
     c1_7->SetRightMargin(0.05);
     c1_7->SetTopMargin(0.05);
     c1_7->SetBottomMargin(bottom_m);
     c1_7->SetFrameFillColor(0);
     hbkg_7->Scale(1.0/10.8);
     hbkg_7->Draw();
     hbkg_7->SetLineWidth(2);
     hbkg_7->SetLineColor(1);
     hbkg_7->SetFillColor(5);
     hbkg_7->SetMarkerStyle(20);
     hbkg_7->SetMarkerSize(0.8);
     hbkg_7->GetXaxis()->SetTitle("cos(#pi^{-})");
     //hbkg_7->GetYaxis()->SetTitle("Events/(0.08)");
     hbkg_7->GetXaxis()->SetTitleSize(title_s);
     hbkg_7->GetYaxis()->SetTitleSize(title_s);
     hbkg_7->GetXaxis()->SetNdivisions(4);
     hbkg_7->GetYaxis()->SetNdivisions(6);
     hbkg_7->GetXaxis()->CenterTitle();
     hbkg_7->GetYaxis()->CenterTitle();
     hbkg_7->GetYaxis()->SetRangeUser(0,1.55*hdata_7->GetMaximum());

   hdata_7->Draw("E same");
   hdata_7->SetLineWidth(2);
   hdata_7->SetMarkerStyle(20);
   hdata_7->SetMarkerSize(0.8);

    hsig_7->Scale(Nsig7/nsig7);
    hsig_7->Add(hbkg_7,1);
    hsig_7->Draw("same");
    hsig_7->SetLineWidth(2);
    hsig_7->SetLineColor(2);
    hsig_7->SetMarkerStyle(20);
    hsig_7->SetMarkerSize(0.8);
    leg = new TLegend(legx0,legy0,0.91,0.89);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    //leg->SetHeader("inclusive MC");
    leg->AddEntry(hdata_7,"data");
    leg->AddEntry(hbkg_7,"inclusive MC","f");
    leg->AddEntry(hsig_7,"signal MC","f");
    leg->Draw();


    c1_7->Modified();
    c->cd();

     c1_8 = new TPad("c1_8", "c1_8",gg,0.01,hh,0.50);
     c1_8->Draw();
     c1_8->cd();
     c1_8->Range(0,0,1,1);
     c1_8->SetFillStyle(4000);
     c1_8->SetBorderMode(0);
     c1_8->SetBorderSize(2);
     c1_8->SetTickx();
     c1_8->SetTicky();
     c1_8->SetLeftMargin(left_m);
     c1_8->SetRightMargin(0.05);
     c1_8->SetTopMargin(0.05);
     c1_8->SetBottomMargin(bottom_m);
     c1_8->SetFrameFillColor(0);
     hbkg_8->Scale(1.0/10.8);
     hbkg_8->Draw();
     hbkg_8->SetLineWidth(2);
     hbkg_8->SetLineColor(1);
     hbkg_8->SetFillColor(5);
     hbkg_8->SetMarkerStyle(20);
     hbkg_8->SetMarkerSize(0.8);
     hbkg_8->GetXaxis()->SetTitle("cos(e^{+})");
     //hbkg_8->GetYaxis()->SetTitle("Events(/0.08)");
     hbkg_8->GetXaxis()->SetTitleSize(title_s);
     hbkg_8->GetYaxis()->SetTitleSize(title_s);
     hbkg_8->GetXaxis()->SetNdivisions(4);
     hbkg_8->GetYaxis()->SetNdivisions(6);
     hbkg_8->GetXaxis()->CenterTitle();
     hbkg_8->GetYaxis()->CenterTitle();
     hbkg_8->GetYaxis()->SetRangeUser(0,1.55*hdata_8->GetMaximum());

   hdata_8->Draw("E same");
   hdata_8->SetLineWidth(2);
   hdata_8->SetMarkerStyle(20);
   hdata_8->SetMarkerSize(0.8);

    hsig_8->Scale(Nsig8/nsig8);
    hsig_8->Add(hbkg_8,1);
    hsig_8->Draw("same");
    hsig_8->SetLineWidth(2);
    hsig_8->SetLineColor(2);
    hsig_8->SetMarkerStyle(20);
    hsig_8->SetMarkerSize(0.8);
    leg = new TLegend(legx0,legy0,0.91,0.89);
    leg->SetBorderSize(0);
    leg->SetTextFont(22);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    //leg->SetHeader("inclusive MC");
    leg->AddEntry(hdata_8,"data");
    leg->AddEntry(hbkg_8,"inclusive MC","f");
    leg->AddEntry(hsig_8,"signal MC","f");
    leg->Draw();

    c1_8->Modified();
    c->cd();



    c->Modified();
    c->cd();
    c->SetSelected(c);

   c->Print("./figure/p_theta.eps");
}
