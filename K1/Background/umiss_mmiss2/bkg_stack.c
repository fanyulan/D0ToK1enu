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

// background analysis in inclsive MC ~ stack 
//nvar0 ~ before ; nvar8 ~ after
using namespace RooFit;
void bkg_stack()
{
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
    
for(int nvar=7; nvar<8; nvar++) {

    TCut cut_sig = "isTrueTag&&isSignal&&abs(KpTrueID)==321&&abs(pi1TrueID)==211&&abs(pi2TrueID)==211&&abs(elecTrueID)==11&&abs(elecMomID)==421";
    TCut cuta = "isbest==1&&type==1&&mBC<1.874&&mBC>1.858";
    TCut cutc = "elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge";
    TCut cut0 = "mode==0&&clveto==1&&deltaE>-0.029&&deltaE<0.027";
    TCut cut1 = "mode==1&&deltaE>-0.069&&deltaE<0.038";
    TCut cut3 = "mode==3&&deltaE>-0.031&&deltaE<0.028";
    TCut cut_K1 = "((mcmodea==666||mcmodea==222)&&mcmode1==103)||((mcmodeb==-666||mcmodeb==-222)&&mcmode2==-103)";
    //TCut cut_K2 = "((mcmodea==666||mcmodea==222)&&mcmode1==104)||((mcmodeb==-666||mcmodeb==-222)&&mcmode2==-104)";
    TCut cut_K2 = "mcmode1==104 || mcmode2==-104";
    
    TCut cut_before          = "abs(Umissfit)<0.2";
    TCut cut_mD_EasPi        = "abs(Umissfit)<0.2&&mD_EasPi<1.81";
    TCut cut_costh_lep_pi    = "abs(Umissfit)<0.2&&mD_EasPi<1.81&&costh_lep_pi<0.94";
    TCut cut_mD_EasPiwPi0    = "abs(Umissfit)<0.2&&mD_EasPi<1.81&&costh_lep_pi<0.94&&mD_EasPiwPi0<1.4";
    TCut cut_costh_miss_neut = "abs(Umissfit)<0.2&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81";
    TCut cut_chi2_eop        = "abs(Umissfit)<0.2&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)";
    TCut cut_m_pipi          = "abs(Umissfit)<0.2&&mD_EasPi<1.8&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)&&m_pipi>0.31&&abs(m_pipi-0.497611)>0.01";
    TCut cut_deltaE_pipzswp  = "abs(Umissfit)<0.2&&mD_EasPi<1.8&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)&&m_pipi>0.31&&abs(m_pipi-0.497611)>0.01&&fabs(deltaE_pipzswp)>0.012";
    
    TCut cut_455 = "(charm==1&&mcmodeb==-455) || (charm==-1&&mcmodea==455)";//kpiw
    TCut cut_103 = "(charm==1&&mcmodeb==-103) || (charm==-1&&mcmodea==103)";//kpipipi
    TCut cut_105 = "(charm==1&&mcmodeb==-105) || (charm==-1&&mcmodea==105)";//kpipipipi0
    TCut cut_715 = "(charm==1&&mcmodeb==-715) || (charm==-1&&mcmodea==715)";//kpienu
    TCut cut_201 = "(charm==1&&mcmodeb==-201) || (charm==-1&&mcmodea==201)";//kpi0enu
    TCut cut_104 = "(charm==1&&mcmodeb==-104) || (charm==-1&&mcmodea==104)";//kpipi0pi0
    TCut cut_102 = "(charm==1&&mcmodeb==-102) || (charm==-1&&mcmodea==102)";//kpipi0
    TCut cut_1032 = "(charm==1&&(mcmodeb==-1032||mcmodeb==-108)) || (charm==-1&&(mcmodea==1032||mcmodea==108))";//kpiw, w->pipi
    double bins = 80, xmin = -0.2, xmax = 0.2;
    TCanvas *c = new TCanvas("c", "umiss2fit_bkg",0,0, 800,600);
    sig        = new TH1F("sig", "",       bins, xmin, xmax);//signal
    //kpiw       = new TH1F("kpiw", "",      bins, xmin, xmax);//455
    //kpiw2      = new TH1F("kpiw2", "",     bins, xmin, xmax);//1032
    kpipipi    = new TH1F("kpipipi", "",   bins, xmin, xmax);//103
    kpipipipi0 = new TH1F("kpipipipi0", "",bins, xmin, xmax);//105
    kpienu     = new TH1F("kpienu", "",    bins, xmin, xmax);//715
    kpi0enu    = new TH1F("kpi0enu", "",   bins, xmin, xmax);//201
    kpipi0pi0  = new TH1F("kpipi0pi0", "", bins, xmin, xmax);//104
    kpipi0     = new TH1F("kpipi0", "",    bins, xmin, xmax);//102
    other      = new TH1F("other", "",     bins, xmin, xmax);//other else_all
  
    TChain t("tagD");
    t.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/new/inclusiveMC/umiss2fit/reduce_truth/root/truth*.root");
        
////////////////////////////before//////////////////////////////////////
    if(nvar==0) { 
        t.Draw("Umissfit>>sig",       (cut0||cut1||cut3)&&cuta&&cutc&&(cut_sig&&cut_K1)&&cut_before);
        //t.Draw("Umissfit>>kpiw",      (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_455&&cut_before);
        //t.Draw("Umissfit>>kpiw2",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_1032&&cut_before); 
        t.Draw("Umissfit>>kpipipi",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_103||cut_1032)&&cut_before);
        t.Draw("Umissfit>>kpipipipi0",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_105||cut_455)&&cut_before);
        t.Draw("Umissfit>>kpienu",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_715&&cut_before);
        t.Draw("Umissfit>>kpi0enu",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_201&&cut_before);
        t.Draw("Umissfit>>kpipi0pi0", (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_104&&cut_before);
        t.Draw("Umissfit>>kpipi0",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_102&&cut_before);
        t.Draw("Umissfit>>other",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_before&&!(cut_455||cut_103||cut_105||cut_715||cut_201||cut_104||cut_1032||cut_102));  }
        
////////////////////////////before/////////////////////////////////////

////////////////////////////cut//////////////////////////////////////

    if(nvar==1) { 
        t.Draw("Umissfit>>sig",       (cut0||cut1||cut3)&&cuta&&cutc&&cut_sig&&cut_K1&&cut_mD_EasPi);
        //t.Draw("Umissfit>>kpiw",      (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_455&&cut_mD_EasPi);
        //t.Draw("Umissfit>>kpiw2",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_1032&&cut_mD_EasPi);
        t.Draw("Umissfit>>kpipipi",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_103||cut_1032)&&cut_mD_EasPi);
        t.Draw("Umissfit>>kpipipipi0",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_105||cut_455)&&cut_mD_EasPi);
        t.Draw("Umissfit>>kpienu",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_715&&cut_mD_EasPi);
        t.Draw("Umissfit>>kpi0enu",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_201&&cut_mD_EasPi);
        t.Draw("Umissfit>>kpipi0pi0", (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_104&&cut_mD_EasPi);
        t.Draw("Umissfit>>kpipi0",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_102&&cut_mD_EasPi);
        t.Draw("Umissfit>>other",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_mD_EasPi&&!(cut_455||cut_103||cut_105||cut_715||cut_201||cut_104||cut_1032||cut_102));  }
        
    if(nvar==3) {
        t.Draw("Umissfit>>sig",       (cut0||cut1||cut3)&&cuta&&cutc&&cut_sig&&cut_K1&&cut_mD_EasPiwPi0);
        //t.Draw("Umissfit>>kpiw",      (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_455&&cut_mD_EasPiwPi0);
        //t.Draw("Umissfit>>kpiw2",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_1032&&cut_mD_EasPiwPi0); 
        t.Draw("Umissfit>>kpipipi",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_103||cut_1032)&&cut_mD_EasPiwPi0);
        t.Draw("Umissfit>>kpipipipi0",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_105||cut_455)&&cut_mD_EasPiwPi0);
        t.Draw("Umissfit>>kpienu",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_715&&cut_mD_EasPiwPi0);
        t.Draw("Umissfit>>kpi0enu",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_201&&cut_mD_EasPiwPi0);
        t.Draw("Umissfit>>kpipi0pi0", (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_104&&cut_mD_EasPiwPi0);
        t.Draw("Umissfit>>kpipi0",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_102&&cut_mD_EasPiwPi0);
        t.Draw("Umissfit>>other",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_mD_EasPiwPi0&&!(cut_455||cut_103||cut_105||cut_715||cut_201||cut_104||cut_1032||cut_102));  }
    if(nvar==2) {
        t.Draw("Umissfit>>sig",       (cut0||cut1||cut3)&&cuta&&cutc&&cut_sig&&cut_K1&&cut_costh_lep_pi);
        //t.Draw("Umissfit>>kpiw",      (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_455&&cut_costh_lep_pi);
        //t.Draw("Umissfit>>kpiw2",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_1032&&cut_costh_lep_pi);
        t.Draw("Umissfit>>kpipipi",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_103||cut_1032)&&cut_costh_lep_pi);
        t.Draw("Umissfit>>kpipipipi0",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_105||cut_455)&&cut_costh_lep_pi);
        t.Draw("Umissfit>>kpienu",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_715&&cut_costh_lep_pi);
        t.Draw("Umissfit>>kpi0enu",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_201&&cut_costh_lep_pi);
        t.Draw("Umissfit>>kpipi0pi0", (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_104&&cut_costh_lep_pi);
        t.Draw("Umissfit>>kpipi0",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_102&&cut_costh_lep_pi);
        t.Draw("Umissfit>>other",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_costh_lep_pi&&!(cut_455||cut_103||cut_105||cut_715||cut_201||cut_104||cut_1032||cut_102));  }
    if(nvar==4) {
        t.Draw("Umissfit>>sig",       (cut0||cut1||cut3)&&cuta&&cutc&&cut_sig&&cut_K1&&cut_costh_miss_neut);
        //t.Draw("Umissfit>>kpiw",      (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_455&&cut_costh_miss_neut);
        //t.Draw("Umissfit>>kpiw2",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_1032&&cut_costh_miss_neut);
        t.Draw("Umissfit>>kpipipi",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_103||cut_1032)&&cut_costh_miss_neut);
        t.Draw("Umissfit>>kpipipipi0",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_105||cut_455)&&cut_costh_miss_neut);
        t.Draw("Umissfit>>kpienu",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_715&&cut_costh_miss_neut);
        t.Draw("Umissfit>>kpi0enu",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_201&&cut_costh_miss_neut);
        t.Draw("Umissfit>>kpipi0pi0", (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_104&&cut_costh_miss_neut);
        t.Draw("Umissfit>>kpipi0",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_102&&cut_costh_miss_neut);
        t.Draw("Umissfit>>other",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_costh_miss_neut&&!(cut_455||cut_103||cut_105||cut_715||cut_201||cut_104||cut_1032||cut_102));  }
    if(nvar==5) {
        t.Draw("Umissfit>>sig",       (cut0||cut1||cut3)&&cuta&&cutc&&cut_sig&&cut_K1&&cut_chi2_eop);
        //t.Draw("Umissfit>>kpiw",      (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_455&&cut_chi2_eop);
        //t.Draw("Umissfit>>kpiw2",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_1032&&cut_chi2_eop);
        t.Draw("Umissfit>>kpipipi",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_103||cut_1032)&&cut_chi2_eop);
        t.Draw("Umissfit>>kpipipipi0",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_105||cut_455)&&cut_chi2_eop);
        t.Draw("Umissfit>>kpienu",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_715&&cut_chi2_eop);
        t.Draw("Umissfit>>kpi0enu",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_201&&cut_chi2_eop);
        t.Draw("Umissfit>>kpipi0pi0", (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_104&&cut_chi2_eop);
        t.Draw("Umissfit>>kpipi0",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_102&&cut_chi2_eop);
        t.Draw("Umissfit>>other",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_chi2_eop&&!(cut_455||cut_103||cut_105||cut_715||cut_201||cut_104||cut_1032||cut_102));  }
    /*if(nvar==6) {
        t.Draw("Umissfit>>sig",       (cut0||cut1||cut3)&&cuta&&cutc&&cut_sig&&cut_eop_E);
        t.Draw("Umissfit>>kpiw",      (cut0||cut1||cut3)&&cuta&&cutc&&!cut_sig&&cut_455&&cut_eop_E);
        t.Draw("Umissfit>>kpipipi",   (cut0||cut1||cut3)&&cuta&&cutc&&!cut_sig&&cut_103&&cut_eop_E);
        t.Draw("Umissfit>>kpipipipi0",(cut0||cut1||cut3)&&cuta&&cutc&&!cut_sig&&cut_105&&cut_eop_E);
        t.Draw("Umissfit>>kpienu",    (cut0||cut1||cut3)&&cuta&&cutc&&!cut_sig&&cut_715&&cut_eop_E);
        t.Draw("Umissfit>>kpi0enu",   (cut0||cut1||cut3)&&cuta&&cutc&&!cut_sig&&cut_201&&cut_eop_E);
        t.Draw("Umissfit>>kpipi0pi0", (cut0||cut1||cut3)&&cuta&&cutc&&!cut_sig&&cut_104&&cut_eop_E);
        t.Draw("Umissfit>>kpipi0",    (cut0||cut1||cut3)&&cuta&&cutc&&!cut_sig&&cut_102&&cut_eop_E);
        t.Draw("Umissfit>>other",     (cut0||cut1||cut3)&&cuta&&cutc&&!cut_sig&&cut_eop_E&&!(cut_455||cut_103||cut_105||cut_715||cut_201||cut_104||cu
t_102));  }
    */
    if(nvar==6) {
        t.Draw("Umissfit>>sig",       (cut0||cut1||cut3)&&cuta&&cutc&&cut_sig&&cut_K1&&cut_m_pipi);
        //t.Draw("Umissfit>>kpiw",      (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_455&&cut_m_pipi);
        //t.Draw("Umissfit>>kpiw2",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_1032&&cut_m_pipi);
        t.Draw("Umissfit>>kpipipi",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_103||cut_1032)&&cut_m_pipi);
        t.Draw("Umissfit>>kpipipipi0",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_105||cut_455)&&cut_m_pipi);
        t.Draw("Umissfit>>kpienu",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_715&&cut_m_pipi);
        t.Draw("Umissfit>>kpi0enu",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_201&&cut_m_pipi);
        t.Draw("Umissfit>>kpipi0pi0", (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_104&&cut_m_pipi);
        t.Draw("Umissfit>>kpipi0",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_102&&cut_m_pipi);
        t.Draw("Umissfit>>other",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_m_pipi&&!(cut_455||cut_103||cut_105||cut_715||cut_201||cut_104||cut_1032||cut_102));  }
//////////////////////////////////cut////////////////////////////////////

///////////////////////////////////after/////////////////////////////////
    if(nvar==7) {
        t.Draw("umiss2fit>>sig",       (cut0||cut1||cut3)&&cuta&&cutc&&cut_sig&&cut_K1&&cut_deltaE_pipzswp);
        //t.Draw("umiss2fit>>kpiw",      (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_455&&cut_deltaE_pipzswp);
        //t.Draw("umiss2fit>>kpiw2",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_1032&&cut_deltaE_pipzswp);
        t.Draw("umiss2fit>>kpipipi",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_103||cut_1032)&&cut_deltaE_pipzswp);
        t.Draw("umiss2fit>>kpipipipi0",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&(cut_105||cut_455)&&cut_deltaE_pipzswp);
        t.Draw("umiss2fit>>kpienu",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_715&&cut_deltaE_pipzswp);
        t.Draw("umiss2fit>>kpi0enu",   (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_201&&cut_deltaE_pipzswp);
        t.Draw("umiss2fit>>kpipi0pi0", (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_104&&cut_deltaE_pipzswp);
        t.Draw("umiss2fit>>kpipi0",    (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_102&&cut_deltaE_pipzswp);
        t.Draw("umiss2fit>>other",     (cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_deltaE_pipzswp&&!(cut_455||cut_103||cut_105||cut_715||cut_201||cut_104||cut_1032||cut_102));  }

///////////////////////////////////after/////////////////////////////////

     //cout<<kpiw->GetEntries()<<"\t"<<kpipipi->GetEntries()<<"\t"<<kpipipipi0->GetEntries()<<"\t"<<kpienu->GetEntries()<<"\t"<<kpi0enu->GetEntries()<<"\t"<<kpipi0pi0->GetEntries()<<"\t"<<kpipi0->GetEntries()<<"\t"<<other->GetEntries()<<endl;
    cout<<sig->GetEntries()<<endl;
    
    //kpiw2->Add(kpipipipi0);
    //kpipipi->Add(kpiw2); 
    kpipipipi0->Add(kpipipi);
    kpienu->Add(kpipipipi0);
    kpi0enu->Add(kpienu);  
    kpipi0pi0->Add(kpi0enu);
    kpipi0->Add(kpipi0pi0);
    other->Add(kpipi0); 
    sig->Add(other);
    
    sig->Draw();
    sig->SetFillColor(kRed-4);
    other->Draw("same");
    other->SetFillColor(kGreen-6);
    kpipi0->Draw("same");
    kpipi0->SetFillColor(kBlue-7);
    kpipi0pi0->Draw("same");
    kpipi0pi0->SetFillColor(kMagenta-7);
    kpi0enu->Draw("same");
    kpi0enu->SetFillColor(kBlue+3);
    kpienu->Draw("same");
    kpienu->SetFillColor(kYellow);
    kpipipipi0->Draw("same");
    kpipipipi0->SetFillColor(kCyan-3);
    kpipipi->Draw("same");
    kpipipi->SetFillColor(kBlue-9);
    //kpiw2->Draw("same");
    //kpiw2->SetFillColor(kBlue-9);     
    //kpiw->Draw("same");
    //kpiw->SetFillColor(kCyan-3);     

    leg = new TLegend(0.16,0.38,0.4,0.88);//(x0,y0,x1,y1)
    leg->SetTextFont(22);
    leg->SetTextSize(0.06);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetHeader("inclusive MC");
    leg->AddEntry(kpipipi,"K#pi#pi#pi","f");
    leg->AddEntry(kpipipipi0,"K#pi#pi#pi#pi^{0}","f");
    leg->AddEntry(kpienu,"K#pie#nu_{e}","f");
    leg->AddEntry(kpi0enu,"K#pi^{0}e#nu_{e}","f");
    leg->AddEntry(kpipi0pi0,"K#pi#pi^{0}#pi^{0}","f");
    leg->AddEntry(kpipi0,"K#pi#pi^{0}","f");
    leg->AddEntry(other,"other bkg","f");
    leg->AddEntry(sig,"signal","f");
    leg->Draw();

    PLabel1X = new TLegend(0.50,0.02,0.60,0.09);
    PLabel1X->SetTextAlign(22);
    PLabel1X->SetTextFont(22);
    PLabel1X->SetBorderSize(0);
    PLabel1X->SetTextSize(0.05);
    PLabel1X->SetTextColor(1);
    PLabel1X->SetLineColor(0);
    PLabel1X->SetFillColor(0);
    PLabel1X->SetLineWidth(0);
    PLabel1X->SetTextAngle(0);
    PLabel1X->SetHeader("U_{missfit} (GeV)");
    PLabel1X->Draw();

    if(nvar==0)  PLabel1Y = new TLegend(0.0,0.50,0.1,0.60);// after cut
    else         PLabel1Y = new TLegend(0.,0.50,0.06,0.60);  //before cut
    PLabel1Y->SetTextAlign(22);
    PLabel1Y->SetTextFont(22);
    PLabel1Y->SetBorderSize(0);
    PLabel1Y->SetTextSize(0.05);
    PLabel1Y->SetTextColor(1);
    PLabel1Y->SetLineColor(0);
    PLabel1Y->SetFillColor(0);
    PLabel1Y->SetLineWidth(0);
    PLabel1Y->SetTextAngle(90);
    PLabel1Y->SetHeader("Events / (0.005 GeV^{2}/#font[12]{c}^{4})");
    PLabel1Y->Draw();
   
    TLatex lt;
    lt.SetNDC();
    lt.SetTextAngle(0);
    lt.SetTextSize(0.06);
    lt.DrawLatex(0.72, 0.82,Form("N_{sig} = %.0f", sig->GetEntries()-other->GetEntries()));
    lt.DrawLatex(0.72, 0.76,Form("N_{bkg} = %.0f", other->GetEntries()));     
      TLatex ltm;
      ltm.SetNDC();
      ltm.SetTextAngle(0);
      ltm.SetTextSize(0.05);

    if(nvar==0)     ltm.DrawLatex(0.60,0.82,"");
    if(nvar==1)     ltm.DrawLatex(0.65,0.82,"(c)");
    if(nvar==2)     ltm.DrawLatex(0.65,0.82,"(c)");
    if(nvar==3)     ltm.DrawLatex(0.65,0.82,"(c)");
    if(nvar==4)     ltm.DrawLatex(0.65,0.82,"(c)");
    if(nvar==5)     ltm.DrawLatex(0.65,0.82,"(d)");
    if(nvar==6)     ltm.DrawLatex(0.65,0.82,"(b)");
    //if(nvar==7)     ltm.DrawLatex(0.65,0.82,"(c)");
 
    c->Modified();
    c->cd();
    c->SetSelected(c);
    //c->Print(Form(".umiss2fit_stack_bkg_%d.eps",nvar));
    c->Print("umiss2fit_stack_bkg.eps");
}
}
