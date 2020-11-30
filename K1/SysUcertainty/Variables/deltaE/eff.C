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

using namespace RooFit ;
using std::vector;

void eff()
{

  //truth~ tagD
  TCut cut_sig = "isTrueTag&&isSignal&&abs(KpTrueID)==321&&abs(pi1TrueID)==211&&abs(pi2TrueID)==211&&abs(elecTrueID)==11&&abs(elecMomID)==421";
  TCut cuta    = "elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge";
  TCut cutb_minus    = "isbest==1&&type==1&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)&&m_pipi>0.31&&abs(m_pipi-0.497611)>0.01&&abs(deltaE_pipzswp)>0.008";
  TCut cutb_plus     = "isbest==1&&type==1&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)&&m_pipi>0.31&&abs(m_pipi-0.497611)>0.01&&abs(deltaE_pipzswp)>0.016";
  TCut cutc    = "abs(Umissfit)<0.1";
  TCut cut0    = "mode==0&&clveto==1&&deltaE>-0.029&&deltaE<0.027";
  TCut cut1    = "mode==1&&deltaE>-0.069&&deltaE<0.038";
  TCut cut3    = "mode==3&&deltaE>-0.031&&deltaE<0.028";
  TCut cut_K1  = "((mcmodea==666||mcmodea==222)&&mcmode1==103)||((mcmodeb==-666||mcmodeb==-222)&&mcmode2==-103)";
  //generate~ truth
  TCut cut_mc0  = "(mcmodea==666||mcmodea==222) && mcmodeb==-101 && mcmode1==103";
  TCut cut_mc01 = "mcmodea==101 && (mcmodeb==-666||mcmodeb==-222) && mcmode2==-103";
  TCut cut_mc1  = "(mcmodea==666||mcmodea==222) && mcmodeb==-102 && mcmode1==103";
  TCut cut_mc11 = "mcmodea==102 && (mcmodeb==-666||mcmodeb==-222) && mcmode2==-103";
  TCut cut_mc3  = "(mcmodea==666||mcmodea==222) && (mcmodeb==-103 || mcmodeb==-108) && mcmode1==103";
  TCut cut_mc31 = "(mcmodea==103 || mcmodea==108) && (mcmodeb==-666||mcmodeb==-222) && mcmode2==-103";
  //cut
  TChain t0("tagD");
  t0.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/BELLE/fit2/red_truth/truth.root");
  for(int m=0; m <2; m++){
  if (m == 0) {
  TTree *T0_tagD = t0.CopyTree(cut_sig&&cuta&&cutb_minus&&cutc&&cut_K1&&cut0);//mode0
  TTree *T1_tagD = t0.CopyTree(cut_sig&&cuta&&cutb_minus&&cutc&&cut_K1&&cut1);
  TTree *T3_tagD = t0.CopyTree(cut_sig&&cuta&&cutb_minus&&cutc&&cut_K1&&cut3);
  }
  if (m == 1) {
  TTree *T0_tagD = t0.CopyTree(cut_sig&&cuta&&cutb_plus&&cutc&&cut_K1&&cut0);//mode0
  TTree *T1_tagD = t0.CopyTree(cut_sig&&cuta&&cutb_plus&&cutc&&cut_K1&&cut1);
  TTree *T3_tagD = t0.CopyTree(cut_sig&&cuta&&cutb_plus&&cutc&&cut_K1&&cut3);
  }
  TChain t1("truth");
  t1.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/BELLE/fit2/DDecay_root/d0/*.root");
  t1.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/BELLE/fit2/KsKpi/d0/*.root");
  TChain t2("truth");
  t2.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/BELLE/fit2/DDecay_root/d0bar/*.root");
  t2.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/BELLE/fit2/KsKpi/d0bar/*.root");
  TTree *T0_truth  = t1.CopyTree(cut_mc0);//mode0
  TTree *T01_truth = t2.CopyTree(cut_mc01);
  TTree *T1_truth = t1.CopyTree(cut_mc1);
  TTree *T11_truth = t2.CopyTree(cut_mc11);
  TTree *T3_truth = t1.CopyTree(cut_mc3);
  TTree *T31_truth = t2.CopyTree(cut_mc31);
  //DT eff
  double N0obs     = T0_tagD->GetEntries(), N1obs = T1_tagD->GetEntries(),    N3obs = T3_tagD->GetEntries();
  double N0obs_err = sqrt(N0obs),           N1obs_err = sqrt(N1obs),          N3obs_err = sqrt(N3obs);
  double N0total   = T0_truth->GetEntries(), N1total = T1_truth->GetEntries(), N3total = T3_truth->GetEntries();
  double N01total  = T01_truth->GetEntries(),N11total= T11_truth->GetEntries(),N31total= T31_truth->GetEntries();
  double DT0_eff   = N0obs / (N0total+N01total), DT1_eff = N1obs / (N1total+N11total),   DT3_eff = N3obs /(N3total+N31total);
  double DT0_eff_err = N0obs_err/(N0total+N01total),   DT1_eff_err = N1obs_err/(N1total+N11total),  DT3_eff_err=N3obs_err/(N3total+N31total);
  //sig eff = ST eff / DT eff.    propagation of err : sig_eff_err = ( ST_eff*DT_eff_err + DT_eff*ST_eff_err ) / (ST_eff*ST_eff)
  /*double ST0_eff = 0.6537, ST1_eff = 0.3467, ST3_eff = 0.3820;// quoted from zhsf
  double ST0_eff_err = 0.0009, ST1_eff_err = 0.0004, ST3_eff_err = 0.0006;
*/double ST0_eff = 0.66089, ST1_eff = 0.349938, ST3_eff = 0.4105;
  double ST0_eff_err = 0.00277811, ST1_eff_err = 0.00158968, ST3_eff_err = 0.0005;
  double SL0_eff = DT0_eff / ST0_eff,  SL1_eff = DT1_eff / ST1_eff,  SL3_eff = DT3_eff / ST3_eff;
  double SL0_eff_err=sqrt(ST0_eff*ST0_eff*DT0_eff_err*DT0_eff_err + DT0_eff*DT0_eff*ST0_eff_err*ST0_eff_err)/(ST0_eff*ST0_eff);
  double SL1_eff_err=sqrt(ST1_eff*ST1_eff*DT1_eff_err*DT1_eff_err + DT1_eff*DT1_eff*ST1_eff_err*ST1_eff_err)/(ST1_eff*ST1_eff);
  double SL3_eff_err=sqrt(ST3_eff*ST3_eff*DT3_eff_err*DT3_eff_err + DT3_eff*DT3_eff*ST3_eff_err*ST3_eff_err)/(ST3_eff*ST3_eff);

  double N_st0 = 542153., N_st1 = 1080690., N_st3 =737036.,  N_ST_total = 2359879;//fanyl
  double eff_ave = (N_st0*SL0_eff + N_st1*SL1_eff + N_st3*SL3_eff) / N_ST_total;
  double eff_ave_err = (N_st0*SL0_eff_err + N_st1*SL1_eff_err + N_st3*SL3_eff_err) / N_ST_total;
  //cout
  cout<<__LINE__<<"**************outup*************"<<endl;
  cout<<"produced number in sigMC: "<<endl;
  cout<<N0total+N01total<<"\t"<<N1total+N11total<<"\t"<<N3total+N31total<<endl;
  cout<<"observed number in sigMC: "<<endl;
  cout<<N0obs<<"\t"<<N1obs<<"\t"<<N3obs<<endl;
  cout<<"DT_eff = "<< "N_obs / N_total"<<endl;
  cout<<"DT0_eff ="<<DT0_eff<<"\t"<<DT0_eff_err<<endl; 
  cout<<"DT1_eff ="<<DT1_eff<<"\t"<<DT1_eff_err<<endl;
  cout<<"DT3_eff ="<<DT3_eff<<"\t"<<DT3_eff_err<<endl;
  cout<<"SL_eff = "<<"DT eff / ST eff"<<endl;
  cout<<"SL0_eff = "<<SL0_eff<<"\t"<<SL0_eff_err<<endl;
  cout<<"SL1_eff = "<<SL1_eff<<"\t"<<SL1_eff_err<<endl;
  cout<<"SL3_eff = "<<SL3_eff<<"\t"<<SL3_eff_err<<endl;
  cout<<"average SL eff = "<<eff_ave<<"\t"<<eff_ave_err<<endl;
   
  //write
  ofstream read;
  if(m == 0)read.open("./dat/eff_minus.txt");
  if(m == 1)read.open("./dat/eff_plus.txt");
  read<<SL0_eff<<"\n";
  read<<SL1_eff<<"\n";
  read<<SL3_eff<<"\n";
  read<<eff_ave;
  

  }
  

}
