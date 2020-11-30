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

void misidentify103()
{

  //truth~ tagD
  TCut cut_sig = "isTrueTag&&isSignal&&abs(KpTrueID)==321&&abs(pi1TrueID)==211&&abs(pi2TrueID)==211&&abs(elecTrueID)==11&&abs(elecMomID)==421";
  TCut cut_K1  = "((mcmodea==666||mcmodea==222)&&mcmode1==103)||((mcmodeb==-666||mcmodeb==-222)&&mcmode2==-103)";
  TCut cuta    = "elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge";
  TCut cutb    = "isbest==1&&type==1&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)&&m_pipi>0.31&&abs(m_pipi-0.497611)>0.01&&abs(deltaE_pipzswp)>0.012";
  TCut cutc    = "abs(Umissfit)<0.1&&!(mcmode1==104||mcmode2==-104)";
  TCut cut0    = "mode==0&&clveto==1&&deltaE>-0.029&&deltaE<0.027";
  TCut cut1    = "mode==1&&deltaE>-0.069&&deltaE<0.038";
  TCut cut3    = "mode==3&&deltaE>-0.031&&deltaE<0.028";
  TCut cut_k3pi= "(charm==1&&mcmodeb==-103)||(charm==-1&&mcmodea==103)"; 
  //generate~ truth
  TCut cut_mc0  = "mcmodeb==-101 && (mcmodea==103||mcmodea==108)" ;
  TCut cut_mc01 = "mcmodea==101 && (mcmodeb==-103||mcmodeb==-108)" ;
  TCut cut_mc1  = "mcmodeb==-102 && (mcmodea==103||mcmodea==108)" ;
  TCut cut_mc11 = "mcmodea==102 && (mcmodeb==-103||mcmodeb==-108)" ;
  TCut cut_mc3  = "(mcmodeb==-103||mcmodeb==-108)&&(mcmodea==103||mcmodea==108)" ;
  TCut cut_mc31 = "(mcmodea==103||mcmodea==108) && (mcmodeb==-103||mcmodeb==-108)" ;
  //
  //cut
  TChain t0("tagD");
  t0.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/K3pi/red_truth/truth*.root");
  TTree *T0_tagD = t0.CopyTree(!(cut_sig&&cut_K1)&&cuta&&cutb&&cutc&&cut0&&cut_k3pi);//mode0
  TTree *T1_tagD = t0.CopyTree(!(cut_sig&&cut_K1)&&cuta&&cutb&&cutc&&cut1&&cut_k3pi);
  TTree *T3_tagD = t0.CopyTree(!(cut_sig&&cut_K1)&&cuta&&cutb&&cutc&&cut3&&cut_k3pi);
  
  //TTree *T0_tagD = t0.CopyTree(cuta&&cutb&&cutc&&cut0&&cut_k3pi);//mode0
  //TTree *T1_tagD = t0.CopyTree(cuta&&cutb&&cutc&&cut1&&cut_k3pi);
  //TTree *T3_tagD = t0.CopyTree(cuta&&cutb&&cutc&&cut3&&cut_k3pi);
  TChain t1("truth");
  TChain t2("truth");
  t1.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/K3pi/DDecay_root/d0/*.root");
  t1.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/K3pi/KsKpi_s/root/D0_*.root");
  t2.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/K3pi/DDecay_root/d0bar/*.root");
  t2.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/K3pi/KsKpi_s/root/D0bar_*.root");
  t2.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/K3pi/KsKpi_d/KsKpi.root");//can not distinguish
  TTree *T0_truth  = t1.CopyTree(cut_mc0);//mode0
  TTree *T01_truth = t2.CopyTree(cut_mc01);
  TTree *T1_truth = t1.CopyTree(cut_mc1);
  TTree *T11_truth = t2.CopyTree(cut_mc11);
  TTree *T3_truth = t1.CopyTree(cut_mc3);
  TTree *T31_truth = t2.CopyTree(cut_mc31);
  //DT eff
  double ST0_eff = 0.66089, ST1_eff = 0.349938, ST3_eff = 0.410512, ST0_eff_err = 0.00277811, ST1_eff_err = 0.00158968, ST3_eff_err = 0.000517714;
  double N0obs     = T0_tagD->GetEntries(),  N1obs = T1_tagD->GetEntries(),    N3obs = T3_tagD->GetEntries(),
         N0obs_err = sqrt(N0obs),            N1obs_err = sqrt(N1obs),          N3obs_err = sqrt(N3obs),
         N0total   = T0_truth->GetEntries(), N1total = T1_truth->GetEntries(), N3total = T3_truth->GetEntries(),
         N01total  = T01_truth->GetEntries(),N11total= T11_truth->GetEntries(),N31total= T31_truth->GetEntries();
  /******DT mis_eff*********/
  double misidentify0_DTeff   = N0obs / (N0total+N01total), 
         misidentify1_DTeff = N1obs / (N1total+N11total), 
         misidentify3_DTeff = N3obs /(N3total+N31total);
  /*******ST mis_eff*******/
  double misidentify0_eff = misidentify0_DTeff/ST0_eff, 
         misidentify1_eff = misidentify1_DTeff/ST1_eff, 
         misidentify3_eff = misidentify3_DTeff/ST3_eff;
  /******DT mis_eff err*******/
  double misidentify0_DTeff_err = N0obs_err/(N0total+N01total),   
         misidentify1_DTeff_err = N1obs_err/(N1total+N11total),  
         misidentify3_DTeff_err=N3obs_err/(N3total+N31total);
  /******ST mis_eff err******/
  double misidentify0_eff_err = Err(misidentify0_DTeff, misidentify0_DTeff_err, ST0_eff, ST0_eff_err), 
         misidentify1_eff_err = Err(misidentify1_DTeff, misidentify1_DTeff_err, ST1_eff, ST1_eff_err), 
         misidentify3_eff_err = Err(misidentify3_DTeff, misidentify3_DTeff_err, ST3_eff, ST3_eff_err);
  /****Nst * bf * mis_eff * correction factor*****/
  double bf_k3pi = 0.0823+0.0033*0.692,    bf_k3pi_err = sqrt(0.0014*0.0014+0.0005*0.0005) , correction = 0.88,
         N_st0 = 6.00716e+06/10.8, N_st1 = 1.13678e+07/10.8, N_st3 =7.71511e+06/10.8;
  double N0 = (N_st0)*bf_k3pi*misidentify0_eff*0.88, 
         N1 = (N_st1)*bf_k3pi*misidentify1_eff*0.88, 
         N3 = (N_st3)*bf_k3pi*misidentify3_eff*0.88;
  double N0_err = N_st0*0.88*err(bf_k3pi, bf_k3pi_err, misidentify0_eff, misidentify0_eff_err), 
         N1_err = N_st1*0.88*err(bf_k3pi, bf_k3pi_err, misidentify1_eff, misidentify1_eff_err), 
         N3_err = N_st3*0.88*err(bf_k3pi, bf_k3pi_err, misidentify3_eff, misidentify3_eff_err);

  //cout
  cout<<"produced number in sigMC: "<<endl;
  cout<<N0total+N01total<<"\t"<<N1total+N11total<<"\t"<<N3total+N31total<<endl;
  cout<<"misidntified number in sigMC: "<<endl;
  cout<<N0obs<<"\t"<<N1obs<<"\t"<<N3obs<<endl;
  cout<<"misidentify eff:  "<<endl;
  cout<<misidentify0_eff<<"+/-"<<misidentify0_eff_err<<endl;
  cout<<misidentify1_eff<<"+/-"<<misidentify1_eff_err<<endl;
  cout<<misidentify3_eff<<"+/-"<<misidentify3_eff_err<<endl;
  cout<<"misidentify number:"<<endl;
  cout<<N0<<"+/-"<<N0_err<<endl;
  cout<<N1<<"+/-"<<N1_err<<endl;
  cout<<N3<<"+/-"<<N3_err<<endl;

  //write
  ofstream out;
  out.open("./dat/misidentify103.txt");//mode0
  out<<"produced number in sigMC: "<< N0total+N01total<<"\t"<<N1total+N11total<<"\t"<<N3total+N31total<<"\n";
  out<<"misidentified number in sigMC: "<< N0obs<<"\t"<<N1obs<<"\t"<<N3obs<<"\n";
  out<<"total mididentified number: "<< N0obs+N1obs+N3obs<<"\n";
  out<<"misidentified efficiency : "<<"\n";
  out<<misidentify0_eff<<"\t"<<misidentify0_eff_err<<"\n";
  out<<misidentify1_eff<<"\t"<<misidentify1_eff_err<<"\n";
  out<<misidentify3_eff<<"\t"<<misidentify3_eff_err<<"\n"; 
  out<<"misidentify number:"<<"\n";
  out<<N0<<"+/-"<<N0_err<<"\n";
  out<<N1<<"+/-"<<N1_err<<"\n";
  out<<N3<<"+/-"<<N3_err;

  ofstream out1;
  out1.open("./dat/misidentify103.dat");
  out1<<N0<<"\n";
  out1<<N1<<"\n";
  out1<<N3<<"\n";


}
// x1*y1
double err(double x1, double x_err, double y1, double y_err)
{
    double c = sqrt(x1*x1*y_err*y_err+y1*y1*x_err*x_err);
    return c;
}

//p1/q1
double Err(double p1, double p_err, double q1, double q_err)
{
    double d = sqrt(p1*p1*q_err*q_err+q1*q1*p_err*p_err)/(q1*q1);
    return d;
}
