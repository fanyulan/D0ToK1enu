#include <fstream>
#include <iostream>
using namespace std;
void sel_tagD()
{
     TChain mc("tagD");
     mc.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/BELLE/fit2/only_truth/truth*.root");
     Int_t mode, clveto, type, isbest, mcmode1, mcmode2, mcmodea, mcmodeb;
     Int_t isTrueTag, isSignal, KpTrueID, pi1TrueID, pi2TrueID, elecTrueID, elecMomID, charm, pion2Charge, pionCharge, kaonCharge, elecCharge;
     Double_t  deltaE, mBC, Umissfit, mD_EasPi, mD_EasPiwPi0, costh_lep_pi, dedxchi_E, eop_E, costh_miss_neut, m_pipi, mhad, deltaE_pipzswp;
     Double_t  p_kaonp3, p_pion2p3, p_pionp3, elecP, cosE;
        mc.SetBranchAddress("mode",  &mode);
        mc.SetBranchAddress("clveto",&clveto);
        mc.SetBranchAddress("type",  &type);
        mc.SetBranchAddress("isbest",&isbest);
        mc.SetBranchAddress("mcmode1",&mcmode1);
        mc.SetBranchAddress("mcmode2",&mcmode2);
        mc.SetBranchAddress("mcmodea",&mcmodea);  
        mc.SetBranchAddress("mcmodeb",&mcmodeb);  
        mc.SetBranchAddress("isTrueTag",  &isTrueTag);
        mc.SetBranchAddress("isSignal",   &isSignal);
        mc.SetBranchAddress("KpTrueID",   &KpTrueID);
        mc.SetBranchAddress("pi1TrueID",  &pi1TrueID);
        mc.SetBranchAddress("pi2TrueID",  &pi2TrueID);
        mc.SetBranchAddress("elecTrueID", &elecTrueID);
        mc.SetBranchAddress("elecMomID",  &elecMomID);
        mc.SetBranchAddress("charm",      &charm);
        mc.SetBranchAddress("pion2Charge",&pion2Charge);
        mc.SetBranchAddress("pionCharge", &pionCharge);
        mc.SetBranchAddress("kaonCharge", &kaonCharge);
        mc.SetBranchAddress("elecCharge", &elecCharge);
        mc.SetBranchAddress("deltaE",         &deltaE);
        mc.SetBranchAddress("mBC",            &mBC);
        mc.SetBranchAddress("Umissfit",       &Umissfit);
        mc.SetBranchAddress("mD_EasPi",       &mD_EasPi);
        mc.SetBranchAddress("mD_EasPiwPi0",   &mD_EasPiwPi0);
        mc.SetBranchAddress("costh_lep_pi",   &costh_lep_pi);
        mc.SetBranchAddress("dedxchi_E",      &dedxchi_E);
        mc.SetBranchAddress("eop_E",          &eop_E);
        mc.SetBranchAddress("costh_miss_neut",&costh_miss_neut);
        mc.SetBranchAddress("m_pipi",         &m_pipi);
        mc.SetBranchAddress("mhad",           &mhad);
        mc.SetBranchAddress("deltaE_pipzswp", &deltaE_pipzswp);
        mc.SetBranchAddress("p_kaonp3", &p_kaonp3);
        mc.SetBranchAddress("p_pion2p3",&p_pion2p3);
        mc.SetBranchAddress("p_pionp3", &p_pionp3);
        mc.SetBranchAddress("elecP",    &elecP);
        mc.SetBranchAddress("cosE",     &cosE);

        Int_t nentries = (Int_t)mc.GetEntries();
	TFile *newfile = new TFile("tagD_pid.root","recreate");
	TTree *newtree = mc.CloneTree(0);
	for (Int_t i=0; i<nentries; i++) {
             mc.GetEntry(i);
     /////////////////////////////////////////////////////////   
     /*
     bool cut0 = "mode==0&&clveto&&deltaE>-0.029&&deltaE<0.027";
     bool cut1 = "mode==1&&deltaE>-0.069&&deltaE<0.038";
     bool cut3 = "mode==3&&deltaE>-0.031&&deltaE<0.028";
     bool cuta = "isbest==1&&type==1&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.8&&mD_EasPiwPi0<1.8&&costh_lep_pi<0.9&&costh_miss_neut<0.85&&((eop_E-0.6)>1./10*dedxchi_E)&& m_pipi>0.3&&fabs(m_pipi-0.497611)>0.005&&fabs(deltaE_pipzswp)>0.02";
     bool cutb = "fabs(Umissfit)<0.1&&p_kaonp3>0&&p_pionp3>0&&p_pion2p3>0";
     bool cutc = "isTrueTag&&isSignal&&abs(KpTrueID)==321&&abs(pi1TrueID)==211&&abs(pi2TrueID)==211&&abs(elecTrueID)==11&&abs(elecMomID)==421&&elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge&&((mcmodea==666||mcmodea==222)&&mcmode1==103)||((mcmodeb==-666||mcmodeb==-222)&&mcmode2==-103)";
 
        if((cut0||cut1||cut3)&&cuta&&cutb&&cutc){
         newtree->Fill();
      }
      *//////////////////////////////////////////////////////
      if((mode==0&&clveto&&deltaE>-0.029&&deltaE<0.027) || (mode==1&&deltaE>-0.069&&deltaE<0.038) || (mode==3&&deltaE>-0.031&&deltaE<0.028)) {
          if(isbest==1&&type==1&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.42)>0.14*dedxchi_E)&& m_pipi>0.31&&fabs(m_pipi-0.497611)>0.01&&fabs(deltaE_pipzswp)>0.014) {
              if(fabs(Umissfit)<0.1&&p_kaonp3>0&&p_pionp3>0&&p_pion2p3>0){
                  if(isTrueTag&&isSignal&&abs(KpTrueID)==321&&abs(pi1TrueID)==211&&abs(pi2TrueID)==211&&abs(elecTrueID)==11&&abs(elecMomID)==421&&elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge&&((mcmodea==666||mcmodea==222)&&mcmode1==103)||((mcmodeb==-666||mcmodeb==-222)&&mcmode2==-103)){ newtree->Fill(); }
                  } 
              }   
         } 
      }
	cout << newtree->GetEntries() << endl;;
	newtree->Write();
	newfile->Close();
}
