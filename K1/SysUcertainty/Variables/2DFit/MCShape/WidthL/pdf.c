#include <fstream>
#include <iostream>
using namespace std;
void pdf(){
for (int nmode=0; nmode<3;nmode++) {
   TChain t("tagD");
   t.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/mc_shape/width_low/red_truth/truth.root");

   int mode,clveto,isbest,type ;
   int mcmode1, mcmode2, mcmodea, mcmodeb;
   double mBC, deltaE, eop_E, mhad ,Umissfit, dedxchi_E, m_pipi, costh_lep_pi, costh_miss_neut, mD_EasPiwPi0, mD_EasPi, deltaE_pipzswp;
   int isTrueTag, isSignal, KpTrueID, pi1TrueID, pi2TrueID, elecTrueID, elecMomID, elecCharge, kaonCharge, pionCharge, pion2Charge, charm;
  
   t.SetBranchAddress("isTrueTag",&isTrueTag);
   t.SetBranchAddress("isSignal",&isSignal);
   t.SetBranchAddress("KpTrueID",&KpTrueID);
   t.SetBranchAddress("pi1TrueID",&pi1TrueID);
   t.SetBranchAddress("pi2TrueID",&pi2TrueID);
   t.SetBranchAddress("elecTrueID",&elecTrueID);
   t.SetBranchAddress("elecMomID",&elecMomID);

   t.SetBranchAddress("mcmode1",&mcmode1);
   t.SetBranchAddress("mcmode2",&mcmode2);
   t.SetBranchAddress("mcmodea",&mcmodea);
   t.SetBranchAddress("mcmodeb",&mcmodeb);

   t.SetBranchAddress("elecCharge",&elecCharge);
   t.SetBranchAddress("kaonCharge",&kaonCharge);
   t.SetBranchAddress("pionCharge",&pionCharge);
   t.SetBranchAddress("pion2Charge",&pion2Charge);
   t.SetBranchAddress("charm",&charm);
 
   t.SetBranchAddress("mode",&mode);
   t.SetBranchAddress("clveto",&clveto);
   t.SetBranchAddress("isbest",&isbest);
   t.SetBranchAddress("type",&type);

   t.SetBranchAddress("mBC",&mBC);
   t.SetBranchAddress("deltaE",&deltaE);
   t.SetBranchAddress("eop_E",&eop_E);
   t.SetBranchAddress("mhad",&mhad);
   t.SetBranchAddress("Umissfit",&Umissfit);
   t.SetBranchAddress("dedxchi_E",&dedxchi_E);
   t.SetBranchAddress("m_pipi",&m_pipi);
   t.SetBranchAddress("costh_lep_pi",&costh_lep_pi);
   t.SetBranchAddress("costh_miss_neut",&costh_miss_neut);
   t.SetBranchAddress("mD_EasPiwPi0",&mD_EasPiwPi0);
   t.SetBranchAddress("mD_EasPi",&mD_EasPi);
   t.SetBranchAddress("deltaE_pipzswp",&deltaE_pipzswp);

   Double_t Umissfit, mhad;
   TTree* sig = new TTree("sig","");   
 
   sig->Branch("Umissfit",&Umissfit); 
   sig->Branch("mhad",&mhad);
    
   int g_index=0;
   for( unsigned int i=0; i<t.GetEntries(); i++ ){
       t.GetEntry(i);
       bool cut_sig = isTrueTag&&isSignal&&abs(KpTrueID)==321&&abs(pi1TrueID)==211&&abs(pi2TrueID)==211&&abs(elecTrueID)==11&&abs(elecMomID)==421&&(((mcmodea==666||mcmodea==222)&&mcmode1==103)||((mcmodeb==-666||mcmodeb==-222)&&mcmode2==-103)); // truth massage in mc
       bool K_2  = mcmode1==104 || mcmode2==-104;
       bool cutc = elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge;
       bool cut0 = mode==0&&clveto==1&&deltaE>-0.029&&deltaE<0.027;
       bool cut1 = mode==1&&deltaE>-0.069&&deltaE<0.038;
       bool cut3 = mode==3&&deltaE>-0.031&&deltaE<0.028;
       bool cuta = isbest==1&&type==1&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)&&m_pipi>0.31&&fabs(m_pipi-0.497611)>0.01&&fabs(deltaE_pipzswp)>0.012;
       bool cutb =fabs(Umissfit)<0.1;
       bool cut_k3pipi0 = (charm==1&&(mcmodeb==-105||mcmodeb==-455))||(charm==-1&&(mcmodea==105||mcmodea==455));
       bool cut_k3pi    = (charm==1&&mcmodeb==-103)||(charm==-1&&mcmodea==103);
       
       //without mhad cut

       if(nmode==0)   bool sig_shape  = cuta&&cutb&&cutc&&cut_sig&&cut0;//mode0
       if(nmode==1)   bool sig_shape  = cuta&&cutb&&cutc&&cut_sig&&cut1;
       if(nmode==2)   bool sig_shape  = cuta&&cutb&&cutc&&cut_sig&&cut3;


    if(sig_shape){ 
      Umissfit = Umissfit;
      mhad = mhad;  
      sig->Fill();  
     }
      
      if (g_index%10000==0) std::cout << g_index << "\t" << i << std::endl;
      g_index++;
   }  
      TFile *f = new TFile(Form("./pdf/sig_%d.root",nmode),"recreate");

      sig->Write();
      f->Close(); 
      cout<<sig->GetEntries()<<endl;
}
}




