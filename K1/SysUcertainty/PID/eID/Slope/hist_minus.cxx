#include <fstream>
#include <iostream>
using namespace std;

void hist_minus()
{
   TChain t("Bhabha");
   t.Add("/scratchfs/bes/fangyl/data1/uncertainty/all/*.root");//data
   //t.Add("/scratchfs/bes/fangyl/inclusive_mc/uncertainty/A-bhabha*.root");//inclusive mc
 
   Double_t the_rec, p_rec,Eemc_rec, eop_E, dedxchi_E;
   Double_t found;
   t.SetBranchAddress("the_rec", &the_rec);
   t.SetBranchAddress("p_rec", &p_rec);
   t.SetBranchAddress("Eemc_rec", &Eemc_rec);
   t.SetBranchAddress("eop_E", &eop_E);
   t.SetBranchAddress("dedxchi_E", &dedxchi_E);
   t.SetBranchAddress("found", &found);
   

   Double_t cos, p;
   Double_t isfound;
   TTree* Bhabha = new TTree("Bhabha","");
   Bhabha->Branch("cos",&cos);
   Bhabha->Branch("p",&p);
   Bhabha->Branch("isfound",&isfound);
   int g_index=0;
   for( unsigned int i=0; i<t.GetEntries(); i++ ){
      t.GetEntry(i);
      bool cut_eop = (eop_E-0.32)>0.16*dedxchi_E;
      bool cut_p   = p_rec>0&&p_rec<2.0;
      bool cut_cos = cos(the_rec)>-1.0&&cos(the_rec)<1.0;
      bool cut = cut_eop&&cut_p&&cut_cos;
  if(cut){
     p = p_rec;
     cos = cos(the_rec);
     isfound = found;
     Bhabha->Fill();
    }    
  }
   TFile *f = new TFile("./dat/data_bhabha.root","recreate");
   //TFile *f = new TFile("./dat/data_mc_bhabha.root","recreate"); 
   Bhabha->Write();
   f->Close();
}
