#include <fstream>
#include <iostream>
using namespace std;
void sel_truth_new()
{

    	TChain mc("truth");
      
        mc.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/root/*.root");
        Int_t mcmode1, mcmode2, mcmodea, mcmodeb;
        Double_t ptkK, ptkPi1, ptkPi, ptkE, cosE;
      
        mc.SetBranchAddress("mcmodea",&mcmodea);
        mc.SetBranchAddress("mcmodeb",&mcmodeb);

        mc.SetBranchAddress("mcmode1",&mcmode1);
        mc.SetBranchAddress("mcmode2",&mcmode2);
       
       Int_t nentries = (Int_t)mc.GetEntries();
        //TFile *newfile = new TFile("truth.root","recreate");
        TFile *newfile = new TFile("truth0.root","recreate");
    	TTree *newtree = mc.CloneTree(0);


    	for (Int_t i=0; i<nentries; i++) {
             mc.GetEntry(i);
             if( ( (mcmodea==666||mcmodea==222)&&mcmode1==103&&(mcmodeb==-101||mcmodeb==-102||mcmodeb==-103||mcmodeb==-1032)  ) ||  ( (mcmodeb==-666||mcmodeb==-222)&&mcmode2==-103&&(mcmodea==101||mcmodea==102||mcmodea==103||mcmodea==1032)  )  ){
             newtree->Fill();
             }

         }
        cout << newtree->GetEntries() << endl;
        newtree->Write();
        newfile->Close();
}
