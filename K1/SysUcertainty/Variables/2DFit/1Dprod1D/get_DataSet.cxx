void get_DataSet()
{
    using namespace RooFit;
 
    TCut normal  = "isbest==1&&type==1&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)&&m_pipi>0.31&&abs(m_pipi-0.497611)>0.01&&abs(deltaE_pipzswp)>0.012&&abs(Umissfit)<0.1&&elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge";
    TCut m0 = "mode==0&&clveto==1&&deltaE>-0.029&&deltaE<0.027";
    TCut m1 = "mode==1&&deltaE>-0.069&&deltaE<0.038";
    TCut m3 = "mode==3&&deltaE>-0.031&&deltaE<0.028";
    TCut cuta = "mhad>0.9&&mhad<1.68";

    //RooRealVar Umissfit("Umissfit","", -0.1, 0.1 );
    RooRealVar mhad("mhad","", 0.9,1.68);
    TChain *ch0 = new TChain("tagD");
    ch0->Add("/scratchfs/bes/fangyl/data1/data/all/*.root");
    TTree* t0 = ch0->CopyTree(normal&&cuta&&(m0||m1||m3));
    
    RooDataSet *ds_data = new RooDataSet("ds_data","",t0,mhad);
    cout<<t0->GetEntries()<<endl; 
    TFile *f0 = new TFile("./figure/DataSet_mhad.root","recreate");
    //TFile *f0 = new TFile("./figure/DataSet_mmiss.root","recreate");
    ds_data->Write();
    f0->Close();
}
