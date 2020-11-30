void preSelect()
{
    TChain ch("tagD");
    ch.Add("/scratchfs/bes/fangyl/data1/data/all/*.root");
    TCut normal  = "isbest==1&&type==1&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)&&m_pipi>0.31&&abs(m_pipi-0.497611)>0.01&&abs(deltaE_pipzswp)>0.012&&abs(Umissfit)<0.1&&elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge";   
    TCut m0 = "mode==0&&clveto==1&&deltaE>-0.029&&deltaE<0.027";
    TCut m1 = "mode==1&&deltaE>-0.069&&deltaE<0.038";
    TCut m3 = "mode==3&&deltaE>-0.031&&deltaE<0.028";
    TCut cuta = "mhad<1.68&&mhad>0.9";
    TTree* t = ch.CopyTree(normal&&cuta&&(m0||m1||m3));
    TFile *f = new TFile("figure/data.root","recreate");
    t->Write();
    f->Close();
    delete f;
    delete t;
}
