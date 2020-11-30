void sel()
{
    using namespace RooFit;
        TCut mode3 = "mode==3 && isbest==1 && mBC>1.8365 && mBC<1.8865 && deltaE>-0.031 && deltaE<0.028";   
        TCut event = "eventNo%10==5";
        TCut bkg   = "isTrueTag!=1";
        RooRealVar mBC("mBC","", 1.8365,1.8865);                                                          
        TChain t("tagD"); 
        t.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/new/inclusiveMC/umiss2fit/reduce_truth/only_truth/root/truth*.root");
        double bins = 200;
        TH1F* h0 = new TH1F("h0","",bins,1.8365,1.8865);
        t.Draw("mBC>>h0",mode3&&event&&bkg);
        RooDataHist* hbkg = new RooDataHist("hbkg","",mBC,h0);
        TFile *f0 = new TFile("bkg3.root","recreate");
        h0->Write();
        f0->Close();
}
