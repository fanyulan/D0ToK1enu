void smear()
{
    using namespace std;
    using namespace RooFit;

    RooRealVar Umissfit("Umissfit","", -0.1,0.1);
    //RooRealVar mhad("mhad","", 0.9,1.68);
    TFile *f0 = TFile::Open("./figure/DataSet_mmiss.root");
    //TFile *f0 = TFile::Open("./figure/DataSet_mhad.root");
    RooDataSet *ds0 = (RooDataSet*) f0->Get("ds_data");
    //debug
    /*TCanvas* c = new TCanvas("c","c", 600, 600 );
    RooPlot* rpl = mhad.frame();
    ds0->plotOn(rpl);
    rpl->Draw();
    c->cd(); 
    cout<<__LINE__<<endl;   
    */
    double smooth = 0.6, nSigma = 1;// 0.6 for Umissfit, 0.9 for mhad
    TChain *ch1 = new TChain("sig");
    ch1->Add("figure/sig.root");
    RooDataSet *sigshape  = new RooDataSet("sigshape", "sigshape", ch1,Umissfit);
    RooKeysPdf sig("sig", "", Umissfit, *sigshape, 0,1);

    RooRealVar mean("mean", "",  1e-4,-1e-4, 1e-4);// (1e-3,-2e-3, 3e-3) Umissfit, (1e-3,-2e-3, 2e-3) mhad 
    RooRealVar sigma("sigma","", 1e-3, 0, 3e-3);   // (1e-3, 0, 3e-3) Umissfit, (1e-3, 0, 2e-3) mhad
    RooGaussian gauss("gauss", "", Umissfit,mean,sigma);
    RooFFTConvPdf psig("psig","",Umissfit,sig,gauss);
        
    TChain* ch_bkg = new TChain("bkg");
    ch_bkg->Add("figure/offical_bkg.root");
    RooDataSet *bkgshape = new RooDataSet("bkgshape",  "bkgshape",  ch_bkg, Umissfit);
    RooKeysPdf pbkg("pbkg","", Umissfit,*bkgshape, nSigma, smooth);

    TChain* ch_bkg_k3pi = new TChain("bkg");
    ch_bkg_k3pi->Add("figure/k3pi_bkg.root");
    RooDataSet *bkgshape_k3pi = new RooDataSet("bkgshape_k3pi", "bkgshape_k3pi", ch_bkg_k3pi, Umissfit);
    RooKeysPdf pbkg_k3pi("pbkg_k3pi","", Umissfit, *bkgshape_k3pi, nSigma, smooth);

    double bfmin = 0, bfmax = 1e-3;
    RooRealVar bf("bf", "", 4e-4, bfmin, bfmax);//bf
    RooRealVar nDz("nDz", "", 2367990);       //Ntag
    RooRealVar sigeff("sigeff", "", 0.126171);//eff
    RooFormulaVar nsig("nsig","@0*@1*@2",RooArgList(nDz, sigeff, bf));
    RooRealVar nbkg("nbkg", "", 90, 40, 150);
    RooRealVar nbkg_k3pi("nbkg_k3pi", "", 20.1);//3.6+7.7+8.8
    RooAddPdf allpdf("allpdf", "", RooArgList(psig,pbkg,pbkg_k3pi), RooArgList(nsig,nbkg,nbkg_k3pi));  

    RooFitResult* res = allpdf.fitTo(*ds0, RooFit::Save(), Extended(kTRUE));
    res->Print("v");
    double nll0 = res->minNll(); 
    
    double Nsig = nsig.getVal(), Nsig_err = bf.getError() * nDz.getVal() * sigeff.getVal();  
    double n_err = sqrt(( nbkg.getError())**2);

    RooPlot* rp = Umissfit.frame();
    ds0->plotOn(rp, MarkerStyle(20), MarkerSize(1.0),MarkerColor(kBlack));
    allpdf.plotOn(rp,LineStyle(1),  LineColor(kBlue),   LineWidth(2));
    double chisq = rp->chiSquare(2);
    allpdf.plotOn(rp,Components(psig),LineStyle(kDashed),LineColor(kRed), LineWidth(2));
    allpdf.plotOn(rp,Components(pbkg),LineStyle(kDashed),LineColor(kBlack), LineWidth(2));
    allpdf.plotOn(rp,Components(pbkg_k3pi),LineStyle(kDashed),LineColor(kGreen), LineWidth(2));
    
    ofstream out("figure/para_mmiss.txt");
    out<<mean.getVal()<<"\n";
    out<<sigma.getVal()<<"\n";

    TCanvas* c0 = new TCanvas("cvs","c", 600, 600 );
    rp->Draw();
    
    TLatex l0;
    l0.SetNDC();
    l0.SetTextSize(0.02);
    l0.SetTextSize(0.03);
    l0.DrawLatex(0.55, 0.83,  Form("nsig = %.1f #pm %.1f",Nsig,Nsig_err));
    l0.DrawLatex(0.55, 0.78,  Form("nbkg = %.1f #pm %.1f",nbkg.getVal()+nbkg_k3pi.getVal(),n_err));

    c0->Modified();
    c0->cd();
   
    cout<<"mean: "<<mean.getVal()<<"+/-"<<mean.getError()<<endl;
    cout<<"sigma: "<<sigma.getVal()<<"+/-"<<sigma.getError()<<endl;
    //bf.setVal(0); bf.setConstant();
    //RooFitResult* res2 = allpdf.fitTo(*ds0, RooFit::Save(), RooFit::Minos(1));
    //double nll1 = res2->minNll();
    //cout<<"nll0: "<<nll0<<" nll1: "<<nll1<<endl;
    //cout<<"significance: "<< sqrt(2*(nll1-nll0))<<endl;
    cout<<"chisq: "<<chisq<<endl;

}
