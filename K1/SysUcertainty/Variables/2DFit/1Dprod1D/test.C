using std::map;
using namespace RooFit;
int test(){
TStyle *bes3Style= new TStyle("bes3","bes3 style");
Int_t icol=0;
bes3Style->SetFrameBorderMode(icol);
bes3Style->SetCanvasBorderMode(icol);
bes3Style->SetPadBorderMode(icol);
bes3Style->SetPadColor(icol);
bes3Style->SetCanvasColor(icol);
bes3Style->SetStatColor(icol);
bes3Style->SetTitleFillColor(icol);
bes3Style->SetPalette(1);
bes3Style->SetPaperSize(TStyle::kUSLetter);
bes3Style->SetPadTopMargin(.12);
bes3Style->SetPadLeftMargin(.15);
bes3Style->SetPadRightMargin(.08);
bes3Style->SetPadBottomMargin(.15);
Int_t font=132;
Double_t tsize=0.05;
bes3Style->SetTextFont(font);
bes3Style->SetTextSize(tsize);
bes3Style->SetLabelSize(tsize,"xyz");
bes3Style->SetLabelOffset(0.01,"xyz");
bes3Style->SetTitleFont(font,"xyz");
bes3Style->SetTitleSize(tsize,"xyz");
bes3Style->SetTitleXOffset(1.);
bes3Style->SetTitleYOffset(1.4);
bes3Style->SetTitleBorderSize(2.);
bes3Style->SetMarkerStyle(0);
bes3Style->SetMarkerSize(0.8);
bes3Style->SetFrameBorderMode (0.);
bes3Style->SetFrameLineWidth  (2.);
bes3Style->SetLineWidth(1.0);
bes3Style->SetHistLineWidth(2.);
bes3Style->SetLineStyleString(2,"[12 12]");
bes3Style->SetErrorX(0.001);
bes3Style->SetOptTitle(0);
bes3Style->SetOptStat(0); 
bes3Style->SetLineStyleString(2,"[30 10]");
bes3Style->SetLineStyleString(3,"[4 8]");
bes3Style->SetLineStyleString(4,"[15 12 4 12]");
bes3Style->SetLineStyleString(5,"[15 15]");
bes3Style->SetLineStyleString(6,"[15 12 4 12 4 12]");
bes3Style->SetOptDate(0);
bes3Style->SetDateY(.98);
bes3Style->SetStripDecimals(kFALSE);
bes3Style->SetEndErrorSize(0.0);
gROOT->SetStyle("bes3");
gROOT->ForceStyle();
gROOT->Reset();
using namespace RooFit;


    TCut m0 ="mode==0";

    RooRealVar mhad("mhad","M_{K^{-}#pi^{+}#pi^{-}}  (GeV/c^2)", 0.9,1.68 );
    RooRealVar Umissfit("Umissfit","U_{miss}  (GeV/c^2)", -0.1, 0.1 );
    RooArgSet m_u( mhad, Umissfit);
    TChain ch("tagD");
    ch.Add("figure/data.root");

    TTree* t0 = ch.CopyTree(m0);

    double smooth = 1.0, nSigma = 0;
    double smooth_k3pi = 1.0, smooth_bkg=1.0;
    //mode0
    double var0=0, var1=0, var2=0, var3=0;
    ifstream par_mmiss("./figure/para_mmiss2.txt");
    par_mmiss>>var0>>var1;
    ifstream par_mhad("./figure/para_mhad.txt");
    par_mhad>>var2>>var3;

    cout<<__LINE__<<" para: "<<var0<<"\t"<<var1<<","<<"\t"<<var2<<"\t"<<var3<<endl;
    RooRealVar mean_mmiss ("mean_mmiss" , "" ,var0);
    RooRealVar mean_mhad  ("mean_mhad"  , "" ,var2);
    RooRealVar sigma_mmiss("sigma_mmiss", "" ,var1);
    RooRealVar sigma_mhad ("sigma_mhad" , "" ,var3);

    RooGaussian gauss_mmiss("gauss_mmiss","",Umissfit, mean_mmiss, sigma_mmiss);
    RooGaussian gauss_mhad ("gauss_mhad", "",mhad,     mean_mhad,  sigma_mhad );

    TChain* ch_sig0 = new TChain("sig");
    ch_sig0->Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/CWR_VSS/update/pdf/dat/sig_0.root");
    RooDataSet *sigshape0_mmiss = new RooDataSet("sigshape0_mmiss", "sigshape0_mmiss", ch_sig0, Umissfit);
    RooDataSet *sigshape0_mhad  = new RooDataSet("sigshape0_mhad",  "sigshape0_mhad",  ch_sig0, mhad);
    RooKeysPdf psig0_mmiss("psig0_mmiss","",Umissfit,*sigshape0_mmiss, nSigma, smooth);
    RooKeysPdf psig0_mhad ("psig0_mhad", "",mhad,    *sigshape0_mhad,  nSigma, smooth);
    RooFFTConvPdf psig0_mmiss_g("psig0_mmiss_g","",Umissfit,psig0_mmiss,gauss_mmiss);
    RooFFTConvPdf psig0_mhad_g("psig0_mhad_g","",mhad,psig0_mhad,gauss_mhad);
    RooProdPdf psig0("psig0","", RooArgSet(psig0_mmiss_g, psig0_mhad_g));
        
    ////test
    //RooPlot *rp  = Umissfit.frame();
    //RooPlot *rp1 = mhad.frame();
    ////sigshape0_mmiss->plotOn(rp);
    ////sigshape0_mhad->plotOn(rp1);
    //psig0_mmiss_g.plotOn(rp, Project(mhad));
    //psig0_mhad_g.plotOn(rp1, Project(Umissfit), LineColor(2));
    //psig0.plotOn(rp1, Project(Umissfit), LineColor(2));
    //psig0.plotOn(rp, Project(Umissfit), LineColor(2));
    //rp->Draw(); 
    //rp1->Draw(); 
    //break;
    ////test

    TChain* ch_bkg0 = new TChain("bkg");
    ch_bkg0->Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/CWR_VSS/update/pdf/dat/offical_bkg_0.root");
    RooDataSet *bkgshape0_mmiss = new RooDataSet("bkgshape0_mmiss", "bkgshape0_mmiss", ch_bkg0, Umissfit);
    RooDataSet *bkgshape0 _mhad = new RooDataSet("bkgshape0_mhad",  "bkgshape0_mhad",  ch_bkg0, mhad);
    RooKeysPdf pbkg0_mmiss("pbkg0_mmiss","",Umissfit,*bkgshape0_mmiss, nSigma, smooth_bkg);
    RooKeysPdf pbkg0_mhad("pbkg0_mhad", "",mhad,    *bkgshape0_mhad,  nSigma, smooth_bkg);
    RooProdPdf pbkg0("pbkg0","", RooArgSet(pbkg0_mmiss, pbkg0_mhad));

    TChain* ch_bkg0_k3pi = new TChain("bkg");
    ch_bkg0_k3pi->Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/CWR_VSS/update/pdf/dat/k3pi_bkg_0.root");
    RooDataSet *bkgshape0_k3pi_mmiss = new RooDataSet("bkgshape0_k3pi_mmiss", "bkgshape0_k3pi_mmiss", ch_bkg0_k3pi,  Umissfit);
    RooDataSet *bkgshape0_k3pi_mhad  = new RooDataSet("bkgshape0_k3pi_mhad",  "bkgshape0_k3pi_mhad" , ch_bkg0_k3pi,  mhad);
    RooKeysPdf pbkg0_k3pi_mmiss("pbkg0_k3pi_mmiss","",Umissfit, *bkgshape0_k3pi_mmiss, nSigma, smooth_k3pi);
    RooKeysPdf pbkg0_k3pi_mhad("pbkg0_k3pi_mhad",  "",mhad,     *bkgshape0_k3pi_mhad,  nSigma, smooth_k3pi);
    RooProdPdf pbkg0_k3pi("pbkg0_k3pi","", RooArgSet(pbkg0_k3pi_mmiss, pbkg0_k3pi_mhad));


  double eff[4];
  for(int i=0;i<4;i++){
  ifstream dat("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/CWR_VSS/update/eff/dat/eff.txt");
  dat>>eff[0]>>eff[1]>>eff[2]>>eff[3];
  }
  double Ntag[3];
  for(int i=0;i<3;i++){
  ifstream ntag("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/mBC/FitData/figure/read.txt");
  ntag>>Ntag[0]>>Ntag[1]>>Ntag[2];
  }
    double bfmin = 0, bfmax = 50;
    RooRealVar bf("bf", "", 25e-4, bfmin, bfmax);//e-5
    RooRealVar nDz_0("nDz_0", "", Ntag[0]);//e+05 data
    RooRealVar sigeff_0("sigeff_0", "", eff[0]);

    RooFormulaVar nsig0("nsig0", "@0*@1*@2", RooArgList(nDz_0, sigeff_0, bf));
    RooRealVar nbkg0("nbkg0", "", 20,1,50); 
    
    double mis[3];
    for(int i=0;i<3;i++){
    ifstream nmis("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/CWR_VSS/update/misidentify/dat/misidentify103.dat");
    nmis>>mis[0]>>mis[1]>>mis[2]; }
    RooRealVar nbkg0_k3pi("nbkg0_k3pi", "", mis[0]);
    
    RooAddPdf allpdf0("allpdf0", "", RooArgList(psig0,pbkg0,pbkg0_k3pi ), RooArgList(nsig0,nbkg0,nbkg0_k3pi));
    //test
    //RooPlot *rp = Umissfit.frame();
    //RooPlot *rp1= mhad.frame();
    //allpdf0.plotOn(rp);
    //rp->Draw();
    //break;
    //test
    RooDataSet* ds0 = new RooDataSet("ds0", "", t0, m_u);
    RooCategory sample("sample","sample") ;
    std::map<std::string,RooDataSet*> dsmap;
    sample.defineType("mode0"); dsmap["mode0"] = ds0;

    RooDataSet combData("combData","combined data",m_u, RooFit::Index(sample),RooFit::Import(dsmap));
    RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
    simPdf.addPdf(allpdf0, "mode0");//(pdf,cat_lable), pdf to be added; name of the category state to be associated to be pdf

    RooFitResult* res = simPdf.fitTo(combData, RooFit::Save(), RooFit::Minos(1));
    res->Print("v");
    double minNll = -(res->minNll());//return minimized -log(L) value
    
    double Nsig0 = nsig0.getVal(), Nsig0_err = bf.getError() * nDz_0.getVal() * sigeff_0.getVal();
    double n0_err = sqrt(( nbkg0.getError())**2);        

    //bf.setVal(0); bf.setConstant();
    //RooFitResult* res2 = simPdf.fitTo(combData, RooFit::Save(), RooFit::Minos(1));
    //minNll += res2->minNll();
    //cout<<"Statisical significance: "<<sqrt(2*minNll)<<endl;
    
     cout<<Nsig0<<"+/-"<<Nsig0_err<<endl;;
/*
    TCanvas *cv = new TCanvas("cv", "", 800,1200);
    cv->Divide(2,3);
    cv->cd(1);
*/

//////////////////////////////////////////////plot and draw///////////////////////////////////////////////
    TCanvas *c1 = new TCanvas("c1", "c1",0,0,600,600);
    c1->Range(0,0,1,1);
    c1->SetBorderSize(2);
    c1->SetLeftMargin(0.15);
    c1->SetRightMargin(0.05);
    c1->SetTopMargin(0.05);
    c1->SetBottomMargin(0.05);
    c1->SetFrameFillColor(0);
    /////////////////////////////Y axis//////////
    leg = new TLegend(0.01,0.47,0.055,0.60);
    leg->SetTextAlign(22);
    leg->SetTextAlign(22);
    leg->SetTextFont(132);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->SetTextColor(1);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetLineWidth(0);
    leg->SetTextAngle(90);
    leg->SetHeader("Events (GeV^{2}/c^{4})");
    leg->Draw();
    /////////////////////////Y axis///////////
    leg = new TLegend(0.51,0.47,0.56,0.60);
    leg->SetTextAlign(22);
    leg->SetTextAlign(22);
    leg->SetTextFont(132);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->SetTextColor(1);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetLineWidth(0);
    leg->SetTextAngle(90);
    leg->SetHeader("Events (GeV/c^{2})");
    leg->Draw();
    /////////////////////////X axis///////////
    leg = new TLegend(0.2,0.015,0.4,0.055);
    leg->SetTextAlign(22);
    leg->SetTextAlign(22);
    leg->SetTextFont(132);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->SetTextColor(1);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetLineWidth(0);
    leg->SetTextAngle(0);
    leg->SetHeader("M_{miss}^{2} (GeV^{2}/c^{4})");
    leg->Draw();
    //////////////////////X axis////////////
    leg = new TLegend(0.66,0.015,0.9,0.055);
    leg->SetTextAlign(22);
    leg->SetTextFont(132);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->SetTextColor(1);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetLineWidth(0);
    leg->SetTextAngle(0);
    leg->SetHeader("M_{K^{-}#pi^{+}#pi^{-}} (GeV/c^{2})");
    leg->Draw();
    ///////////////////////mode0//////////
    RooPlot *frame = Umissfit.frame(40);
    RooPlot *frame1 = mhad.frame(26);
//    frame0->SetTitle("mode0");
    combData.plotOn(frame,RooFit::Cut("sample==sample::mode0")) ;
    simPdf.plotOn(frame, RooFit::Slice(sample,"mode0"), RooFit::ProjWData(sample, combData));
//    simPdf.plotOn(frame0, RooFit::Slice(sample,"mode0"), RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(*psig0,*pomegasig)),RooFit::LineStyle(kDashed), RooFit::LineColor(kBlue));
    simPdf.plotOn(frame, RooFit::Slice(sample,"mode0"), RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(psig0)),RooFit::LineStyle(kDashed+1), RooFit::LineColor(kRed));
    simPdf.plotOn(frame, RooFit::Slice(sample,"mode0"), RooFit::ProjWData(sample, combData), RooFit::Components(pbkg0),RooFit::LineStyle(kDashed), RooFit::LineColor(kBlack));
    simPdf.plotOn(frame, RooFit::Slice(sample,"mode0"), RooFit::ProjWData(sample, combData), RooFit::Components(pbkg0_k3pi),RooFit::LineStyle(kDashed), RooFit::LineColor(kCyan-3));
//    simPdf.plotOn(frame0, RooFit::Slice(sample,"mode0"), RooFit::ProjWData(sample, combData), RooFit::Components(*pomegasig),RooFit::LineStyle(kDashed), RooFit::LineColor(kGreen+3));
    frame->GetYaxis()->SetRangeUser(0,15);
    TPad *c1_1 = new TPad("c1_1", "c1_1",0.02,0.62,0.52,1.0);
    c1_1->Draw();
    c1_1->cd();
    c1_1->Range(0,0,1,1);
    c1_1->SetFillStyle(4000);
    c1_1->SetBorderMode(0);
    c1_1->SetBorderSize(2);
    c1_1->SetTickx();
    c1_1->SetTicky();
    c1_1->SetLeftMargin(0.20);
    c1_1->SetRightMargin(0.05);
    c1_1->SetTopMargin(0.1);
    c1_1->SetBottomMargin(0.15);
    c1_1->SetFrameFillColor(0);
    
    frame->Draw();
    frame->GetYaxis()->SetNdivisions(505);
    frame->GetXaxis()->SetNdivisions(5);
    frame->GetYaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetNoExponent(kFALSE);
    frame->GetYaxis()->SetTitleSize(0);
    frame->GetXaxis()->SetTitleSize(0);
    frame->GetXaxis()->SetLabelSize(0.0);
    
    TGaxis *axis1 = new TGaxis(-0.1,0,-0.1,12,0,12,510,"");
    axis1->SetName("axis1");
    axis1->SetLabelSize(0.1);
    axis1->SetLabelOffset(0.01);
    axis1->SetNdivisions(505);
    axis1->Draw();

 
    c1_1->Modified();
    c1->cd();
  
    ///////////
    combData.plotOn(frame1,RooFit::Cut("sample==sample::mode0")) ;
    simPdf.plotOn(frame1, RooFit::Slice(sample,"mode0"), RooFit::ProjWData(sample, combData));
    simPdf.plotOn(frame1, RooFit::Slice(sample,"mode0"), RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(psig0)),RooFit::LineStyle(kDashed+1), RooFit::LineColor(kRed));
    simPdf.plotOn(frame1, RooFit::Slice(sample,"mode0"), RooFit::ProjWData(sample, combData), RooFit::Components(pbkg0),RooFit::LineStyle(kDashed), RooFit::LineColor(kBlack));
    simPdf.plotOn(frame1, RooFit::Slice(sample,"mode0"), RooFit::ProjWData(sample, combData), RooFit::Components(pbkg0_k3pi),RooFit::LineStyle(kDashed), RooFit::LineColor(kCyan-3));
    frame1->GetYaxis()->SetRangeUser(0,12.5);
    
   c1_2 = new TPad("c1_2", "c1_2",0.5,0.62,1,1.0);
   c1_2->Draw();
   c1_2->cd();
   c1_2->Range(0,0,1,1);
   c1_2->SetFillStyle(4000);
   c1_2->SetBorderMode(0);
   c1_2->SetBorderSize(2);
   c1_2->SetTickx();
   c1_2->SetTicky();
   c1_2->SetLeftMargin(0.20);
   c1_2->SetRightMargin(0.05);
   c1_2->SetTopMargin(0.1);
   c1_2->SetBottomMargin(0.15);
   c1_2->SetFrameFillColor(0);

   frame1->Draw();
   frame1->GetYaxis()->SetNdivisions(505);
   frame1->GetXaxis()->SetNdivisions(5);
   frame1->GetYaxis()->SetLabelSize(0.0);
   frame1->GetYaxis()->SetMoreLogLabels(2);
   frame1->GetYaxis()->SetNoExponent(kFALSE);
   frame1->GetYaxis()->SetTitleSize(0);
   frame1->GetXaxis()->SetTitleSize(0);
   frame1->GetXaxis()->SetLabelSize(0.0);

    TGaxis *axis1 = new TGaxis(0.9,0,0.9,13,0,13,510,"");
    axis1->SetName("axis1");
    axis1->SetLabelSize(0.1);
    axis1->SetLabelOffset(0.01);
    axis1->SetNdivisions(505);
    axis1->Draw();

    TLatex l0;
    l0.SetNDC();
    l0.SetTextSize(0.035);
    l0.SetTextSize(0.06);
    l0.DrawLatex(0.55, 0.83,  Form("nsig = %.1f #pm %.1f",Nsig0,Nsig0_err));
    l0.DrawLatex(0.55, 0.78,  Form("nbkg = %.1f #pm %.1f",nbkg0.getVal()+nbkg0_k3pi.getVal(),n0_err));
 
    c1_2->Modified();
    c1->cd();

}
