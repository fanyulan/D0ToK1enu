using std::map;
using namespace RooFit;
int sim_smear(){
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
    TCut m1 ="mode==1";
    TCut m3 ="mode==3";

    RooRealVar mhad("mhad","M_{K^{-}#pi^{+}#pi^{-}}  (GeV/c^2)", 0.9,1.68 );
    RooRealVar Umissfit("Umissfit","U_{miss}  (GeV/c^2)", -0.1, 0.1 );
    RooArgSet m_u( mhad, Umissfit);
    TChain ch("tagD");
    ch.Add("figure/data.root");

    TTree* t0 = ch.CopyTree(m0);
    TTree* t1 = ch.CopyTree(m1);
    TTree* t3 = ch.CopyTree(m3);

    cout<<t0->GetEntries()<<t1->GetEntries()<<t3->GetEntries()<<endl;
    double smooth = 0.6, nSigma = 0;
    double smooth_k3pi = 1.0, smooth_bkg=1.0;
    //mode0
    double var0=0, var1=0, var2=0, var3=0;
    ifstream par_mmiss("./figure/para_mmiss.txt");
    par_mmiss>>var0>>var1;
    ifstream par_mhad("./figure/para_mhad.txt");
    par_mhad>>var2>>var3;

    cout<<__LINE__<<" para: "<<var0<<"\t"<<var1<<","<<"\t"<<var2<<"\t"<<var3<<endl;
    RooRealVar mean_mmiss ("mean_mmiss" , "" ,var0);
    RooRealVar mean_mhad  ("mean_mhad"  , "" ,var2);
    RooRealVar sigma_mmiss("sigma_mmiss", "" ,var1);
    RooRealVar sigma_mhad ("sigma_mhad" , "" ,var3);

    RooGaussian gauss_mmiss("gauss_mmiss","",Umissfit, mean_mmiss, sigma_mmiss);
    RooGaussian gauss_mhad ("gauss_mhad","", mhad, mean_mhad, sigma_mhad);

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
    //RooPlot *rp = Umissfit.frame();
    //RooPlot *rp1= mhad.frame();
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
    //mode1 
    TChain* ch_sig1 = new TChain("sig");
    ch_sig1->Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/CWR_VSS/update/pdf/dat/sig_1.root");
    RooDataSet *sigshape1_mmiss = new RooDataSet("sigshape1_mmiss", "sigshape1_mmiss", ch_sig1, Umissfit);
    RooDataSet *sigshape1_mhad  = new RooDataSet("sigshape1_mhad",  "sigshape1_mhad",  ch_sig1, mhad);
    RooKeysPdf psig1_mmiss("psig1_mmiss","",Umissfit,*sigshape1_mmiss, nSigma, smooth);
    RooKeysPdf psig1_mhad ("psig1_mhad", "",mhad,    *sigshape1_mhad,  nSigma, smooth);
    RooFFTConvPdf psig1_mmiss_g("psig1_mmiss_g","",Umissfit,psig1_mmiss,gauss_mmiss);
    RooFFTConvPdf psig1_mhad_g ("psig1_mhad_g","",mhad,     psig1_mhad, gauss_mhad);
    RooProdPdf psig1("psig1","", RooArgSet(psig1_mmiss_g, psig1_mhad_g));

    TChain* ch_bkg1 = new TChain("bkg");
    ch_bkg1->Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/CWR_VSS/update/pdf/dat/offical_bkg_1.root");
    RooDataSet *bkgshape1_mmiss = new RooDataSet("bkgshape1_mmiss", "bkgshape1_mmiss", ch_bkg1, Umissfit);
    RooDataSet *bkgshape1 _mhad = new RooDataSet("bkgshape1_mhad",  "bkgshape1_mhad",  ch_bkg1, mhad);
    RooKeysPdf pbkg1_mmiss("pbkg1_mmiss","",Umissfit,*bkgshape1_mmiss, nSigma, smooth_bkg);
    RooKeysPdf pbkg1_mhad("pbkg1_mhad", "" ,mhad,    *bkgshape1_mhad,  nSigma, smooth_bkg);
    RooProdPdf pbkg1("pbkg1","", RooArgSet(pbkg1_mmiss, pbkg1_mhad));

    TChain* ch_bkg1_k3pi = new TChain("bkg");
    ch_bkg1_k3pi->Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/CWR_VSS/update/pdf/dat/k3pi_bkg_1.root");
    RooDataSet *bkgshape1_k3pi_mmiss = new RooDataSet("bkgshape1_k3pi_mmiss", "bkgshape1_k3pi_mmiss", ch_bkg1_k3pi, Umissfit);
    RooDataSet *bkgshape1_k3pi_mhad  = new RooDataSet("bkgshape1_k3pi_mhad",  "bkgshape1_k3pi_mhad" , ch_bkg1_k3pi, mhad);
    RooKeysPdf pbkg1_k3pi_mmiss("pbkg1_k3pi_mmiss","",Umissfit, *bkgshape1_k3pi_mmiss, nSigma, smooth_k3pi);
    RooKeysPdf pbkg1_k3pi_mhad("pbkg1_k3pi_mhad",  "",mhad,     *bkgshape1_k3pi_mhad,  nSigma, smooth_k3pi);
    RooProdPdf pbkg1_k3pi("pbkg1_k3pi","", RooArgSet(pbkg1_k3pi_mmiss, pbkg1_k3pi_mhad));
    //mode3
    TChain* ch_sig3 = new TChain("sig");
    ch_sig3->Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/CWR_VSS/update/pdf/dat/sig_2.root");
    RooDataSet *sigshape3_mmiss = new RooDataSet("sigshape3_mmiss", "sigshape3_mmiss", ch_sig3, Umissfit);
    RooDataSet *sigshape3_mhad  = new RooDataSet("sigshape3_mhad",  "sigshape3_mhad",  ch_sig3, mhad);
    RooKeysPdf psig3_mmiss("psig3_mmiss","",Umissfit,*sigshape3_mmiss,  nSigma, smooth);
    RooKeysPdf psig3_mhad ("psig3_mhad", "",mhad,    *sigshape3_mhad,   nSigma, smooth);
    RooFFTConvPdf psig3_mmiss_g("psig3_mmiss_g","",Umissfit, psig3_mmiss,gauss_mmiss);
    RooFFTConvPdf psig3_mhad_g ("psig3_mhad_g","", mhad,     psig3_mhad, gauss_mhad);
    RooProdPdf psig3("psig3","", RooArgSet(psig3_mmiss_g, psig3_mhad_g));

    TChain* ch_bkg3 = new TChain("bkg");
    ch_bkg3->Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/CWR_VSS/update/pdf/dat/offical_bkg_2.root");
    RooDataSet *bkgshape3_mmiss = new RooDataSet("bkgshape3_mmiss", "bkgshape3_mmiss", ch_bkg3, Umissfit);
    RooDataSet *bkgshape3 _mhad = new RooDataSet("bkgshape3_mhad",  "bkgshape3_mhad",  ch_bkg3, mhad);
    RooKeysPdf pbkg3_mmiss("pbkg3_mmiss","",Umissfit,*bkgshape3_mmiss, nSigma, smooth_bkg);
    RooKeysPdf pbkg3_mhad("pbkg3_mhad", "",mhad,    *bkgshape3_mhad,  nSigma, smooth_bkg);
    RooProdPdf pbkg3("pbkg3","", RooArgSet(pbkg3_mmiss, pbkg3_mhad));

    TChain* ch_bkg3_k3pi = new TChain("bkg");
    ch_bkg3_k3pi->Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/CWR_VSS/update/pdf/dat/k3pi_bkg_2.root");
    RooDataSet *bkgshape3_k3pi_mmiss = new RooDataSet("bkgshape3_k3pi_mmiss", "bkgshape3_k3pi_mmiss", ch_bkg3_k3pi, Umissfit);
    RooDataSet *bkgshape3_k3pi_mhad  = new RooDataSet("bkgshape3_k3pi_mhad",  "bkgshape3_k3pi_mhad" , ch_bkg3_k3pi,  mhad);
    RooKeysPdf pbkg3_k3pi_mmiss("pbkg3_k3pi_mmiss","",Umissfit, *bkgshape3_k3pi_mmiss, nSigma, smooth_k3pi);
    RooKeysPdf pbkg3_k3pi_mhad("pbkg3_k3pi_mhad",  "",mhad,     *bkgshape3_k3pi_mhad,  nSigma, smooth_k3pi);
    RooProdPdf pbkg3_k3pi("pbkg3_k3pi","", RooArgSet(pbkg3_k3pi_mmiss, pbkg3_k3pi_mhad));


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
    RooRealVar bf("bf", "", 25, bfmin, bfmax);//e-5
    RooRealVar nDz_0("nDz_0", "", Ntag[0]);//e+05 data
    RooRealVar nDz_1("nDz_1", "", Ntag[1]);
    RooRealVar nDz_3("nDz_3", "", Ntag[2]);
    RooRealVar sigeff_0("sigeff_0", "", eff[0]);
    RooRealVar sigeff_1("sigeff_1", "", eff[1]);
    RooRealVar sigeff_3("sigeff_3", "", eff[2]);
    RooRealVar sigratio("sigratio", "", eff[3]);

    RooFormulaVar nsig0("nsig0", "@0*@1*@2", RooArgList(nDz_0, sigeff_0, bf));
    RooFormulaVar nsig1("nsig1", "@0*@1*@2", RooArgList(nDz_1, sigeff_1, bf));
    RooFormulaVar nsig3("nsig3", "@0*@1*@2", RooArgList(nDz_3, sigeff_3, bf));
    RooRealVar nbkg0("nbkg0", "", 20,1,50); 
    RooRealVar nbkg1("nbkg1", "", 40,1,100);
    RooRealVar nbkg3("nbkg3", "", 30,1,100);
    
    double mis[3];
    for(int i=0;i<3;i++){
    ifstream nmis("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/CWR_VSS/update/misidentify/dat/misidentify103.dat");
    nmis>>mis[0]>>mis[1]>>mis[2]; }
    RooRealVar nbkg0_k3pi("nbkg0_k3pi", "", mis[0]);
    RooRealVar nbkg1_k3pi("nbkg1_k3pi", "", mis[1]);
    RooRealVar nbkg3_k3pi("nbkg3_k3pi", "", mis[2]);
    
    RooAddPdf allpdf0("allpdf0", "", RooArgList(psig0,pbkg0,pbkg0_k3pi ), RooArgList(nsig0,nbkg0,nbkg0_k3pi));
    RooAddPdf allpdf1("allpdf1", "", RooArgList(psig1,pbkg1,pbkg1_k3pi ), RooArgList(nsig1,nbkg1,nbkg1_k3pi));
    RooAddPdf allpdf3("allpdf3", "", RooArgList(psig3,pbkg3,pbkg3_k3pi ), RooArgList(nsig3,nbkg3,nbkg3_k3pi));

    RooDataSet* ds0 = new RooDataSet("ds0", "", t0, m_u);
    RooDataSet* ds1 = new RooDataSet("ds1", "", t1, m_u);
    RooDataSet* ds3 = new RooDataSet("ds3", "", t3, m_u);
    RooCategory sample("sample","sample") ;
    std::map<std::string,RooDataSet*> dsmap;
    sample.defineType("mode0"); dsmap["mode0"] = ds0;
    sample.defineType("mode1"); dsmap["mode1"] = ds1;
    sample.defineType("mode3"); dsmap["mode3"] = ds3;

    RooDataSet combData("combData","combined data",m_u, RooFit::Index(sample),RooFit::Import(dsmap));
    RooSimultaneous simPdf("simPdf","simultaneous pdf",sample) ;
    simPdf.addPdf(allpdf0, "mode0");//(pdf,cat_lable), pdf to be added; name of the category state to be associated to be pdf
    simPdf.addPdf(allpdf1, "mode1");
    simPdf.addPdf(allpdf3, "mode3");

    RooFitResult* res = simPdf.fitTo(combData, RooFit::Save(), RooFit::Minos(1));
    res->Print("v");
    double minNll = -(res->minNll());//return minimized -log(L) value
    
    double Nsig0 = nsig0.getVal(), Nsig0_err = bf.getError() * nDz_0.getVal() * sigeff_0.getVal();
    double Nsig1 = nsig1.getVal(), Nsig1_err = bf.getError() * nDz_1.getVal() * sigeff_1.getVal();
    double Nsig3 = nsig3.getVal(), Nsig3_err = bf.getError() * nDz_3.getVal() * sigeff_3.getVal();
    double n0_err = sqrt(( nbkg0.getError())**2);        
    double n1_err = sqrt(( nbkg1.getError())**2);        
    double n3_err = sqrt(( nbkg3.getError())**2);        

    //bf.setVal(0); bf.setConstant();
    //RooFitResult* res2 = simPdf.fitTo(combData, RooFit::Save(), RooFit::Minos(1));
    //minNll += res2->minNll();
    //cout<<"Statisical significance: "<<sqrt(2*minNll)<<endl;
    
     cout<<Nsig0<<"+/-"<<Nsig0_err<<endl;;
     cout<<Nsig1<<"+/-"<<Nsig1_err<<endl;;
     cout<<Nsig3<<"+/-"<<Nsig3_err<<endl;;
     cout<<"total fit: "<<Nsig0+Nsig1+Nsig3<<"+/-"<<Nsig1_err/Nsig1*(Nsig0+Nsig1+Nsig3)<<endl;
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
    frame->GetYaxis()->SetNdivisions(4);
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
    axis1->SetNdivisions(4);
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
   frame1->GetYaxis()->SetNdivisions(4);
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
    axis1->SetNdivisions(4);
    axis1->Draw();

    TLatex l0;
    l0.SetNDC();
    l0.SetTextSize(0.035);
    l0.SetTextSize(0.06);
    l0.DrawLatex(0.55, 0.83,  Form("nsig = %.1f #pm %.1f",Nsig0,Nsig0_err));
    l0.DrawLatex(0.55, 0.78,  Form("nbkg = %.1f #pm %.1f",nbkg0.getVal()+nbkg0_k3pi.getVal(),n0_err));
 
    c1_2->Modified();
    c1->cd();

    ///////////////////////mode1//////////////
    RooPlot *frame = Umissfit.frame(40);
    RooPlot *frame1 = mhad.frame(26);

    combData.plotOn(frame,RooFit::Cut("sample==sample::mode1")) ;
    simPdf.plotOn(frame, RooFit::Slice(sample,"mode1"), RooFit::ProjWData(sample, combData));
    simPdf.plotOn(frame, RooFit::Slice(sample,"mode1"), RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(psig1)),RooFit::LineStyle(kDashed+1), RooFit::LineColor(kRed));
    simPdf.plotOn(frame, RooFit::Slice(sample,"mode1"), RooFit::ProjWData(sample, combData), RooFit::Components(pbkg1),RooFit::LineStyle(kDashed), RooFit::LineColor(kBlack));
    simPdf.plotOn(frame, RooFit::Slice(sample,"mode1"), RooFit::ProjWData(sample, combData), RooFit::Components(pbkg1_k3pi),RooFit::LineStyle(kDashed), RooFit::LineColor(kCyan-3));
    frame->GetYaxis()->SetRangeUser(0,30);

   c1_3 = new TPad("c1_3", "c1_3",0.02,0.33,0.52,0.71);
   c1_3->Draw();
   c1_3->cd();
   c1_3->Range(0,0,1,1);
   c1_3->SetFillStyle(4000);
   c1_3->SetBorderMode(0);
   c1_3->SetBorderSize(2);
   c1_3->SetTickx();
   c1_3->SetTicky();
   c1_3->SetLeftMargin(0.20);
   c1_3->SetRightMargin(0.05);
   c1_3->SetTopMargin(0.1);
   c1_3->SetBottomMargin(0.15);
   c1_3->SetFrameFillColor(0);

   frame->Draw();
   frame->GetYaxis()->SetNdivisions(504);
   frame->GetXaxis()->SetNdivisions(5);
   frame->GetYaxis()->SetLabelSize(0.0);
   frame->GetYaxis()->SetNoExponent(kFALSE);
   frame->GetYaxis()->SetTitleSize(0);
   frame->GetXaxis()->SetTitleSize(0);
   frame->GetXaxis()->SetLabelSize(0.0);

   TGaxis *axis1 = new TGaxis(-0.1,0,-0.1,28,0,28,510,"");
   axis1->SetName("axis1");
   axis1->SetLabelSize(0.1);
   axis1->SetLabelOffset(0.01);
   axis1->SetNdivisions(4);
   axis1->Draw();

   c1_3->Modified();
   c1->cd();
 
    combData.plotOn(frame1,RooFit::Cut("sample==sample::mode1")) ;
    simPdf.plotOn(frame1, RooFit::Slice(sample,"mode1"), RooFit::ProjWData(sample, combData));
    simPdf.plotOn(frame1, RooFit::Slice(sample,"mode1"), RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(psig1)),RooFit::LineStyle(kDashed+1), RooFit::LineColor(kRed));
    simPdf.plotOn(frame1, RooFit::Slice(sample,"mode1"), RooFit::ProjWData(sample, combData), RooFit::Components(pbkg1),RooFit::LineStyle(kDashed), RooFit::LineColor(kBlack));
    simPdf.plotOn(frame1, RooFit::Slice(sample,"mode1"), RooFit::ProjWData(sample, combData), RooFit::Components(pbkg1_k3pi),RooFit::LineStyle(kDashed), RooFit::LineColor(kCyan-3));

   frame1->GetYaxis()->SetRangeUser(0,24);
   TPad *c1_4 = new TPad("c1_4", "c1_4",0.5,0.33,1,0.71);
   c1_4->Draw();
   c1_4->cd();
   c1_4->Range(0,0,1,1);
   c1_4->SetFillStyle(4000);
   c1_4->SetBorderMode(0);
   c1_4->SetBorderSize(2);
   c1_4->SetTickx();
   c1_4->SetTicky();
   c1_4->SetLeftMargin(0.20);
   c1_4->SetRightMargin(0.05);
   c1_4->SetTopMargin(0.1);
   c1_4->SetBottomMargin(0.15);
   c1_4->SetFrameFillColor(0);
   
   frame1->Draw();
   frame1->GetYaxis()->SetNdivisions(504);
   frame1->GetXaxis()->SetNdivisions(5);
   frame1->GetYaxis()->SetLabelSize(0);
   frame1->GetYaxis()->SetNoExponent(kFALSE);
   frame1->GetYaxis()->SetTitleSize(0);
   frame1->GetXaxis()->SetTitleSize(0);
   frame1->GetXaxis()->SetLabelSize(0.0);

    TGaxis *axis1 = new TGaxis(0.9,0,0.9,24,0,24,510,"");
    axis1->SetName("axis1");
    axis1->SetLabelSize(0.1);
    axis1->SetLabelOffset(0.01);
    axis1->SetNdivisions(4);
    axis1->Draw();
    TLatex l1;
    l1.SetNDC();
    l1.SetTextSize(0.035);
    l1.SetTextSize(0.06);
    l1.DrawLatex(0.55, 0.83,  Form("nsig = %.1f #pm %.1f",Nsig1,Nsig1_err));
    l1.DrawLatex(0.55, 0.78,  Form("nbkg = %.1f #pm %.1f",nbkg1.getVal()+nbkg1_k3pi.getVal(),n1_err));
  
    c1_4->Modified();
    c1->cd();

    ////////////////////////mode3///////////////////
    RooPlot *frame = Umissfit.frame(40);
    RooPlot *frame1 = mhad.frame(26);
    
    combData.plotOn(frame,RooFit::Cut("sample==sample::mode3")) ;
    simPdf.plotOn(frame, RooFit::Slice(sample,"mode3"), RooFit::ProjWData(sample, combData));
    simPdf.plotOn(frame, RooFit::Slice(sample,"mode3"), RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(psig3)),RooFit::LineStyle(kDashed+1), RooFit::LineColor(kRed));
    simPdf.plotOn(frame, RooFit::Slice(sample,"mode3"), RooFit::ProjWData(sample, combData), RooFit::Components(pbkg3),RooFit::LineStyle(kDashed), RooFit::LineColor(kBlack));
    simPdf.plotOn(frame, RooFit::Slice(sample,"mode3"), RooFit::ProjWData(sample, combData), RooFit::Components(pbkg3_k3pi),RooFit::LineStyle(kDashed), RooFit::LineColor(kCyan-3));

    frame->GetYaxis()->SetRangeUser(0,20);

   c1_5 = new TPad("c1_5", "c1_5",0.02,0.045,0.52,0.425);
   c1_5->Draw();
   c1_5->cd();
   c1_5->Range(0,0,1,1);
   c1_5->SetFillStyle(4000);
   c1_5->SetBorderMode(0);
   c1_5->SetBorderSize(2);
   c1_5->SetTickx();
   c1_5->SetTicky();
   c1_5->SetLeftMargin(0.20);
   c1_5->SetRightMargin(0.05);
   c1_5->SetTopMargin(0.1);
   c1_5->SetBottomMargin(0.15);
   c1_5->SetFrameFillColor(0);
    
   frame->Draw();
   frame->GetYaxis()->SetNdivisions(4);
   frame->GetXaxis()->SetNdivisions(5);
   frame->GetYaxis()->SetLabelSize(0.0);
   frame->GetYaxis()->SetNoExponent(kFALSE);
   frame->GetYaxis()->SetTitleSize(0.0);
   frame->GetXaxis()->SetTitleSize(0.0);
   frame->GetXaxis()->SetLabelSize(0.0);

    TGaxis *axis1 = new TGaxis(-0.1,0,-0.1,20,0,20,510,"");
    axis1->SetName("axis1");
    axis1->SetLabelSize(0.1);
    axis1->SetLabelOffset(0.01);
    axis1->SetNdivisions(4);
    axis1->Draw();
    TGaxis *axis2 = new TGaxis(-0.1,0,0.1,0,-0.1,0.1,510,"");
   axis2->SetName("axis2");
   axis2->SetLabelSize(0.1);
   axis2->SetNdivisions(503);
   axis2->Draw();

   c1_5->Modified();
   c1->cd();

    combData.plotOn(frame1,RooFit::Cut("sample==sample::mode3")) ;
    simPdf.plotOn(frame1, RooFit::Slice(sample,"mode3"), RooFit::ProjWData(sample, combData));
    simPdf.plotOn(frame1, RooFit::Slice(sample,"mode3"), RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(psig3)),RooFit::LineStyle(kDashed+1), RooFit::LineColor(kRed));
    simPdf.plotOn(frame1, RooFit::Slice(sample,"mode3"), RooFit::ProjWData(sample, combData), RooFit::Components(pbkg3),RooFit::LineStyle(kDashed), RooFit::LineColor(kBlack));
    simPdf.plotOn(frame1, RooFit::Slice(sample,"mode3"), RooFit::ProjWData(sample, combData), RooFit::Components(pbkg3_k3pi),RooFit::LineStyle(kDashed), RooFit::LineColor(kCyan-3));

   frame1->GetYaxis()->SetRangeUser(0,15);

   c1_6 = new TPad("c1_6", "c1_6",0.5,0.045,1,0.425);
   c1_6->Draw();
   c1_6->cd();
   c1_6->Range(0,0,1,1);
   c1_6->SetFillStyle(4000);
   c1_6->SetBorderMode(0);
   c1_6->SetBorderSize(2);
   c1_6->SetTickx();
   c1_6->SetTicky();
   c1_6->SetLeftMargin(0.20);
   c1_6->SetRightMargin(0.05);
   c1_6->SetTopMargin(0.1);
   c1_6->SetBottomMargin(0.15);
   c1_6->SetFrameFillColor(0);

   frame1->Draw();
   frame1->GetYaxis()->SetNdivisions(504);
   frame1->GetXaxis()->SetNdivisions(5);
   frame1->GetYaxis()->SetLabelSize(0.0);
   frame1->GetYaxis()->SetNoExponent(kFALSE);
   frame1->GetYaxis()->SetTitleSize(0.0);
   frame1->GetXaxis()->SetTitleSize(0.0);
   frame1->GetXaxis()->SetLabelSize(0.0);

      TGaxis *axis1 = new TGaxis(0.9,0,0.9,14,0,14,510,"");
    axis1->SetName("axis1");
    axis1->SetLabelSize(0.1);
    axis1->SetLabelOffset(0.01);
    axis1->SetNdivisions(4);
    axis1->Draw();
   TGaxis *axis2 = new TGaxis(0.9,0,1.68,0,0.9,1.68,510,"");
   axis2->SetName("axis2");
   axis2->SetLabelSize(0.1);
   axis2->SetNdivisions(503);
   axis2->Draw();


    TLatex l3;
    l3.SetNDC();
    l3.SetTextSize(0.035);
    l3.SetTextSize(0.06);
    l3.DrawLatex(0.55, 0.83,  Form("nsig = %.1f #pm %.1f",Nsig3,Nsig3_err));
    l3.DrawLatex(0.55, 0.78,  Form("nbkg = %.1f #pm %.1f",nbkg3.getVal()+nbkg3_k3pi.getVal(),n3_err));

    c1_6->Modified();
    c1->cd();

    c1->Modified();
    c1->cd();
    c1->SetSelected(c1);
///////////////////////////////////////////////////////////////////modes///////////////////
   
    TCanvas *cv2 = new TCanvas("cv2", "", 1200,800);
    cv2->Divide(2,1); 
    RooPlot *framex = Umissfit.frame(40);
    RooPlot *framey = mhad.frame(26);
     
    cv2->cd(1);
    combData.plotOn(framex);
    simPdf.plotOn(framex, RooFit::ProjWData(sample, combData));
    simPdf.plotOn(framex, RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(psig0,psig1,psig3)),RooFit::LineStyle(kDashed+1), RooFit::LineColor(kRed));
    simPdf.plotOn(framex, RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(pbkg0,pbkg1,pbkg3)),RooFit::LineStyle(kDashed), RooFit::LineColor(kBlack));
    simPdf.plotOn(framex, RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(pbkg0_k3pi,pbkg1_k3pi,pbkg3_k3pi)),RooFit::LineStyle(kDashed), RooFit::LineColor(kCyan-3));
    framex->SetXTitle("U_{missfit} (GeV^{2})");
    framex->SetTitle("simultaneous fit");
    framex->Draw();

    cv2->cd(2);
    combData.plotOn(framey);
    simPdf.plotOn(framey, RooFit::ProjWData(sample, combData));
    simPdf.plotOn(framey, RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(psig0,psig1,psig3)),RooFit::LineStyle(kDashed+1), RooFit::LineColor(kRed));
    simPdf.plotOn(framey, RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(pbkg0,pbkg1,pbkg3)),RooFit::LineStyle(2), RooFit::LineColor(kBlack));
    simPdf.plotOn(framey, RooFit::ProjWData(sample, combData), RooFit::Components(RooArgSet(pbkg0_k3pi,pbkg1_k3pi,pbkg3_k3pi)),RooFit::LineStyle(2), RooFit::LineColor(kCyan-3));
    framey->SetXTitle("M_{K#pi#pi} (GeV^{2})");
    framey->SetTitle("simultaneous fit");
    framey->Draw();
    
    RooRealVar nsig_sp_m0("nsig_sp_m0", "", nsig0.getVal());
    RooAddPdf allpdf_sp_m0("allpdf_sp_m0", "", RooArgList(psig0,pbkg0,pbkg0_k3pi), RooArgList(nsig_sp_m0,nbkg0,nbkg0_k3pi));
    RooRealVar nsig_sp_m1("nsig_sp_m1", "", nsig1.getVal());
    RooAddPdf allpdf_sp_m1("allpdf_sp_m1", "", RooArgList(psig1,pbkg1,pbkg1_k3pi), RooArgList(nsig_sp_m1,nbkg1,nbkg1_k3pi));
    RooRealVar nsig_sp_m3("nsig_sp_m3", "", nsig3.getVal());
    RooAddPdf allpdf_sp_m3("allpdf_sp_m3", "", RooArgList(psig3,pbkg3,pbkg3_k3pi), RooArgList(nsig_sp_m3,nbkg3,nbkg3_k3pi));
    //mode0
    RooStats::SPlot* splot = new RooStats::SPlot("splot","splot",*ds0,&allpdf_sp_m0,RooArgList(nsig_sp_m0,nbkg0,nbkg0_k3pi));
    splot->Print();
    RooDataSet * dataw_s = new RooDataSet(ds0->GetName(),ds0->GetTitle(),ds0,*ds0->get(),0,Form("%s_sw", nsig_sp_m0.GetName())) ;
    TH1* hmhad0 = dataw_s->createHistogram(Form("hmhad_%s",nsig_sp_m0.GetName()), mhad, RooFit::Cut(""));// ?
    delete splot; 
    //mode1
    splot = new RooStats::SPlot("splot","splot",*ds1,&allpdf_sp_m1,RooArgList(nsig_sp_m1,nbkg1,nbkg1_k3pi));
    splot->Print();
    delete dataw_s; 
    dataw_s = new RooDataSet(ds1->GetName(),ds1->GetTitle(),ds1,*ds1->get(),0,Form("%s_sw", nsig_sp_m1.GetName())) ;
    TH1* hmhad1 = dataw_s->createHistogram(Form("hmhad_%s",nsig_sp_m1.GetName()), mhad, RooFit::Cut(""));
    delete splot;
    //mode3
    splot = new RooStats::SPlot("splot","splot",*ds3,&allpdf_sp_m3,RooArgList(nsig_sp_m3,nbkg3,nbkg3_k3pi));
    splot->Print();
    delete dataw_s; 
    dataw_s = new RooDataSet(ds3->GetName(),ds3->GetTitle(),ds3,*ds3->get(),0,Form("%s_sw", nsig_sp_m3.GetName())) ;
    TH1* hmhad3 = dataw_s->createHistogram(Form("hmhad_%s",nsig_sp_m3.GetName()), mhad, RooFit::Cut(""));

    TH1* hmhadall = (TH1*)hmhad0->Clone("hmhadall");
    hmhadall->Add(hmhad1);
    hmhadall->Add(hmhad3);
    hmhadall->SetMinimum(0);
    hmhadall->SetStats(0);
    hmhadall->SetTitle("");
    hmhadall->SetXTitle("M(K#pi#pi) [GeV]");

     c1->Print("figure/data2D_modes.eps");
    cv2->Print("figure/data2D.eps");
    ofstream out;
    out.open("figure/data2D.dat");
    out<<"bf : "<<bf.getVal()<<"\t"<<bf.getError()<<"\n";
    out<<"mode0 yields : "<<Nsig0<<"\t"<<Nsig0_err<<"\n";
    out<<"mode1 yields : "<<Nsig1<<"\t"<<Nsig1_err<<"\n";
    out<<"mode3 yields : "<<Nsig3<<"\t"<<Nsig3_err<<"\n";
    out<<"mode0 bkg    : "<<nbkg0.getVal()+nbkg0_k3pi.getVal()<<"\t"<<n0_err<<"\n";
    out<<"mode1 bkg    : "<<nbkg1.getVal()+nbkg1_k3pi.getVal()<<"\t"<<n1_err<<"\n";
    out<<"mode3 bkg    : "<<nbkg3.getVal()+nbkg3_k3pi.getVal()<<"\t"<<n3_err<<"\n";
    out<<"total fit    : "<<Nsig0+Nsig1+Nsig3<<"+/-"<<Nsig1_err*(Nsig0+Nsig1+Nsig3)/Nsig1;
 
     double subdecay = 0.32852;
     double bf_1  = bf.getVal(),   bf_1_err  = bf.getError();
     double bf_K1 = bf_1/subdecay, bf_K1_err = bf_1_err/subdecay; 
     ofstream out_1("figure/bf.txt");
     out_1<<setprecision(9)<<bf_K1<<"\n";
     out_1<<setprecision(9)<<bf_K1_err<<"\n";

     //bf.setVal(0); bf.setConstant();
     //RooFitResult* res2 = simPdf.fitTo(combData, RooFit::Save(), RooFit::Minos(1));
     //minNll += res2->minNll();
     //cout<<"Statisical significance: "<<sqrt(2*minNll)<<endl;
     hmhadall->Draw();
    

}
