#include<cstdio>
#include<string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <iomanip>
        
 #include "RooAbsReal.h"
 #include "RooFit.h"
 #include "RooRealVar.h"
 #include "TMath.h"
using namespace std;

void mc_1_mode3_deltaE(){
    gSystem->Load("libRooFit");
    using namespace RooFit;    

    gStyle->SetLabelSize(0.07,"xyz");
    gStyle->SetNdivisions(507,"xyz");
    gROOT->SetStyle("Plain");
    gStyle->SetLabelSize(0.06,"xyz");
    gStyle->SetNdivisions(405,"xyz");

    gStyle->SetPadTopMargin(.10);
    gStyle->SetPadLeftMargin(.20);
    gStyle->SetPadRightMargin(.05);
    gStyle->SetPadBottomMargin(.15);
    gStyle->SetTitleSize(0.06,"xyz");
    gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(0.5);
    gStyle->SetTitle("");
    gStyle->SetOptStat(1);
    gStyle->SetEndErrorSize(0);
    using namespace RooFit;

for (int nmode=3; nmode<4; nmode++) {
    if (nmode == 2) continue;   
    RooRealVar mBC("mBC","", 1.8365,1.8865); 
    TFile *f = TFile::Open("root/et1_dh3_mc.root");
    TH1F *h1 = (TH1F*) f->Get("h0");
    RooDataHist* hdata = new RooDataHist("hdata","",mBC,h1);
  
    if(nmode==0)     TChain sig0("tag0");
    if(nmode==0)     sig0.Add("/besfs/groups/psipp/psippgroup/public/qusq/Phix/script/0_shape_mbc/shape_tag0.root");
    if(nmode==1)     TChain sig1("tag1");
    if(nmode==1)     sig1.Add("/besfs/groups/psipp/psippgroup/public/qusq/Phix/script/0_shape_mbc/shape_tag1.root"); 
    if(nmode==3)     TChain sig3("tag3");
    if(nmode==3)     sig3.Add("/besfs/groups/psipp/psippgroup/public/qusq/Phix/script/0_shape_mbc/shape_tag3.root");
    TH1F* h2 = new TH1F("h2","",200,1.8365,1.8865);
    if(nmode==0) sig0.Draw("mbc>>h2"); 
    if(nmode==1) sig1.Draw("mbc>>h2"); 
    if(nmode==3) sig3.Draw("mbc>>h2"); 
    RooDataHist* hshape = new RooDataHist("hshape","",mBC,h2); 
    RooHistPdf sigshape("sigshape","",mBC,*hshape,2);
   if(nmode==0)    RooRealVar  mean1( "mean1" , "", -3e-3,3e-3);
   if(nmode==0)    RooRealVar  mean2( "mean2" , "", -1e-3,1e-3);
   if(nmode==0)    RooRealVar  sigma1("sigma1", "", 0.,4e-3);
   if(nmode==0)    RooRealVar  sigma2("sigma2", "", 0.,3e-3);

   if(nmode==1)    RooRealVar  mean1( "mean1" , "", -1e-3,1e-3);  
   if(nmode==1)    RooRealVar  mean2( "mean2" , "", -1e-3,1e-3);
   if(nmode==1)    RooRealVar  sigma1("sigma1", "", 0.,1e-3);
   if(nmode==1)    RooRealVar  sigma2("sigma2", "", 0.,2e-3);
   
   if(nmode==3)    RooRealVar  mean1( "mean1" , "", -2e-3,2e-3);
   if(nmode==3)    RooRealVar  mean2( "mean2" , "", -2e-3,2e-3);
   if(nmode==3)    RooRealVar  sigma1("sigma1", "", 0.,1e-2);
   if(nmode==3)    RooRealVar  sigma2("sigma2", "", 0.,3e-3);
                                                             
    RooGaussian gauss1("gauss1", "", mBC, mean1, sigma1 );
    RooGaussian gauss2("gauss2", "", mBC, mean2, sigma2 );
    RooRealVar  frac1("frac1","", 0.0, 1.0) ;
    RooAddModel gauss("gauss","",RooArgList(gauss1,  gauss2),frac1) ;
    RooFFTConvPdf sigpdf("sigpdf", "", mBC,sigshape, gauss);
    
    RooRealVar endpt("endpt","",1.8865);
    RooRealVar kappa("kappa","", -50., 30.);
    RooRealVar arg2  ("arg2", "", 0.3, 10.); 
    RooArgusBG bkg("bkg","",mBC,endpt,kappa,arg2);
    
    double sum_data = hdata.sumEntries();   
    RooRealVar Nsig("Nsig","", 0.9*sum_data, 0.0, sum_data);
    RooRealVar Nbkg("Nbkg","", 0.1*sum_data, 0.0, sum_data);

    RooAddPdf  model("model","", RooArgList(sigpdf,bkg), RooArgList(Nsig,Nbkg));
    /**************************optimize, written by zhangsf*******************************/
    RooAbsReal *nll = model.createNLL(*hdata);
    RooMinimizer m(*nll);
    RooFitResult *rfr;

    Int_t status = 1;
    Int_t covQual = 0;
    Double_t edm = 100.;
    Double_t chi2 = 999.;
    RooFitResult *rfr = model.fitTo(*hdata,Save(),Extended(kTRUE),Warnings(kFALSE),PrintEvalErrors(-1));
    status = rfr->status();
    edm = rfr->edm();
    covQual = rfr->covQual();
    Int_t times = 0;

    int seedn;
    ifstream seedx;
    seedx.open("../seed/seed1.dat");
    for (int i=0; i<1; i++) { seedx >> seedn; }
       int seed1=seedn*11;
       int seed2=seedn*22;
       int seed3=seedn*33;
       int seed4=seedn*44;
       int seed5=seedn*55;
       int seed6=seedn*66;
       int seed7=seedn*77;
       TRandom3 r1(seed1);
       TRandom3 r2(seed2);
       TRandom3 r3(seed3);
       TRandom3 r4(seed4);
       TRandom3 r5(seed5);
       TRandom3 r6(seed6);
       TRandom3 r7(seed7);
/*    
    RooPlot* rp1 = mBC.frame();
    hdata->plotOn(rp1, LineWidth(3),Name("h_data"));
    model.plotOn (rp1, LineWidth(2),Name("curve"));
*/
    while ((edm>1||status!=0||covQual!=3/*||chi2>3*/)&&times<1000){
     Double_t m_mean1  = r1.Uniform(-3e-3,3e-3); 
     Double_t m_mean2  = r1.Uniform(-3e-3,3e-3);
     Double_t m_sigma1 = r1.Uniform(0,1e-2);
     Double_t m_sigma2 = r1.Uniform(0,3e-3);
     Double_t m_frac1  = r1.Uniform(0,1);     
     Double_t m_kappa  = r1.Uniform(-50.,30.);
     Double_t m_arg2   = r1.Uniform(0.3, 10.);
     mean1.setVal(m_mean1);
     mean2.setVal(m_mean2);
     sigma1.setVal(m_sigma1);
     sigma2.setVal(m_sigma2);
     frac1.setVal(m_frac1);
     kappa.setVal(m_kappa);
     arg2.setVal(m_arg2);
     
     m.migrad();
     rfr = m.save("fitresult","fitresult");
     status = rfr->status();
     if (status != 0) {
              cout << "Wrong status from Migard." << endl;
              continue;
      }
      m.hesse();
      rfr = m.save("fitresult","fitresult");
 

     status = rfr->status();
     edm = rfr->edm(); //?
     covQual = rfr->covQual();//?
     //model.plotOn(rp1, LineWidth(2),Name("curve"));
     //chi2 = rp1->chiSquare(9);
     rfr->printValue(cout);
     times++;
     }   
     /*********************************end of the optimize*************************************/

     //RooFitResult* rfr = model.fitTo(*hdata, RooFit::Save(), Extended(kTRUE));
     //rfr->Print("v");
      
     RooPlot* rp1 = mBC.frame();
     rp1->GetXaxis()->CenterTitle( kTRUE );
     rp1->GetXaxis()->SetLabelSize(0.05);
     rp1->GetYaxis()->SetLabelSize(0.045);
     rp1->GetXaxis()->SetTitleSize(0.04);
     rp1->GetYaxis()->SetTitleSize(0.04);
     rp1->GetXaxis()->SetTitleOffset(0.9);
     rp1->GetYaxis()->SetTitleOffset(1.5);
     rp1->GetYaxis()->SetNdivisions(505);
     rp1->SetTitleOffset( 1.35, "x" );
     rp1->SetXTitle( "M_{BC} (GeV/c^{2})" );
     rp1->GetYaxis()->CenterTitle( kTRUE );
     rp1->GetYaxis()->SetNdivisions(505);
     rp1->GetXaxis()->SetNdivisions(505);
     double x = rp1->GetMaximum();
     rp1->SetAxisRange(0,1.1*x, "Y");
     rp1->SetTitleOffset(1.5, "y");   
     rp1->SetYTitle( "Event (/0.5 MeV/#font[12]{c}^{2})" );
 
     hdata->plotOn(rp1,MarkerStyle(20), MarkerSize(1.0),MarkerColor(kBlack) , Name("h_data"));
     model.plotOn(rp1, LineWidth(2), Name("curve"));
     double chi2 = rp1->chiSquare();
     model.plotOn(rp1, Components(sigpdf), LineColor(4),LineWidth(1),LineStyle(2));
     //double chi2 = rp1->chiSquare();
     model.plotOn(rp1, Components(bkg), LineColor(2),LineWidth(1),LineStyle(2));
    /*************************************************calculate the chisq****************************************************/
    
    RooHist* h_data = rp1->getHist("h_data");
    RooCurve* curve = rp1->getCurve("curve");
    Double_t xstart, xstop, y;
    curve->GetPoint(0, xstart, y);             //first bins ~ get x and y value of point number i(or say, the ith bin, i = 0,1,2...)
    curve->GetPoint(curve->GetN()-1, xstop, y);//last bins ~ GetN(), get the number of bins(pionts) 
    RooHist *hist1 = new RooHist(h_data->getNominalBinWidth());// (average bin width)
    RooHist *hist2 = new RooHist(h_data->getNominalBinWidth());
    for (Int_t i=0; i<h_data->GetN(); i++) { 
      Double_t x, point;
      h_data->GetPoint(i, x, point);          //piont, the y value in ith bins
      if (x<xstart||x>xstop) continue;
      if (point==0) continue;
      Double_t yy = point - curve->interpolate(x);//linearly interpolated value of curve at x
      Double_t dl = h_data->GetErrorYlow(i);      //TGraphErrors class
      Double_t dh = h_data->GetErrorYhigh(i);
      hist2->addBinWithError(x,point,dl,dh);      
      Double_t norm = (yy>0?dl:dh);              //interesting!
      if (norm==0) continue;    // ???           
      yy /= norm;               //yy = (y_data - y_pdf) / data_err
      dl /= norm;               //dl = (data_low_err) / data_err   ??
      dh /= norm;               //dh = (data_high_err) / data_err  ??
      hist1->addBinWithError(x,yy,dl,dh);
    }
   
   /********************************************end of calculate the chisq***********************************************/
   
  //plot on 
    RooPlot* pull = mBC.frame(Title("Pull distribution")) ;
    pull->addPlotable(hist1,"P") ;
    pull->GetYaxis()->SetRangeUser(-5,5);
    pull->GetYaxis()->SetLabelSize(0.10);
    pull->GetYaxis()->SetTickLength(0);
    pull->Draw();
   
    TCanvas* cvs = new TCanvas("cvs","c", 800, 600 );
    cvs->cd(1);
    TPad *pad1 = new TPad("pad1","pad1",0.01,0.20,0.99,0.99);
    TPad *pad2 = new TPad("pad2","pad2",0.01,0.01,0.99,0.23);
    pad1->Draw(); 
    pad2->Draw();
     
    pad1->cd(); 
    pad1->SetLeftMargin(0.15);
    pad1->SetRightMargin(0.05);
        
    rp1->Draw();
    rp1->SetMarkerSize(0.5);

     TLatex lt1;
     lt1.SetNDC();
     lt1.SetTextAngle(0);
     lt1.SetTextSize(0.06);
     
     lt1.DrawLatex(0.22, 0.81,  "#bar{D}^{0}#rightarrowK^{+}#pi^{-}#pi^{+}#pi^{-}");//mode3
     //lt1.DrawLatex(0.22, 0.76,  Form("#chi^{2}_{ndof} = %.1d", rp1->chiSquare()));
     lt1.DrawLatex(0.22, 0.76,  "#chi^{2}_{ndof} =2.7");

     TArrow *ar1 = new TArrow(1.858, 0.5*h1->GetMaximum(),1.858,0.1*h1->GetMaximum() ,0.01,"|>");
     TArrow *ar2 = new TArrow(1.874, 0.5*h1->GetMaximum(),1.874,0.1*h1->GetMaximum() ,0.01,"|>");
 
     ar1->SetLineColor(kBlue);
     ar1->SetFillColor(kBlue);
     ar1->SetLineWidth(2);
     ar1->Draw();
     ar2->SetLineColor(kBlue);
     ar2->SetFillColor(kBlue);
     ar2->SetLineWidth(2);
     ar2->Draw();

    hist2->Draw("PEsame");
    pad2->cd();
    pad2->SetLeftMargin(0.15);
    pad2->SetRightMargin(0.05);
    pull->Draw();
    

    /********************************************number yield**************************************************/ 
    mBC.setRange("signal",1.858,1.874);
    RooAbsReal* intall_signal =  sigpdf.createIntegral(mBC, NormSet(mBC),Range("signal"));
    double N_ratio     =  intall_signal->getVal();    //ratio
    double N_fit       =  (Nsig.getVal())*N_ratio;
    double N_fit_err =  (Nsig.getError())*N_ratio; 
    double N_bkg      = (Nbkg.getVal())*N_ratio;
    
    double N_prod3 = 2055000.0;//200,0000+5,5000
    double N_dcs3  = 2*N_prod3*2.65e-4*0.382, N_107 = 2*N_prod3*2.17e-3*0.692*0.0319;
    if (nmode==3) double eff = (N_fit-N_dcs3-N_107)/N_prod3;

    /************************************************end********************************************************/
    ofstream fit_data;
    fit_data.open(Form("./figure/mode3/mc_1_%d_mc.dat",nmode), ios::app);
    fit_data<<"fit yield: "<<N_fit<<"\t"<<N_fit_err<<"\n";
    fit_data<<"bkg: "<< N_bkg<<"\n";
    fit_data<<"chisq: "<<chi2<<"\n";
    fit_data<<"eff: "<<eff<<" err:"<<N_fit_err/N_prod3<<"\n";
    fit_data<<"--------------------------------"<<"\n";
    fit_data<<"-----------parameter------------"<<"\n";
    fit_data<<"mean1: "<<mean1.getVal()<<"  sigma1: "<<sigma1.getVal()<<"\n";
    fit_data<<"mean2: "<<mean2.getVal()<<"  sigma2: "<<sigma2.getVal()<<"\n";
    fit_data<<"frac:  "<<frac1.getVal()<<"\n";
    fit_data<<"kappa: "<<kappa.getVal()<<"  arg2: "<<arg2.getVal()<<"\n";

    


    //std::ofstream rs("./dat/mBC_1D/mode0.dat", ios::app);
    //rs<<std::setiosflags(ios::left)<<"fit among signal region: "<<setw(10)<<std::setprecision(8)<<N_fit<<"+/-"<<N_fit_err<<" ; "<<endl;
    	    
    /*TLatex lt;
      lt.SetNDC();
      lt.SetTextAngle(0);
      lt.SetTextSize(0.05);
    //pad1->cd();
      lt.DrawLatex(0.68, 0.83,  Form("%.0f #pm %.0f",N_fit,N_fit_err));
    //lt.DrawLatex(0.70, 0.78,  Form("#chi^{2} = %.2f",chi2));    
    */

     cvs->Modified();
     cvs->cd();
     cvs->SetSelected(cvs);

     cvs->Print(Form("figure/mode3/mc_1_%d_mc.eps",nmode),"recreate");    
     cout<<"chi2 "<<chi2<<endl;
     cout<<"eff "<<eff<<endl;
     }
}
