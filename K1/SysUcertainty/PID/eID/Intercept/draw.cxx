#include <iostream.h>
#include <iomanip.h>
using namespace RooFit;
using namespace std;
void draw()
{
//draw
gROOT->SetStyle("Plain");
gStyle->SetMarkerSize(0.9);
gStyle->SetTitleOffset(1.2,"x");
gStyle->SetTitleOffset(1.6,"yz");
gStyle->SetTitleSize(0.05,"xyz");
gStyle->SetLabelSize(0.05,"xyz");
gStyle->SetNdivisions(505,"xyz");
gStyle->SetPadTopMargin(.10);
gStyle->SetPadLeftMargin(.20);
gStyle->SetPadRightMargin(.10);
gStyle->SetPadBottomMargin(.15);
gStyle->SetPalette(1);
gStyle->SetOptTitle(1);
gStyle->SetFillColor(4000);

 Double_t m_cos[10]  = {-0.865,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.865};
 Double_t m_ecos[10] = {0.065,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.065};
 Double_t m_p[5] = {0.15,0.25,0.35,0.45,0.55};
 Double_t m_ep[5]= {0.05,0.05,0.05,0.05,0.05};
 //Double_t m_p[19]     = {0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.25,1.35,1.45,1.55,1.65,1.75,1.85,1.95};
 //Double_t m_ep[19]    = {0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05}; 
 
 ifstream openfile;
 //openfile.open("./dat/minus_delta_eff_1.txt");
 openfile.open("./dat/plus_delta_eff.txt");
   
 TGraph2D *dt1 = new TGraph2D();
 TGraph2D *dt2 = new TGraph2D();
 TGraph2D *dt3 = new TGraph2D();
 TGraph2D *dt4 = new TGraph2D();

 Int_t k = 0;
 for (Int_t i = 0; i<5; i++) {
  for (Int_t j=0; j<10; j++) {
    double eff_data, e_eff_data, eff_mc, e_eff_mc, delta, diff, e_delta;
    openfile>>eff_data;
    openfile>>e_eff_data;
    openfile>>eff_mc;
    openfile>>e_eff_mc;
    openfile>>delta;
    openfile>>e_delta;
    openfile>>diff;

    dt1->SetPoint(k,m_p[i],m_cos[j],eff_data);
    dt2->SetPoint(k,m_p[i],m_cos[j],eff_mc);
    dt3->SetPoint(k,m_p[i],m_cos[j],delta);
    dt4->SetPoint(k,m_p[i],m_cos[j],diff);
    k++;
}}
//draw
gStyle->SetMarkerSize(0.9);
gROOT->SetStyle("Plain");
gStyle->SetLabelSize(0.06,"xyz");
gStyle->SetNdivisions(405,"xyz");
gStyle->SetPadTopMargin(.05);
gStyle->SetPadLeftMargin(.20);
gStyle->SetPadRightMargin(.05);
gStyle->SetPadBottomMargin(.20);
gStyle->SetTitleSize(0.06,"xyz");
gStyle->SetOptTitle(0);
gStyle->SetMarkerSize(0.5);
gStyle->SetTitle("");
gStyle->SetOptStat(0);
gStyle->SetEndErrorSize(0);
using namespace RooFit;

   TCanvas *resf1= new TCanvas("resf1"," ",500,500);
   resf1->Divide(1,1);
   resf1->cd(1);
   resf1->cd(1)->SetGrid();
   dt1->SetMarkerColor(2);
   dt1->SetMarkerStyle(22);
   dt1->SetTitle(";P (GeV/c);cos#theta;#epsilon_{data}");
   dt1->Draw("surf1");
   //resf1->Print("./figure/minus_eff_data_1.eps");
   resf1->Print("./figure/plus_eff_data1.eps");

   TCanvas *resf2= new TCanvas("resf2"," ",500,500);
   resf2->Divide(1,1);
   resf2->cd(1);
   resf2->cd(1)->SetGrid();
   dt2->SetMarkerColor(2);
   dt2->SetMarkerStyle(22);
   dt2->SetTitle(";P (GeV/c);cos#theta;#epsilon_{MC}");
   dt2->Draw("surf1");
   //resf2->Print("./figure/minus_eff_mc_1.eps");
   resf2->Print("./figure/plus_eff_mc1.eps");

   TCanvas *resf3= new TCanvas("resf3"," ",500,500);
   resf3->SetTitle("");
   resf3->Divide(1,1);
   resf3->cd(1);
   resf3->cd(1)->SetGrid();
   dt3->SetMarkerColor(2);
   dt3->SetMarkerStyle(22);
   dt3->SetTitle(";P (GeV/c);cos#theta;#Delta=#epsilon_{data}/#epsilon_{MC}");
   dt3->Draw("surf1");
   //resf3->Print("./figure/minus_delta_1.eps");
   resf3->Print("./figure/plus_delta1.eps");

   TCanvas *resf4= new TCanvas("resf4"," ",500,500);
   resf4->SetTitle("");
   resf4->Divide(1,1);
   resf4->cd(1);
   resf4->cd(1)->SetGrid();
   dt4->SetMarkerColor(2);
   dt4->SetMarkerStyle(22);
   dt4->SetTitle(";P (GeV/c);cos#theta;#Delta=#epsilon_{data}/#epsilon_{MC}-1");
   dt4->Draw("surf1");
   //resf4->Print("./figure/minus_diff_1.eps");
   resf4->Print("./figure/plus_diff1.eps");
}
