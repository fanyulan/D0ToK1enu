void u_m(){
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
bes3Style->SetPadTopMargin(.01);
bes3Style->SetPadLeftMargin(.15);
bes3Style->SetPadRightMargin(.05);
bes3Style->SetPadBottomMargin(.15);
Int_t font=132;      //times new roman, reg
Double_t tsize=0.05;
bes3Style->SetTextFont(font);
bes3Style->SetTextSize(tsize);

bes3Style->SetLabelSize(tsize,"xyz");
bes3Style->SetLabelOffset(0.01,"xyz");

bes3Style->SetTitleFont(font,"xyz");
bes3Style->SetTitleSize(tsize,"xyz");
bes3Style->SetTitleXOffset(1.);
bes3Style->SetTitleYOffset(1.4); //offset the title of y axis a bit
bes3Style->SetTitleBorderSize(2.);
bes3Style->SetMarkerStyle(0);
bes3Style->SetMarkerSize(0.8);
bes3Style->SetFrameBorderMode (0.);
bes3Style->SetFrameLineWidth  (2.);
bes3Style->SetLineWidth(1.0);
bes3Style->SetHistLineWidth(2.);
bes3Style->SetLineStyleString(2,"[12 12]");
bes3Style->SetErrorX(0.001);
bes3Style->SetOptTitle(0);     //no title box
bes3Style->SetOptStat(0);    //no stat info
bes3Style->SetLineStyleString(2,"[30 10]");
bes3Style->SetLineStyleString(3,"[4 8]");
bes3Style->SetLineStyleString(4,"[15 12 4 12]");
bes3Style->SetLineStyleString(5,"[15 15]");
bes3Style->SetLineStyleString(6,"[15 12 4 12 4 12]");
bes3Style->SetOptDate(0);
bes3Style->SetDateY(.98);
bes3Style->SetStripDecimals(kFALSE);

bes3Style->SetEndErrorSize(0.0); //make the end of error bar longer 

gROOT->SetStyle("bes3");
gROOT->ForceStyle();


	gSystem->Load("libRooFit");
	using namespace RooFit;
    gStyle->SetLabelSize(0.07,"xyz");
	gStyle->SetNdivisions(507,"xyz");
           

    gROOT->SetStyle("Plain");
    gStyle->SetLabelSize(0.06,"xyz");
    gStyle->SetNdivisions(405,"xyz");
    gStyle->SetPadTopMargin(.10);
    gStyle->SetPadLeftMargin(.15);
    gStyle->SetPadRightMargin(.05);
    gStyle->SetPadBottomMargin(.15);
    gStyle->SetTitleSize(0.06,"xyz");
    gStyle->SetOptTitle(0);
    gStyle->SetMarkerSize(0.5);
    gStyle->SetTitle("");
    gStyle->SetOptStat(0);
    gStyle->SetEndErrorSize(0);

const double xoff = 1.;
const double yoff = 1.7;

   TChain ch("tagD");
   ch.Add("/scratchfs/bes/fangyl/data1/data/all/*.root");
   TCut cuta = "isbest==1&&type==1&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)&&m_pipi>0.31&&abs(m_pipi-0.497611)>0.01&&abs(deltaE_pipzswp)>0.012";
   TCut cutb = "abs(Umissfit)<0.1&&mhad<1.68&&mhad>0.9";
   TCut cutc = "elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge";
   TCut cut0 = "mode==0&&clveto==1&&deltaE>-0.029&&deltaE<0.027";
   TCut cut1 = "mode==1&&deltaE>-0.069&&deltaE<0.038";
   TCut cut3 = "mode==3&&deltaE>-0.031&&deltaE<0.028";

   double bins = 100;
   TH2F* h = new TH2F("h","",bins,-0.1,0.1,bins,0.9,1.68);//m_miss2(x), mKpipi(y)
   ch.Draw("mhad:Umissfit>>h", cuta&&cutb&&cutc&&(cut0||cut1||cut3));  // Atention  (y,x)
   

   TCanvas *c=new TCanvas("c","",0,0,600,600);
   double nh_u=h->GetEntries();
   TPad *p1 =new TPad("c","",0.,0.0,1,1);
   p1->Draw("(Umissfit:mhad>>h");
   p1->cd();
   p1->SetFillStyle(4000);
   p1->Range(0,0,1,1);
   p1->SetLeftMargin(0.15);
   p1->SetRightMargin(0.05);
   p1->SetTopMargin(0.05);
   p1->SetBottomMargin(0.13);
   p1->SetFrameFillColor(1);
  
   //h->Draw("COL");//COLZ
   h->Draw("p");//COLZ
   h->SetMarkerStyle(20);
   h->GetXaxis()->CenterTitle( kTRUE );
   h->SetTitleOffset( xoff, "x" );
   h->SetXTitle("M_{miss}^{2}  /(2 MeV^{2}/#font[12]{c}^{4})");
   h->GetYaxis()->CenterTitle( kTRUE );
   h->GetXaxis()->SetLabelSize(0.05);
   h->GetYaxis()->SetLabelSize(0.05);
   h->GetXaxis()->SetTitleSize(0.0);
   h->GetYaxis()->SetTitleSize(0.0);
   h->GetYaxis()->SetNdivisions(505);
   h->GetXaxis()->SetNdivisions(505);
   h->SetTitleOffset( yoff, "y");
   h->SetYTitle("M_{K^{-}#pi^{+}#pi^{-}} /(7.8 MeV/#font[12]{c}^{2})");
   h->SetStats(false);

    leg = new TLegend(0.52,0.01,0.64,0.05);
    leg->SetTextAlign(22);
    leg->SetTextAlign(22);
    leg->SetTextFont(22);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05);
    leg->SetTextColor(1);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetLineWidth(0);
    leg->SetTextAngle(0);
    leg->SetHeader("M_{miss}^{2} /(2 MeV^{2}/#font[12]{c}^{4})");
    leg->Draw();
    leg = new TLegend(0.015,0.45,0.055,0.60);
    leg->SetTextAlign(22);
    leg->SetTextAlign(22);
    leg->SetTextFont(22);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.05);
    leg->SetTextColor(1);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetLineWidth(0);
    leg->SetTextAngle(90);
    leg->SetHeader("M_{K^{-}#pi^{+}#pi^{-}} /(7.8 MeV/#font[12]{c}^{2})");
    leg->Draw();
    leg1=new TLatex();
    leg1->SetTextFont(22);
    leg1->SetTextColor(1);
    leg1->SetTextSize(0.05);
    leg1->SetTextAlign(1);

    TLatex lt;
    lt.SetNDC();
    lt.SetTextAngle(0);
    lt.SetTextSize(0.06);
    //lt.DrawLatex(0.5, 0.92,"signal MC");
    //lt.DrawLatex(0.75, 0.85,"(a)");

    c->cd();
    c->Print("./Dalitz.eps");
//    c->Print("./Dalitz_inclusiveMC.eps");


  
}



