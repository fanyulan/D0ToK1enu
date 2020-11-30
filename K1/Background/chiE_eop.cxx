void chiE_eop(){
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

const double xoff = 1.2;
const double yoff = 1.7;

   TChain ch("tagD");
   ch.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/sigMC/update/BELLE/fit2/red_truth/truth.root");
 
   TCut cut_sig = "isTrueTag&&isSignal&&abs(KpTrueID)==321&&abs(pi1TrueID)==211&&abs(pi2TrueID)==211&&abs(elecTrueID)==11&&abs(elecMomID)==421";
   TCut cuta    = "isbest==1&&type==1&&mBC<1.874&&mBC>1.858&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&m_pipi>0.31&&abs(m_pipi-0.497611)>0.01&&abs(deltaE_pipzswp)>0.013";
   TCut cutb    = "abs(Umissfit)<0.2";
   TCut cutc = "elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge";
   TCut cut0 = "mode==0&&clveto==1&&deltaE>-0.029&&deltaE<0.027";
   TCut cut1 = "mode==1&&deltaE>-0.069&&deltaE<0.038";
   TCut cut3 = "mode==3&&deltaE>-0.031&&deltaE<0.028";
   TCut cut_K1 = "((mcmodea==666||mcmodea==222)&&mcmode1==103)||((mcmodeb==-666||mcmodeb==-222)&&mcmode2==-103)";
   TCut cut_K2 = "mcmode1==104 || mcmode2==-104";

   double bins = 100;
   TH2F* h = new TH2F("h","",bins,0.0,4.0,bins,0.0,1.5);//dedx_chiE(x), eop(y)
   ch.Draw("eop_E:dedxchi_E>>h", /*cut_sig&&cut_K1&&*/!cut_K2&&cuta&&cutb&&cutc&&(cut0||cut1||cut3)  );  // Atention  (y,x)
   
   TF1 *f_relevence = new TF1("f_relevence","relevence(x,y)", 0.0, 3.5);
   f_relevence->SetLineColor(8);
   TF1 *f_relevence_factor = new TF1("f_relevence_factor","relevence_factor(x,y)", 0.0, 3.5);
   f_relevence_factor->SetLineColor(2);

   TCanvas *c=new TCanvas("c","",0,0,800,600);
   double nh_u=h->GetEntries();
   TPad *p1 =new TPad("c","",0.,0.0,1,1);
   p1->Draw("(eop_E:dedxchi_E>>h");
   p1->cd();
   p1->SetFillStyle(4000);
   p1->Range(0,0,1,1);
   p1->SetLeftMargin(0.15);
   p1->SetRightMargin(0.1);
   p1->SetTopMargin(0.1);
   p1->SetBottomMargin(0.15);
   p1->SetFrameFillColor(1);
  
   h->Draw("COLZ");//COLZ
   h->GetXaxis()->CenterTitle( kTRUE );
   h->SetTitleOffset( xoff, "x" );
   h->SetXTitle("#chi_{dE/dx}^{2}");
   h->GetYaxis()->CenterTitle( kTRUE );
   h->GetXaxis()->SetLabelSize(0.05);
   h->GetYaxis()->SetLabelSize(0.05);
   h->GetXaxis()->SetTitleSize(0.0);
   h->GetYaxis()->SetTitleSize(0.0);
   h->GetYaxis()->SetNdivisions(505);
   h->GetXaxis()->SetNdivisions(505);
   h->SetTitleOffset( yoff, "y");
   h->SetYTitle("E/p");
   h->SetStats(false);

    leg = new TLegend(0.52,0.04,0.64,0.08);
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
    leg->SetHeader("#chi_{dE/dx}^{2}");
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
    leg->SetHeader("E/p");
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
    lt.DrawLatex(0.5, 0.92,"signal MC");
    lt.DrawLatex(0.75, 0.85,"(b)");
    //f_relevence->Draw("same");
    f_relevence_factor->Draw("same");

    c->cd();
    c->Print("Compare.eps");
    //c->Print("Dalitz.eps");


  
}


    Double_t relevence(Double_t x, Double_t y)
    {
       // if(x>=0&&x<=2.2) return y = 2./20*x+0.6; // x~dedx,  y~eop
       if(x>=0&&x<=4) return y = 0.14*x+0.42; // x~dedx,  y~eop
    }
    
    Double_t relevence_factor(Double_t x, Double_t y)
    {
       // if(x>=0&&x<=2.2) return y = 2./20*x+0.6; // x~dedx,  y~eop
       if(x>=0&&x<=4) return y = 0.18*x+0.32; // x~dedx,  y~eop
    }

