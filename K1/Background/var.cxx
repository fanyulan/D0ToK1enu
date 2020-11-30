void var()
{


TStyle *bes3Style= new TStyle("bes3","bes3 style");

Int_t icol=0;
bes3Style->SetFrameBorderMode(icol);
bes3Style->SetCanvasBorderMode(icol);
bes3Style->SetPadBorderMode(icol);
bes3Style->SetPadColor(icol);
bes3Style->SetCanvasColor(icol);
bes3Style->SetStatColor(icol);
bes3Style->SetTitleFillColor(icol);
bes3Style->SetPalette(1); // set a good color palette 

bes3Style->SetPaperSize(TStyle::kUSLetter);

//bes3Style->SetPadTopMargin(.12);
//bes3Style->SetPadLeftMargin(.15);
//bes3Style->SetPadRightMargin(.08);
//bes3Style->SetPadBottomMargin(.15);

bes3Style->SetPadTopMargin(.02);
bes3Style->SetPadLeftMargin(.12);
bes3Style->SetPadRightMargin(.03);
bes3Style->SetPadBottomMargin(.125);

Int_t font=132;      //times new roman, reg
Double_t tsize=0.05; //should be set between 0.03-0.05, is in units of "% of pad"

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
bes3Style->SetLineStyleString(2,"[12 12]"); // postscript dashes

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

// Begin plots
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
///////////////////////////////////////////////////////Kso pi0 pi0 eta
for(int nvar=0; nvar<7; nvar++) {
   TCut cut_sig = "isTrueTag&&isSignal&&abs(KpTrueID)==321&&abs(pi1TrueID)==211&&abs(pi2TrueID)==211&&abs(elecTrueID)==11&&abs(elecMomID)==421";
   TCut cuta = "isbest==1&&type==1&&mBC<1.874&&mBC>1.858";
   TCut cutc = "elecCharge!=0&&(elecCharge+kaonCharge)==0&&(pionCharge+pion2Charge)==0&&charm==kaonCharge";
   TCut cut0 = "mode==0&&clveto==1&&deltaE>-0.029&&deltaE<0.027";
   TCut cut1 = "mode==1&&deltaE>-0.069&&deltaE<0.038";
   TCut cut3 = "mode==3&&deltaE>-0.031&&deltaE<0.028";
   TCut cut_K1 = "((mcmodea==666||mcmodea==222)&&mcmode1==103)||((mcmodeb==-666||mcmodeb==-222)&&mcmode2==-103)";
   TCut cut_K2 = "mcmode1==104 || mcmode2==-104";
   
   TCut cut_mD_EasPi        = "abs(Umissfit)<0.2";
   TCut cut_costh_lep_pi    = "abs(Umissfit)<0.2&&mD_EasPi<1.81";
   TCut cut_mD_EasPiwPi0    = "abs(Umissfit)<0.2&&mD_EasPi<1.81&&costh_lep_pi<0.94";
   TCut cut_costh_miss_neut = "abs(Umissfit)<0.2&&mD_EasPi<1.81&&costh_lep_pi<0.94&&mD_EasPiwPi0<1.4";
   TCut cut_m_pipi          = "abs(Umissfit)<0.2&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.18*dedxchi_E)"; 
   TCut cut_deltaE_pipzswp  = "abs(Umissfit)<0.2&&mD_EasPi<1.81&&mD_EasPiwPi0<1.4&&costh_lep_pi<0.94&&costh_miss_neut<0.81&&((eop_E-0.32)>0.88*dedxchi_E)&& m_pipi>0.31&&abs(m_pipi-0.497611)>0.01";
   TCut cut_eop_chiE = "abs(Umissfit)<0.2&&mD_EasPi<1.81&&costh_lep_pi<0.94&&mD_EasPiwPi0<1.4&&costh_miss_neut<0.81";

   TCut cut_103 = "(charm==1&&(mcmodeb==-103||mcmodeb==-1032||mcmodeb==-108)) || (charm==-1&&(mcmodea==103||mcmodea==1032||mcmodea==108))";
   TCut cut_102 = "(charm==1&&(mcmodeb==-102||mcmodeb==-104)) || (charm==-1&&(mcmodea==102||mcmodea==104))";//Kpipi0, Kpipi0pi0
   TCut cut_105 = "(charm==1&&(mcmodeb==-455||mcmodeb==-105)) || (charm==-1&&(mcmodea==455||mcmodea==103))"; 
   TCut cut_201 = "(charm==1&&mcmodeb==-201) || (charm==-1&&mcmodea==201)";//Kpi0ev
   TCut cut_715 = "(charm==1&&(mcmodeb==-716||mcmodeb==-715||mcmodeb==-203||mcmodeb==-894||mcmodeb==-603||mcmodeb==588)) || (charm==-1&&(mcmodea==716||mcmodea==715||mcmodea==203||mcmodea==894||mcmodea==603||mcmodea==588))";
   double bins = 100;
   if(nvar==0){     TH1F *h_sig = new TH1F("h_sig","",bins,1.1,2.1);//mD_EasPi
                    TH1F *h_bkg_main  = new TH1F("h_bkg_main","",bins,1.1,2.1);
                    TH1F *h_bkg_other = new TH1F("h_bkg_other","",bins,1.1,2.1); 
                    TH1F *h_data= new TH1F("h_data","",bins,1.1,2.1);
                    } 
                    
   if(nvar==2){     TH1F *h_sig = new TH1F("h_sig","",bins,-1.2,2.2);//mD_EasPiPi0
                    TH1F *h_bkg_main = new TH1F("h_bkg_main","",bins,-1.2,2.2); 
                    TH1F *h_bkg_other = new TH1F("h_bkg_other","",bins,-1.2,2.2); 
                    TH1F *h_data= new TH1F("h_data","",bins,-1.2,2.2);
                    
                    }
   if(nvar==1||nvar==3){    TH1F *h_sig = new TH1F("h_sig","",bins,-1,1);//costh_lep_pi, costh_miss_neut
                            TH1F *h_bkg_main = new TH1F("h_bkg_main","",bins,-1,1); 
                            TH1F *h_bkg_other = new TH1F("h_bkg_other","",bins,-1,1); 
                            TH1F *h_data= new TH1F("h_data","",bins,-1,1);
                             }
   //if(nvar==4){     TH1F *h_sig = new TH1F("h_sig","",bins,0,6);//dedxchi_E
   //                 TH1F *h_bkg = new TH1F("h_bkg","",bins,0,6);            }
   if(nvar==4){     TH1F *h_sig = new TH1F("h_sig","",bins,0.2,0.9);//m_pipi
                    TH1F *h_bkg_main = new TH1F("h_bkg_main","",bins,0.2,0.9);    
                    TH1F *h_bkg_other = new TH1F("h_bkg_other","",bins,0.2,0.9);    
                    TH1F *h_data= new TH1F("h_data","",bins,0.2,0.9);
                    }
   if(nvar==5){     TH1F *h_sig = new TH1F("h_sig","",bins,-0.8,0.4);//deltaE_pipzswp
                    TH1F *h_bkg_main = new TH1F("h_bkg_main","",bins,-0.8,0.4);      
                    TH1F *h_bkg_other = new TH1F("h_bkg_other","",bins,-0.8,0.4);      
                    TH1F *h_data= new TH1F("h_data","",bins,-0.8,0.4);
                    }
   if(nvar==6){     TH1F *h_sig = new TH1F("h_sig","",bins,0,1.2);//eop
                    TH1F *h_bkg_k3pi = new TH1F("h_bkg_k3pi","",bins,0,1.2);      
                    TH1F *h_bkg_k3pipi0 = new TH1F("h_bkg_k3pipi0","",bins,0,1.2);      
                    TH1F *h_bkg_other = new TH1F("h_bkg_other","",bins,0,1.2);      
                    TH1F *h_data= new TH1F("h_data","",bins,0,1.2);
                    }

   TChain t("tagD");
   t.Add("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/new/inclusiveMC/umiss2fit/reduce_truth/root/truth*.root");      
   
   TChain tdata("tagD");  
   tdata.Add("/scratchfs/bes/fangyl/data1/data/all/*.root");

   if(nvar==0){t.Draw("mD_EasPi>>h_sig",(cut0||cut1||cut3)&&cuta&&cutc&&cut_sig&&cut_mD_EasPi&&cut_K1);//sig 
               t.Draw("mD_EasPi>>h_bkg_main",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_mD_EasPi&&cut_103);  
               t.Draw("mD_EasPi>>h_bkg_other",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_mD_EasPi&&!cut_103);   
               tdata.Draw("mD_EasPi>>h_data",(cut0||cut1||cut3)&&cuta&&cutc&&cut_mD_EasPi);               
               }
   if(nvar==2){t.Draw("mD_EasPiwPi0>>h_sig",(cut0||cut1||cut3)&&cuta&&cutc&&(cut_sig&&cut_K1)&&cut_mD_EasPiwPi0);//sig
               t.Draw("mD_EasPiwPi0>>h_bkg_main",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_mD_EasPiwPi0&&cut_105); 
               t.Draw("mD_EasPiwPi0>>h_bkg_other",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_mD_EasPiwPi0&&!cut_105); 
               tdata.Draw("mD_EasPiwPi0>>h_data",(cut0||cut1||cut3)&&cuta&&cutc&&cut_mD_EasPiwPi0);
               }
   if(nvar==1){t.Draw("costh_lep_pi>>h_sig",(cut0||cut1||cut3)&&cuta&&cutc&&(cut_sig&&cut_K1)&&cut_costh_lep_pi);//sig
               t.Draw("costh_lep_pi>>h_bkg_main",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_costh_lep_pi&&cut_102); 
               t.Draw("costh_lep_pi>>h_bkg_other",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_costh_lep_pi&&!cut_102); 
               tdata.Draw("costh_lep_pi>>h_data",(cut0||cut1||cut3)&&cuta&&cutc&&cut_costh_lep_pi);
               }
   if(nvar==3){t.Draw("costh_miss_neut>>h_sig",(cut0||cut1||cut3)&&cuta&&cutc&&(cut_sig&&cut_K1)&&cut_costh_miss_neut);//sig
               t.Draw("costh_miss_neut>>h_bkg_main",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_costh_miss_neut&&cut_105); 
               t.Draw("costh_miss_neut>>h_bkg_other",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_costh_miss_neut&&!cut_105); 
               tdata.Draw("costh_miss_neut>>h_data",(cut0||cut1||cut3)&&cuta&&cutc&&cut_costh_miss_neut);
               }
   //if(nvar==4){t.Draw("dedxchi_E>>h_sig",(cut0||cut1||cut3)&&cuta&&cutc&&cut_sig&&cut_dedxchi_E);//sig
   //            t.Draw("dedxchi_E>>h_bkg",(cut0||cut1||cut3)&&cuta&&cutc&&!cut_sig&&cut_dedxchi_E);             }
   //if(nvar==5){t.Draw("eop_E>>h_sig",(cut0||cut1||cut3)&&cuta&&cutc&&cut_sig&&cut_eop_E);//sig
   //            t.Draw("eop_E>>h_bkg",(cut0||cut1||cut3)&&cuta&&cutc&&!cut_sig&&cut_eop_E);                     }
   if(nvar==4){t.Draw("m_pipi>>h_sig",(cut0||cut1||cut3)&&cuta&&cutc&&(cut_sig&&cut_K1)&&cut_m_pipi);//sig
               t.Draw("m_pipi>>h_bkg_main",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_m_pipi&&cut_201);
               t.Draw("m_pipi>>h_bkg_other",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_m_pipi&&!cut_201);
               tdata.Draw("m_pipi>>h_data",(cut0||cut1||cut3)&&cuta&&cutc&&cut_m_pipi);
               }
   if(nvar==5){t.Draw("deltaE_pipzswp>>h_sig",(cut0||cut1||cut3)&&cuta&&cutc&&(cut_sig&&cut_K1)&&cut_deltaE_pipzswp);//sig
               t.Draw("deltaE_pipzswp>>h_bkg_main",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_deltaE_pipzswp&&cut_715); 
               t.Draw("deltaE_pipzswp>>h_bkg_other",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_deltaE_pipzswp&&!cut_715); 
               tdata.Draw("deltaE_pipzswp>>h_data",(cut0||cut1||cut3)&&cuta&&cutc&&cut_deltaE_pipzswp);
               }
  
   if(nvar==6){t.Draw("eop_E>>h_sig",(cut0||cut1||cut3)&&cuta&&cutc&&(cut_sig&&cut_K1)&&cut_eop_chiE);//sig
               t.Draw("eop_E>>h_bkg_k3pi",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_eop_chiE&&cut_103);  
               t.Draw("eop_E>>h_bkg_k3pipi0",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_eop_chiE&&cut_105);  
               t.Draw("eop_E>>h_bkg_other",(cut0||cut1||cut3)&&cuta&&cutc&&!(cut_sig&&cut_K1)&&!cut_K2&&cut_eop_chiE&&!cut_105&&!cut_103); 
               tdata.Draw("eop_E>>h_data",(cut0||cut1||cut3)&&cuta&&cutc&&cut_eop_chiE);
               }


     double ndata1= h_data->GetEntries();
     double nsig = h_sig->GetEntries();
     cout<<nsig<<endl;
     //double ratio = 1756./1275; 
     //cout<<ratio<<endl;
     TCanvas *c=new TCanvas("c","",0,0,800,600);     
     TPad *p =new TPad("c1","",0.,0.0,1,1);
     p->Draw();
     p->cd();
     p->SetFillStyle(4000);
     p->Range(0,0,1,1);
     p->SetLeftMargin(0.15);
     p->SetRightMargin(0.06);
     p->SetTopMargin(0.1);
     p->SetBottomMargin(0.15);
     p->SetFrameFillColor(1);

//if(nvar==5){  
//     h_sig->SetLineColor(kRed);
//     h_sig->Draw();
//     h_bkg->SetLineColor(kBlue);
//     h_bkg->Draw("same");      }
//else {
     
     
     //h_bkg_other->Scale(1.0/10.8);
     //h_bkg_main->Scale(1.0/10.8);
     //h_sig->Scale(1.0/10.8);
     //h_bkg_main->Add(h_bkg_other);
     //h_sig->Add(h_bkg_main);
     if(nvar!=6){
        double nbkg_main = h_bkg_main->GetEntries();
        double nbkg_other = h_bkg_other->GetEntries();
        double ndata = h_data->GetEntries(); 
        double nsig = h_sig->GetEntries();
        double ntot = nbkg_main+nbkg_other+nsig;
        h_bkg_main->Scale(ndata/ntot);
        h_bkg_other->Scale(ndata/ntot);
        h_sig-> Scale(ndata/ntot);

        h_bkg_main->Add(h_bkg_other,1);
        h_bkg_main->Draw();
        h_bkg_main->SetLineWidth(2);
        h_bkg_main->SetFillColor(kYellow);
        h_bkg_main->SetLineColor(0);
        h_bkg_main->GetYaxis()->SetRangeUser(0,1.55*h_data->GetMaximum());
        h_bkg_other->Draw("same");
        h_bkg_other->SetLineWidth(2);
        h_bkg_other->SetFillColor(kGray);
        h_bkg_other->SetLineColor(0);
        h_bkg_other->SetMarkerStyle(20);
        h_bkg_other->SetMarkerSize(0.8);
        h_bkg_other->GetYaxis()->SetRangeUser(0,1.55*h_data->GetMaximum());
        
        h_sig->Add(h_bkg_main,1);
        h_sig->SetLineColor(kRed);
        h_sig->Draw("same");
     }
     else{    
        double nbkg_k3pi = h_bkg_k3pi->GetEntries();
        double nbkg_k3pipi0 = h_bkg_k3pipi0->GetEntries();
        double nbkg_other = h_bkg_other->GetEntries();
        double ndata = h_data->GetEntries(); 
        double nsig = h_sig->GetEntries();
        double ntot = nbkg_k3pi+nbkg_k3pipi0+nbkg_other+nsig;
        h_bkg_k3pi->Scale(ndata/ntot);
        h_bkg_k3pipi0->Scale(ndata/ntot);
        h_bkg_other->Scale(ndata/ntot);
        h_sig-> Scale(ndata/ntot);

        h_bkg_k3pipi0->Add(h_bkg_other,1);
        h_bkg_k3pi->Add(h_bkg_k3pipi0,1);

        h_bkg_k3pi->Draw();
        h_bkg_k3pi->SetFillColor(kYellow);
        h_bkg_k3pi->SetLineColor(0);
        h_bkg_k3pi->SetLineWidth(2);
        h_bkg_k3pi->GetYaxis()->SetRangeUser(0,1.55*h_data->GetMaximum());
        h_bkg_k3pipi0->Draw("same");
        h_bkg_k3pipi0->SetFillColor(kCyan-3);
        h_bkg_k3pipi0->SetLineColor(0);
        
        h_bkg_other->Draw("same");
        h_bkg_other->SetLineWidth(2);
        h_bkg_other->SetFillColor(kGray);
        h_bkg_other->SetLineColor(0);
        h_bkg_other->SetMarkerStyle(20);
        h_bkg_other->SetMarkerSize(0.8);
        h_bkg_other->GetYaxis()->SetRangeUser(0,1.55*h_data->GetMaximum());
        
        h_sig->Add(h_bkg_k3pi,1);
        h_sig->SetLineColor(kRed);
        h_sig->Draw("same");

     } 
         
     h_sig->GetXaxis()->CenterTitle( kTRUE );
     h_sig->SetTitleOffset( xoff, "x" );
     h_sig->GetYaxis()->CenterTitle( kTRUE );
     h_sig->GetXaxis()->SetLabelSize(0.05);
     h_sig->GetYaxis()->SetLabelSize(0.05);
     h_sig->GetXaxis()->SetTitleSize(0.0);
     h_sig->GetYaxis()->SetTitleSize(0.0);
     h_sig->GetYaxis()->SetNdivisions(505);
     h_sig->GetXaxis()->SetNdivisions(505);
     h_sig->SetTitleOffset( yoff, "y");
     
     h_data->Draw("E same");
     h_data->SetLineWidth(2);
     h_data->SetMarkerStyle(20);
     h_data->SetMarkerSize(0.8);
     
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

if(nvar==0)    leg->SetHeader("M_{K^{-}#pi^{+}#pi^{-}e(#pi)^{+}} (GeV/#font[12]{c}^{2})");//mD
if(nvar==2)     leg->SetHeader("M_{K^{-}#pi^{+}#pi^{-}e(#pi)^{+}#pi^{0}} (GeV/#font[12]{c}^{2})");
if(nvar==1)     leg->SetHeader("cos(e^{+}#pi^{-})");
if(nvar==3)     leg->SetHeader("cos(miss_#gamma))");
//if(nvar==4)     leg->SetHeader("#chi_{dE/dx}^{2}");
//if(nvar==5)     leg->SetHeader("E/p");
if(nvar==4)     leg->SetHeader("M_{#pi^{+}#pi^{-}} (GeV/#font[12]{c}^{2})");
if(nvar==5)     leg->SetHeader("#Delta E_{(K#pi)_{tag}#pi_{sig}} (GeV)");
if(nvar==6)     leg->SetHeader("E/p");;
     leg->Draw();

     leg1 = new TLegend(0.015,0.45,0.055,0.60);
     leg1->SetTextAlign(22);
     leg1->SetTextFont(22);

     leg1->SetBorderSize(0);
     leg1->SetTextSize(0.05);
     leg1->SetTextColor(1);
     leg1->SetLineColor(0);
     leg1->SetFillColor(0);
     leg1->SetLineWidth(0);
     leg1->SetTextAngle(90); 
if(nvar==0)     leg1->SetHeader("Events / (0.01Gev)/c^{2}");//bin size
if(nvar==1||nvar==3)     leg1->SetHeader("Events / (0.02)");
if(nvar==2)     leg1->SetHeader("Events / (0.034Gev)/c^{2}"); 
//if(nvar==4)     leg1->SetHeader("Events (/0.06)");
//if(nvar==5)     leg1->SetHeader("Events (/0.016)");
if(nvar==4)     leg1->SetHeader("Events / (0.007Gev/c^{2})");
if(nvar==5)     leg1->SetHeader("Events / (0.012Gev)");
if(nvar==6)     leg1->SetHeader("Events / (0.012)");

     leg1->Draw();


if(nvar==0)     TArrow *ar1 = new TArrow(1.81,0.6*h_bkg_other->GetMaximum(),1.81,0.1*h_bkg_other->GetMaximum(),0.01,"|>");//mD1 
if(nvar==1)     TArrow *ar1 = new TArrow(0.94,0.6*h_bkg_other->GetMaximum(),0.94,0.1*h_bkg_other->GetMaximum(),0.01,"|>");//cos
if(nvar==2)     TArrow *ar1 = new TArrow(1.4,0.6*h_bkg_other->GetMaximum(),1.4,0.1*h_bkg_other->GetMaximum(),0.01,"|>");//mD1 
if(nvar==3)     TArrow *ar1 = new TArrow(0.81,0.8*h_bkg_other->GetMaximum(),0.81,0.1*h_bkg_other->GetMaximum(),0.01,"|>");//cos
//if(nvar==4)     TArrow *ar1 = new TArrow(2.6,0.8*h_bkg_other->GetMaximum(),2.6,0.1*h_bkg_other->GetMaximum(),0.01,"|>");//dedx
//if(nvar==5)     TArrow *ar1 = new TArrow(0.65,0.9*h_bkg_other->GetMaximum(),0.65,0.1*h_bkg_other->GetMaximum(),0.01,"|>");//eop
if(nvar==4){    TArrow *ar1 = new TArrow(0.31,0.9*h_bkg_other->GetMaximum(),0.31,0.2*h_bkg_other->GetMaximum(),0.01,"|>");//mpipi  
                TArrow *ar2 = new TArrow(0.487611,0.9*h_bkg_other->GetMaximum(),0.487611,0.1*h_bkg_other->GetMaximum(),0.01,"|>");//mpipi
                TArrow *ar3 = new TArrow(0.507611,0.9*h_bkg_other->GetMaximum(),0.507611,0.1*h_bkg_other->GetMaximum(),0.01,"|>");//mpipi
                TArrow *ar1_1 = new TArrow(0.4,0.7*h_bkg_other->GetMaximum(),0.3,0.5*h_bkg_other->GetMaximum(),0.01,"|>");//mpipi ~ 
                TArrow *ar2_1 = new TArrow(0.6,0.6*h_bkg_other->GetMaximum(),0.495,0.4*h_bkg_other->GetMaximum(),0.01,"|>");   }//mpipi ~ 
if(nvar==5){    TArrow *ar5 = new TArrow(0.014,0.9*h_bkg_other->GetMaximum(),0.014,0.1*h_bkg_other->GetMaximum(),0.01,"|>");//deltaEkpipi
                TArrow *ar4 = new TArrow(-0.014,0.9*h_bkg_other->GetMaximum(),-0.014,0.1*h_bkg_other->GetMaximum(),0.01,"|>"); }//deltaEkpipi                
if(nvar==6)     TArrow *ar1 = new TArrow(0.42,0.7*h_bkg_other->GetMaximum(),0.42,0.04*h_bkg_other->GetMaximum(),0.01,"|>");//mD1 

if(nvar!=5&&nvar!=6)  { 
     ar1->SetLineColor(kBlack);
     ar1->SetFillColor(kBlack);
     ar1->SetLineWidth(2);
     ar1->Draw();
     }
if(nvar==6)  { 
     ar1->SetLineColor(kBlue);
     ar1->SetFillColor(kBlue);
     ar1->SetLineWidth(2);
     ar1->Draw();
     }
if(nvar==4){
     ar1->SetLineColor(kBlue);
     ar1->SetFillColor(kBlue);
     ar1->SetLineWidth(2);
     ar1->Draw();
     ar2->SetLineColor(kBlue);
     ar2->SetFillColor(kBlue);
     ar2->SetLineWidth(2);
     ar2->Draw();
     ar3->SetLineColor(kBlue);
     ar3->SetFillColor(kBlue);
     ar3->SetLineWidth(2);
     ar3->Draw();
     ar1_1->SetLineColor(kGreen+3);
     ar1_1->SetFillColor(kGreen+3);
     ar1_1->SetLineWidth(2);
     ar1_1->Draw();
     ar2_1->SetLineColor(kGreen+3);//green
     ar2_1->SetFillColor(kGreen+3);
     ar2_1->SetLineWidth(2);
     ar2_1->Draw();
     }
     if(nvar==5) {
     ar5->SetLineColor(kBlue);
     ar5->SetFillColor(kBlue);
     ar5->SetLineWidth(2);
     ar5->Draw();
     ar4->SetLineColor(kBlue);
     ar4->SetFillColor(kBlue);
     ar4->SetLineWidth(2);
     ar4->Draw();
     }
     if(nvar!=5&&nvar!=1){ leg = new TLegend(0.60,0.7,0.91,0.89);}//(x0.y0.x1,y1)
     if(nvar==5){ leg = new TLegend(0.30,0.7,0.61,0.89);}//(x0.y0.x1,y1)
     if(nvar==1){ leg = new TLegend(0.30,0.65,0.67,0.89);}//(x0.y0.x1,y1)
     
     leg->SetBorderSize(0);
     leg->SetTextFont(22);
     leg->SetLineColor(0);
     leg->SetFillColor(0);
     leg->AddEntry(h_data,"data");
     if(nvar==0)    leg->AddEntry(h_bkg_main,"D^{0}#rightarrowK^{-}#pi^{+}#pi^{-}#pi^{+} background","f");
     if(nvar==1)    leg->AddEntry(h_bkg_main,"D^{0}#rightarrowK^{-}#pi^{+}#pi^{0} and D^{0}#rightarrowK^{-}#pi^{+}#pi^{0}#pi^{0} background","f");
     if(nvar==2||nvar==3)    leg->AddEntry(h_bkg_main,"D^{0}#rightarrowK^{-}#pi^{+}#pi^{-}#pi^{+}#pi^{0} background","f");
     if(nvar==4)     leg->AddEntry(h_bkg_main,"D^{0}#rightarrowK^{-}#pi^{0}e^{+}#nu_{e} background","f");
     if(nvar==5)     leg->AddEntry(h_bkg_main,"D^{+}#rightarrowK^{-}#pi^{+}X background","f");
     if(nvar==6)    {leg->AddEntry(h_bkg_k3pi,"D^{0}#rightarrowK^{-}#pi^{+}#pi^{-}#pi^{+} background","f");
                     leg->AddEntry(h_bkg_k3pipi0,"D^{0}#rightarrowK^{-}#pi^{+}#pi^{-}#pi^{+}#pi^{0} background","f");}
     leg->AddEntry(h_bkg_other,"other background","f");
     leg->AddEntry(h_sig,"signal","f");
     leg->Draw();

//     TLatex lt;
//     lt.SetNDC();
//     lt.SetTextAngle(0);
//     lt.SetTextSize(0.05);
//     lt.DrawLatex(0.3, 0.82,  Form("N_{sig} = %.0f ",h_sig->GetEntries()-h_bkg->GetEntries()));
//     lt.DrawLatex(0.3, 0.76,  Form("N_{bkg} = %.0f ",h_bkg->GetEntries()));
    // lt.DrawLatex(0.72, 0.92,  "D^{0}#rightarrowK^{-}#pi^{+}#pi^{-}e^{+}#nu_{e}");
      
if(nvar==4){
     TLatex ltm;
     ltm.SetNDC();
     ltm.SetTextAngle(0);
     ltm.SetTextSize(0.05);
     ltm.DrawLatex(0.3, 0.68,"#pi^{0}#rightarrow#gamma e^{+}e^{-}");//mpipi ~
     ltm.DrawLatex(0.6, 0.6,"K_{s}");//mpipi ~
     }
     TLatex ltn;
     ltn.SetNDC();
     ltn.SetTextAngle(0);
     ltn.SetTextSize(0.05);
if(nvar==0)     ltn.DrawLatex(0.2,  0.82,"(b)");
if(nvar==1)     ltn.DrawLatex(0.2,  0.82,"(b)");
if(nvar==2)     ltn.DrawLatex(0.22, 0.82,"(b)");
if(nvar==3)     ltn.DrawLatex(0.2,  0.82,"(b)");
if(nvar==4)     ltn.DrawLatex(0.4,  0.82,"(a)");
if(nvar==5)     ltn.DrawLatex(0.2,  0.82,"(b)");
if(nvar==6)     ltn.DrawLatex(0.2,  0.82,"(c)");

     c->cd();

     c->SaveAs(Form("./var%d.eps",nvar));

}
}
