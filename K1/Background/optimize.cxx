double value_x[20];
double value_y[20];
double value_max1[20];
double value_max2[20];
double D2_fig[200][200];
double D2_sb[200][200];
int num;
double cut1[20][200],cut2[20][200],s[20][200][200];
double sb_cut,sb_max;
TString CUT_final;
TString CUT1[20];

void optimize(){
gSystem->Load("libRooFit");
using namespace RooFit;
gStyle->SetLabelSize(0.07,"xyz");
gStyle->SetNdivisions(507,"xyz");
gROOT->SetStyle("Plain");
gStyle->SetLabelSize(0.06,"xyz");
gStyle->SetNdivisions(405,"xyz");
gStyle->SetPadTopMargin(.10);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadRightMargin(.05);
gStyle->SetPadBottomMargin(.15);
gStyle->SetTitleSize(0.06,"xyz");
gStyle->SetOptTitle(0);
gStyle->SetMarkerSize(0.5);
gStyle->SetTitle("");
gStyle->SetOptStat(0);
gStyle->SetEndErrorSize(0);


TFile *f1=TFile::Open("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/bkg_analysis/opt/figure/sig_opt.root");///sig
TFile *f2=TFile::Open("/besfs/groups/psipp/psippgroup/public/fangyl/cocktail/script/bkg_analysis/opt/figure/bkg_opt.root");///bkg

Optimize ttest;
ttest.SetTreeName("tagD"); //tree name
ttest.AddFile(f1,f2);//
ttest.AddCut("mD_EasPi","M_{K^{-}#pi^{+}#pi^{-}e(#pi)^{+}}(GeV/#font[12]{c}^{2})",50,1.5,2,"<");///name of parameter; X title; bin; cut_min; cut_max; m:<, p:>
ttest.AddCut("costh_lep_pi","cos(e^{+},#pi^{-})",50,0.5,1,"<");///name of parameter; X title; bin; cut_min; cut_max; m:<, p:>
ttest.AddCut("mD_EasPiwPi0","M_{K^{-}#pi^{-}#pi^{+}e(#pi)^{+}#pi^{0}}(GeV/#font[12]{c}^{2})",50,1.4,2,"<");
ttest.AddCut("costh_miss_neut","cos(miss,#gamma)",50,0.5,1,"<");
//ttest.AddCut("eop_E","dedxchi_E","Slope","Intercept",50,0,1,50,0,1,">");
ttest.AddCut("m_pipi","M_{#pi^{+}#pi^{-}} (GeV/#font[12]{c}^{2})",50,0.29,0.39,">");
ttest.AddCut("abs(deltaE_pipzswp)","#Delta E_{(K#pi)_{tag}#pi_{sig}}",50,0,0.1,">");
ttest.AddCut("eop_E","dedxchi_E","Slope","Intercept",50,0,1,50,0,1,">");

//ttest.SetXmax(0,0.5);      
//ttest.SetArrowRange(0,0.95,0.99);
//ttest.SetArrowRange(1,0.9982,0.999);
//ttest.SetArrowRange(2,0.9993,0.9998);
//ttest.SetArrowColor(0,4);
//ttest.SetArrowColor(1,5);
//ttest.SetArrowColor(2,6);
//ttest.SetArrowColor(3,7);
//ttest.SetAxisRange(0,100,300);
//ttest.SetAxisRange(1,150,250);
//ttest.SetAxisRange(2,-1,400);

ttest.opt();
ttest.Draw(0,0);
ttest.Draw(1,0);
ttest.Draw(2,0);
ttest.Draw(3,0);
ttest.Draw(4,0);
ttest.Draw(5,0);
ttest.Draw(6,0);
ttest.Draw(6,1);

TCanvas *resf1= new TCanvas("resf1"," ",800,600);
TGraph2D *dt1 = new TGraph2D();

 int l=0;
 for(int j=0;j<100;j++){
   for(int k=0;k<100;k++){
       dt1->SetPoint(l,D2_fig[j][0],D2_fig[0][k],D2_sb[j][k]);
       l++;
   }
  }
    dt1->SetMarkerSize(0.5);
    dt1->SetMarkerColor(2);
    dt1->SetMarkerStyle(22);
    
    dt1->SetTitle(";E/p;#chi_{dE/dx}^{2};S/#sqrt{S+B}");
    dt1->Draw("surf1");
    
    TLatex lx;
    lx.SetNDC();
    lx.SetTextAngle(0);
    lx.SetTextSize(0.05);
    lx.DrawLatex(0.8,0.8,"(a)");
    
    resf1->Print("./figure/test.eps");
}

class Optimize{
public:
void ~Optimize();
void AddFile(TFile *file1,TFile *file2);
void SetTreeName(char *tree);
void AddCut(string cut,string name1,int bin,double cut_min,double cut_max,string pm);
void AddCut(string cut1,string cut2,string name1,string name2,int bin1,double cut_min1,double cut_max1,int bin2,double cut_min2,double cut_max2,string pm);
void SetXmax(int N,double x1);
void SetArrowRange(int A,double A1,double A2);
void SetAxisRange(int M,double M1,double M2);
void SetArrowColor(int B,int B1);
void fun_opt1(int i,int bin_temp);
void opt();
void Draw();
protected:
TTree *_sig;
TTree *_bkg;
char *_tree;
std::vector<TString> _title1;
std::vector<TString> _title2;
std::vector<string> _cut1;
std::vector<string> _cut2;
std::vector<string> _pm;
std::vector<double> _cut_min1;
std::vector<double> _cut_max1;
std::vector<double> _cut_min2;
std::vector<double> _cut_max2;
std::vector<int> _bin1;
std::vector<int> _bin2;
double _D_cut1[20][200];
double _D_cut2[20][200];
double _D_s[20][200][200];
double _max1[20];
double _max2[20];
double _each_max[20];
std::vector<double> _x_max;
std::vector<int> _N;
std::vector<int> _A;
std::vector<double> _A_min;
std::vector<double> _A_max;
std::vector<int> _M;
std::vector<double> _M_min;
std::vector<double> _M_max;
std::vector<int> _B;
std::vector<int> _B_col;
};

void Optimize::AddFile(TFile *file1,TFile *file2){
_sig=(TTree*)file1->Get(_tree);
_bkg=(TTree*)file2->Get(_tree);
}
void Optimize::SetXmax(int N,double x1){
_x_max.push_back(x1);
_N.push_back(N);
}

void Optimize::SetTreeName(char *tree){
_tree=tree;
}

void Optimize::AddCut(string cut,string name1,int bin,double cut_min,double cut_max,string pm){
_cut1.push_back(cut);
_cut2.push_back("");
_pm.push_back(pm);
_cut_min1.push_back(cut_min);
_cut_max1.push_back(cut_max);
_cut_min2.push_back(0);
_cut_max2.push_back(0);
_bin1.push_back(bin);
_bin2.push_back(0);
_title1.push_back(name1);
_title2.push_back("");
}
void Optimize::AddCut(string cut1,string cut2,string name1,string name2,int bin1,double cut_min1,double cut_max1,int bin2,double cut_min2,double cut_max2,string pm){
_cut1.push_back(cut1);
_cut2.push_back(cut2);
_pm.push_back(pm);
_cut_min1.push_back(cut_min1);
_cut_max1.push_back(cut_max1);
_cut_min2.push_back(cut_min2);
_cut_max2.push_back(cut_max2);
_bin1.push_back(bin1);
_bin2.push_back(bin2);
_title1.push_back(name1);
_title2.push_back(name2);
}

void Optimize::SetArrowRange(int A,double A1,double A2){
_A_min.push_back(A1);
_A_max.push_back(A2);
_A.push_back(A);
}
void Optimize::SetAxisRange(int M,double M1,double M2){
_M_min.push_back(M1);
_M_max.push_back(M2);
_M.push_back(M);
}
void Optimize::SetArrowColor(int B,int B1){
_B_col.push_back(B1);
_B.push_back(B);
}

void Optimize::opt(){
num=_cut1.size();
for(int i=0;i<20;i++){
for(int j=0;j<200;j++){
for(int k=0;k<200;k++){
	cut1[i][j]=0.;
	cut2[i][k]=0.;
	s[i][j][k]=0.;
	}
}
}
sb_cut=0,sb_max=0;
memset(value_x,0,20);
memset(value_y,0,20);
memset(value_max1,0,20);
memset(value_max2,0,20);
//double each_max[20];memset(each_max,0,20);
//double cut[20];memset(cut,0,20);



for(int i=0;i<num;i++){
if(_pm[i]=="<"){value_x[i]=999;value_y[i]=999;}
else {value_x[i]=-999;value_y[i]=-999;}
}

for(int i=0;i<num;i++){

stringstream tem;
if(_bin2[i]==0)tem<<_cut1[i]<<_pm[i]<<value_x[i];
if(_bin2[i]!=0)tem<<"("<<_cut1[i]<<_pm[i]<<value_x[i]<<"*"<<_cut2[i]<<"+("<<value_y[i]<<"))";
CUT1[i]=tem.str();
}

for(int x=0;x<5;x++){


for(int i=0;i<num;i++){ 
int bin_temp1=_bin1[i];	
int bin_temp2=_bin2[i];	
sb_cut=0;

if(_bin2[i]==0)fun_opt1(i,bin_temp1);
if(_bin2[i]!=0)fun_opt2(i,bin_temp1,bin_temp2);

stringstream tem1;
if(_bin2[i]==0)tem1<<_cut1[i]<<_pm[i]<<value_max1[i];
if(_bin2[i]!=0)tem1<<"("<<_cut1[i]<<_pm[i]<<value_max1[i]<<"*"<<_cut2[i]<<"+("<<value_max2[i]<<"))";
CUT1[i]=tem1.str();
}
}

for(int i=0;i<num;i++){
for(int j=0;j<=_bin1[i];j++){
for(int k=0;k<=_bin2[i];k++){
	_D_cut1[i][j]=cut1[i][j];
	_D_cut2[i][k]=cut2[i][k];
	_D_s[i][j][k]=s[i][j][k];
}
}
}

for(int i=0;i<num;i++){
ofstream out1(Form("figure/cut%d.dat",i));
if(_bin2[i]==0)out1<<Form("%.2f",value_max1[i])<<endl;
if(_bin2[i]!=0)out1<<Form("%.2f	%.2f",value_max1[i],value_max2[i])<<endl;
out1.close();
_max1[i]=value_max1[i];
_max2[i]=value_max2[i];
}
cout<<"best cut: "<<CUT_final<<endl;
ofstream out_cut("./figure/cut_out.dat");
out_cut<<CUT_final<<"\n";

}

void Optimize::Draw(int i,int x){
	gROOT->SetStyle("Plain");
	gStyle->SetLabelSize(0.06,"xyz");
	gStyle->SetNdivisions(405,"xyz");
	gStyle->SetPadTopMargin(0.10);
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.05);
	gStyle->SetPadBottomMargin(0.15);
	gStyle->SetTitleSize(0.06,"xyz");
	gStyle->SetOptTitle(0);
	gStyle->SetMarkerSize(0.5);
	gStyle->SetTitle("");
	gStyle->SetOptStat(0);
	gStyle->SetEndErrorSize(0);
	using namespace RooFit;

	TCanvas *c1 = new TCanvas("c1","",0,0,800,600);
	//c1->Range(0,0,1,1);
	//c1->SetBorderSize(2);
	//c1->SetLeftMargin(0.05);
	//c1->SetRightMargin(0.05);
	//c1->SetTopMargin(0.05);
	//c1->SetBottomMargin(0.05);
	//c1->SetFrameFillColor(0);

	TPad *c1_1 = new TPad("c1_1", "",0.,0.,1,1);
	c1_1->Draw();
	c1_1->cd();
	c1_1->Range(0,0,1,1);
	c1_1->SetFillStyle(4000);
	//c1_1->SetBorderMode(0);
	//c1_1->SetBorderSize(2);
	//c1_1->SetTickx();
	//c1_1->SetTicky();
	c1_1->SetLeftMargin(0.15);
	c1_1->SetRightMargin(0.06);
	c1_1->SetTopMargin(0.1);
	c1_1->SetBottomMargin(0.15);
	c1_1->SetFrameFillColor(0);

double cut_tem[200]={0};
double s_tem[200]={0};

int tem_x,tem_y;
if(_pm[i]=="<"){
	tem_x=(_cut_max1[i]-value_max1[i])/((_cut_max1[i]-_cut_min1[i])/_bin1[i]);
	tem_y=(_cut_max2[i]-value_max2[i])/((_cut_max2[i]-_cut_min2[i])/_bin2[i]);
}
if(_pm[i]==">"){
	tem_x=(value_max1[i]-_cut_min1[i])/((_cut_max1[i]-_cut_min1[i])/_bin1[i]);
	tem_y=(value_max2[i]-_cut_min2[i])/((_cut_max2[i]-_cut_min2[i])/_bin2[i]);
}

if(x==0){
for(int j=0;j<=_bin1[i];j++){
cut_tem[j]=_D_cut1[i][j];
if(_bin2[i]==0)s_tem[j]=_D_s[i][j][0];
if(_bin2[i]!=0)s_tem[j]=_D_s[i][j][tem_y];
}
}
if(x==1){
for(int j=0;j<=_bin2[i];j++){
cut_tem[j]=_D_cut2[i][j];
s_tem[j]=_D_s[i][tem_x][j];
}
}

double r_max1[20],r_min1[20];for(int xx=0;xx<20;xx++){r_max1[xx]=-1;r_min1[xx]=-1;}
for(int tem=0;tem<_M.size();tem++){
r_min1[_M[tem]]=_M_min[tem];
r_max1[_M[tem]]=_M_max[tem];
}

if(x==0)gr1 = new TGraph(_bin1[i]+1,cut_tem,s_tem);
if(x==1)gr1 = new TGraph(_bin2[i]+1,cut_tem,s_tem);
gr1->SetLineColor(2);
gr1->SetLineWidth(4);
gr1->SetMarkerColor(4);
gr1->SetMarkerStyle(21);
if(r_max1[i]!=-1)gr1->SetMaximum(r_max1[i]);
if(r_min1[i]!=-1)gr1->SetMinimum(r_min1[i]);
gr1->Draw("AP");

for(int tem=0;tem<_N.size();tem++){
_max1[_N[tem]]=_x_max[tem];
}

double r_tem1[20],r_tem2[20];for(int xx=0;xx<20;xx++){r_tem1[xx]=0.8;r_tem2[xx]=0.9;}
for(int tem=0;tem<_A.size();tem++){
r_tem1[_A[tem]]=_A_min[tem];
r_tem2[_A[tem]]=_A_max[tem];
}

int c_tem[20];for(int xx=0;xx<20;xx++){c_tem[xx]=2;}
for(int tem=0;tem<_B.size();tem++){
c_tem[_B[tem]]=_B_col[tem];
}

TArrow *ar1;
TArrow *ar1_1;

//if(x==0)	  ar1 = new TArrow(_max1[i], r_tem1[i]*_each_max[i], _max1[i], r_tem2[i]*_each_max[i], 0.02, "|>");
//if(x==1)	  ar1 = new TArrow(_max2[i], r_tem1[i]*_each_max[i], _max2[i], r_tem2[i]*_each_max[i], 0.02, "|>");
//cout<<i<<"\t"<<_max[i]<<"\t"<<r_tem1[i]<<"\t"<<r_tem2[i]<<"\t"<<_each_max[i]<<endl;
if(i==0)  ar1 = new TArrow(1.81,12.0,1.81,20.0,0.02,"|>");
if(i==1)  ar1 = new TArrow(0.94,25.0,0.94,26.0,0.02,"|>");
if(i==2)  ar1 = new TArrow(1.40,26.0,1.40,26.4,0.02,"|>");
if(i==3)  ar1 = new TArrow(0.81,24.0,0.81,26.0,0.02,"|>");
if(i==4)  ar1 = new TArrow(0.31,19.7,0.31,20.4,0.02,"|>");
if(i==5)  ar1 = new TArrow(0.012,25.8,0.012,26.4,0.02,"|>");
if(i==6)  ar1 = new TArrow(0.012,25.5,0.012,26.3,0.02,"|>");
//          ar1_1=new TArrow(0.42,10,0.31,20,0.02,"|>");

//if(i!=6){  
              ar1->SetFillColor(c_tem[i]);
	          ar1->SetLineColor(c_tem[i]);
	          ar1->SetLineWidth(3);
	          ar1->SetFillStyle(1001);
	          ar1->Draw();      
	          c1_1->Modified();
	          c1->cd();
//         }
/*if(i==6){     ar1->SetFillColor(c_tem[i]);
	          ar1->SetLineColor(c_tem[i]);
	          ar1->SetLineWidth(3);
	          ar1->SetFillStyle(1001);
	          ar1->Draw();      
	          ar1_1->SetFillColor(c_tem[i]);
              ar1_1->SetLineColor(c_tem[i]);
              ar1_1->SetLineWidth(3);
              ar1_1->SetFillStyle(1001);
              ar1_1->Draw();      

              c1_1->Modified();
	          c1->cd();
         }
*/   
      /*******x title********/
      double x0=0.52, y0=0.04, x1=0.64, y1=0.08;
	  PLabel1X = new TLegend(x0,y0,x1,y1);
	  PLabel1X->SetTextAlign(22);
	  PLabel1X->SetTextFont(22);
	  PLabel1X->SetBorderSize(0);
	  PLabel1X->SetTextSize(0.05);
	  PLabel1X->SetTextColor(1);
	  PLabel1X->SetLineColor(0);
	  PLabel1X->SetFillColor(0);
	  PLabel1X->SetLineWidth(0);
	  PLabel1X->SetTextAngle(0);
      if(x==0)	  PLabel1X->SetHeader(_title1[i]);
      if(x==1)	  PLabel1X->SetHeader(_title2[i]);
	  PLabel1X->Draw();
      /**********y title***********/
      double yx0=0.015,yy0=0.45, yx1=0.055, yy1=0.6;
	  PLabel1Y = new TLegend(yx0,yy0,yx1,yy1);
	  PLabel1Y->SetTextAlign(22);
	  PLabel1Y->SetTextFont(22);
	  PLabel1Y->SetBorderSize(0);
	  PLabel1Y->SetTextSize(0.05);
	  PLabel1Y->SetTextColor(1);
	  PLabel1Y->SetLineColor(0);
	  PLabel1Y->SetFillColor(0);
	  PLabel1Y->SetLineWidth(0);
	  PLabel1Y->SetTextAngle(90);
	  PLabel1Y->SetHeader("S/#sqrt{S+B}");
	  PLabel1Y->Draw();

      /**********var<requiment***********/
      stringstream temx;
      temx<<_title1[i]<<_pm[i]<<_max1[i];
      TString temxx;
      temxx=temx.str();
      if(i<6){
      if(i!=4) {double zx0=0.27, zx1=0.47, zy0=0.8, zy1=0.9;}
      if(i==4) {double zx0=0.6, zx1=0.8, zy0=0.8, zy1=0.9;}
      if(i==5) {double zx0=0.47, zx1=0.7, zy0=0.8, zy1=0.9;}
	  PLabel1Z = new TLegend(zx0,zy0,zx1,zy1);
	  PLabel1Z->SetTextAlign(22);
	  PLabel1Z->SetTextFont(22);
	  PLabel1Z->SetBorderSize(0);
	  PLabel1Z->SetTextSize(0.05);
	  PLabel1Z->SetTextColor(1);
	  PLabel1Z->SetLineColor(0);
	  PLabel1Z->SetFillColor(0);
	  PLabel1Z->SetLineWidth(0);
	  PLabel1Z->SetTextAngle(0);
	  //PLabel1Z->SetHeader(temxx);
 	  if(i==0) {PLabel1Z->SetHeader("M_{K^{-}#pi^{+}#pi^{-}e(#pi)^{+}}<1.81 (GeV/#font[12]{c}^{2})");}
 	  if(i==1) {PLabel1Z->SetHeader("cos(e^{+},#pi^{-})<0.94");}
 	  if(i==2) {PLabel1Z->SetHeader("M_{K^{-}#pi^{-}#pi^{+}e(#pi)^{+}#pi^{0}}<1.4 (GeV/#font[12]{c}^{2}");}
 	  if(i==3) {PLabel1Z->SetHeader("cos(miss,#gamma)<0.81");}
 	  if(i==4) {PLabel1Z->SetHeader("M_{#pi^{+}#pi^{-}}<0.31 (GeV/#font[12]{c}^{2})");}
 	  if(i==5) {PLabel1Z->SetHeader("|#Delta E_{(K#pi)_{tag}#pi_{sig}}|<0.012 (GeV)");}
      PLabel1Z->Draw();
      }
       TLatex lx;
       lx.SetNDC();
       lx.SetTextAngle(0);
       lx.SetTextSize(0.05);
       if(i==0||i==3||i==5) lx.DrawLatex(0.8,0.8,"(a)");
       if(i==1) lx.DrawLatex(0.8,0.85,"(a)");
       if(i==4) lx.DrawLatex(0.8,0.75,"(a)");
       if(i==5) lx.DrawLatex(0.4,0.8,"(a)");
       





	  c1->Update();
	  c1->Modified();
	  c1->Print(Form("./figure/cut%d%d.eps",i,x));
	  delete ar1;
	  delete c1;
}

void Optimize::~Optimize(){
delete _sig;
delete _bkg;
delete _tree;
}

void Optimize::fun_opt1(int i,int bin_temp){
for(int j=0;j<=bin_temp;j++){

if(_pm[i]=="<"){
	value_x[i]=_cut_max1[i]-j*((_cut_max1[i]-_cut_min1[i])/bin_temp);
}
else {
	value_x[i]=_cut_min1[i]+j*((_cut_max1[i]-_cut_min1[i])/bin_temp);;
}

stringstream tem;
tem<<_cut1[i]<<_pm[i]<<value_x[i];
CUT1[i]=tem.str();


TString all;
for(int xx=0;xx<num-1;xx++)all+=CUT1[xx]+"&&";
all+=CUT1[num-1];

double factor = 1.416;
double nsig=_sig->GetEntries(all);
double nsig_revised = nsig*factor;
double nbkg=_bkg->GetEntries(all);
double sb=0;
if(nsig_revised+nbkg!=0)sb=nsig_revised/sqrt(nsig_revised+nbkg);
if(nsig_revised+nbkg==0)sb=0;
cut1[i][j]=value_x[i];
cut2[i][j]=0;
s[i][j][0]=sb;

//cout<<__LINE__<<"\t"<<all<<"\t"<<sb<<"\t"<<nsig<<"\t"<<nbkg<<endl;
//cout<<sb<<endl;
if(sb>=sb_cut){
	sb_cut=sb;
	sb_max=sb;
	if(_pm[i]=="<")value_max1[i]=_cut_max1[i]-j*((_cut_max1[i]-_cut_min1[i])/bin_temp);
	if(_pm[i]==">")value_max1[i]=_cut_min1[i]+j*((_cut_max1[i]-_cut_min1[i])/bin_temp);
	CUT_final=all;
	_each_max[i]=sb_max;
}
}
}

void Optimize::fun_opt2(int i,int bin_temp1,int bin_temp2){
for(int j=0;j<=bin_temp1;j++){
for(int k=0;k<=bin_temp2;k++){

if(_pm[i]=="<"){
	value_x[i]=_cut_max1[i]-j*((_cut_max1[i]-_cut_min1[i])/bin_temp1);
	value_y[i]=_cut_max2[i]-k*((_cut_max2[i]-_cut_min2[i])/bin_temp2);
}
else {
	value_x[i]=_cut_min1[i]+j*((_cut_max1[i]-_cut_min1[i])/bin_temp1);
	value_y[i]=_cut_min2[i]+k*((_cut_max2[i]-_cut_min2[i])/bin_temp2);
}

stringstream tem;
tem<<"("<<_cut1[i]<<_pm[i]<<value_x[i]<<"*"<<_cut2[i]<<"+("<<value_y[i]<<"))";
CUT1[i]=tem.str();


TString all;
for(int xx=0;xx<num-1;xx++)all+=CUT1[xx]+"&&";
all+=CUT1[num-1];

double factor = 1.416;
double nsig=_sig->GetEntries(all);
double nsig_revised = nsig*factor;
double nbkg=_bkg->GetEntries(all);
double sb=0;
if(nsig_revised+nbkg!=0)sb=nsig_revised/sqrt(nsig_revised+nbkg);
if(nsig_revised+nbkg==0)sb=0;
cut1[i][j]=value_x[i];
cut2[i][k]=value_y[i];
s[i][j][k]=sb;

//cout<<all<<"\t"<<sb<<"\t"<<nsig<<"\t"<<nbkg<<"\t"<<endl;
D2_sb[j][k]=sb;
D2_fig[j][0]=value_x[i];
D2_fig[0][k]=value_y[i];
//cout<<sb<<endl;
if(sb>=sb_cut){
	sb_cut=sb;
	sb_max=sb;
	if(_pm[i]=="<"){
		value_max1[i]=_cut_max1[i]-j*((_cut_max1[i]-_cut_min1[i])/bin_temp1);
		value_max2[i]=_cut_max2[i]-k*((_cut_max2[i]-_cut_min2[i])/bin_temp2);
	}
	if(_pm[i]==">"){
		value_max1[i]=_cut_min1[i]+j*((_cut_max1[i]-_cut_min1[i])/bin_temp1);
		value_max2[i]=_cut_min2[i]+k*((_cut_max2[i]-_cut_min2[i])/bin_temp2);
	}
	CUT_final=all;
	_each_max[i]=sb_max;
}
}
}
}


