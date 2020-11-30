void draw(){
    	using namespace RooFit;
    	gSystem->Load("libRooFit");
        gROOT->Reset();
        gROOT->SetStyle("Plain");
        gStyle->SetOptTitle(0);
        gStyle->SetOptStat(0);
        gStyle->SetOptFit(0);
        gStyle->SetLabelSize(0.07,"xyz");
        gStyle->SetNdivisions(507,"xyz");        
        gROOT->SetStyle("Plain");
        gStyle->SetLabelSize(0.06,"xyz");
        gStyle->SetNdivisions(405,"xyz");
        gStyle->SetPadTopMargin(.10);
        gStyle->SetPadLeftMargin(.10);
        gStyle->SetPadRightMargin(.05);
        gStyle->SetPadBottomMargin(.10);
        gStyle->SetTitleSize(0.06,"xyz");
        gStyle->SetOptTitle(0);
        gStyle->SetMarkerSize(0.5);
        gStyle->SetTitle("");
        gStyle->SetOptStat(0);
        gStyle->SetEndErrorSize(0);


        ifstream BF("./figure/bf_new1.txt");
    	
        TCanvas *c = new TCanvas("c", "c",0,0,800,600);
	    leg = new TLegend(0.0,0.0,0.98,0.05);
        leg->SetTextAlign(22);
        leg->SetTextAlign(22);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.025);
        leg->SetTextColor(1);
        leg->SetLineColor(0);
        leg->SetFillColor(0);
        leg->SetLineWidth(0);
        leg->SetTextAngle(0);
        leg->SetHeader("Sample");
        leg->Draw();
        ///////
        TF1 *f_input = new TF1("f_input","FInput(x,y)", 0.0, 11.0);                                                                               
        f_input->SetLineColor(2);
        //////
        double bf_input[10]={0};
        bf_input[0]=1;
        bf_input[1]=2;
        bf_input[2]=3;
        bf_input[3]=4;
        bf_input[4]=5;
        bf_input[5]=6;
        bf_input[6]=7;
        bf_input[7]=8;
        bf_input[8]=9; 
        bf_input[9]=10;

	    TGraphErrors *g1=new TGraphErrors();
        double sum = 0, average = 0.;
    	double x1, y1, z1;
        double k1 = 0.282448; 
	    int i=0;
     	while(BF>>y1>>z1){
        g1->SetPoint(i, bf_input[i], y1*0.32852/k1);
        g1->SetPointError(i, 0.0, z1*0.32852/k1);
        cout<<z1<<endl;
        sum+=y1; 
        i++;
        }
        average = sum/10;
        //std::cout<<average<<std::endl;

        g1->GetXaxis()->SetRangeUser(0, 11);
        g1->GetYaxis()->SetRangeUser(0, 1.6e-03);
        g1->GetXaxis()->SetNdivisions(505);
        g1->GetYaxis()->SetNdivisions(505);
        g1->GetXaxis()->SetLabelSize(0.04);
        g1->GetYaxis()->SetLabelSize(0.04);

	    g1->GetXaxis()->SetTitleSize(0.02);
        g1->GetYaxis()->SetTitle(" Out BF");
        g1->GetYaxis()->SetTitleSize(0.05);
        g1->GetXaxis()->CenterTitle();
        g1->GetYaxis()->CenterTitle();
        g1->GetYaxis()->SetTitleOffset(0.9);
        g1->GetXaxis()->SetTitleOffset(1.0);
        g1->SetMarkerStyle(20);
        g1->SetMarkerSize(0.9);
        g1->SetMarkerColor(kBlue);
        g1->Draw("AP");
        g1->Fit("pol0");

        //box_1 = new TBox(0.1, 0.00076, 10.9, 0.00076);
        //box_1->SetFillColor(kRed);
        //box_1->Draw("same");
        /*	TLatex lt;
        lt.SetNDC();
        lt.SetTextAngle(0);
        lt.SetTextSize(0.06);
        lt.DrawLatex(0.2, 0.8,  Form("D^{0}#rightarrow K^{-}#pi^{+}#eta"));
        */
        TLatex lt1;
        lt1.SetNDC();
        lt1.SetTextAngle(0);
        lt1.SetTextSize(0.045);
     
        lt1.DrawLatex(0.53, 0.2,   Form("#color[2]{Input:  7.60#times10^{-4}}"));
        lt1.DrawLatex(0.53, 0.12,  Form("#color[4]{Output:(7.33#pm0.35)#times10^{-4}}"));
//        lt1.DrawLatex(0.53, 0.2,   Form("#color[2]{Difference: %%}"));

        f_input->Draw("same");       
        c->cd();
        c->Print("./figure/draw.eps");
}

    double FInput(double x, double y)
    {
         if(x>0&&x<11) double y = 0.00076;//0.00076*0.282448
         return y;
    }
