#include <iomanip>
void ePID_plus()
{

   TCut cut_f  = "abs(isfound - 1.0) < 0.0001";     
   //TCut cut_f  = "isfound==1"; 
   TCut cut_uf = "abs(isfound - 0.0) < 0.0001";
   //TCut cut_uf = "isfound==0";
   //TCut cut_eop= "(eop_E-0.6)>1./10*dedxchi_E"; //0.8
   TCut cut_p[19]={"p>-0.1&&p<0.2","p>=0.2&&p<0.3","p>=0.3&&p<0.4","p>=0.4&&p<0.5","p>=0.5&&p<0.6","p>=0.6&&p<0.7","p>=0.7&&p<0.8","p>=0.8&&p<0.9","p>=0.9&&p<1.0","p>=1.0&&p<1.1","p>=1.1&&p<1.2","p>=1.2&&p<1.3","p>=1.3&&p<1.4","p>=1.4&&p<1.5","p>=1.5&&p<1.6","p>=1.6&&p<1.7","p>=1.7&&p<1.8","p>=1.8&&p<1.9","p>=1.9&&p<2.0"}; 
   TCut cut_cos[10]={"cos>=-0.93&&cos<-0.8","cos>=-0.8&&cos<-0.6","cos>=-0.6&&cos<-0.4","cos>=-0.4&&cos<-0.2","cos>=-0.2&&cos<0.0","cos>=0.0&&cos<0.2","cos>=0.2&&cos<0.4","cos>=0.4&&cos<0.6","cos>=0.6&&cos<0.8","cos>=0.8&&cos<0.93"};

   TChain mc("Bhabha");
   mc.Add("./dat/plus_mc_bhabha.root");
   TChain data("Bhabha");
   data.Add("./dat/plus_data_bhabha.root");
   ofstream out("./dat/plus_delta_eff1.txt"); 

   Double_t nf_data[19][10], nuf_data[19][10], nf_mc[19][10], nuf_mc[19][10], eff_data[19][10], eff_mc[19][10],delta_eff[19][10], eff_data_err[19][10], eff_mc_err[19][10], delta_eff_err[19][10], diff[19][10];
   for (int i=0;i<19;i++){ 
       for(int j=0;j<10;j++){
       //std::cout<<__LINE__<<"\t"<<i<<"\t"<<j<<std::endl;
       nf_data[i][j]   = data.GetEntries(cut_f&&cut_p[i]&&cut_cos[j]);
       nuf_data[i][j]  = data.GetEntries(cut_uf&&cut_p[i]&&cut_cos[j]);
       nf_mc[i][j]     = mc.GetEntries(cut_f&&cut_p[i]&&cut_cos[j]); 
       nuf_mc[i][j]    = mc.GetEntries(cut_uf&&cut_p[i]&&cut_cos[j]);
       eff_data[i][j]  = nf_data[i][j]/(nuf_data[i][j]+nf_data[i][j]);
       //std::cout<<__LINE__<<"\t"<<i<<"\t"<<j<<std::endl;
       eff_mc[i][j]    = nf_mc[i][j]/(nuf_mc[i][j]+nf_mc[i][j]);       
       //std::cout<<__LINE__<<"\t"<<i<<"\t"<<j<<std::endl;
       delta_eff[i][j] = eff_data[i][j]/eff_mc[i][j];   
       diff[i][j] = delta_eff[i][j]-1.;
       //err
       eff_data_err[i][j]  = err_single(nf_data[i][j],nuf_data[i][j]);
       eff_mc_err[i][j]    = err_single(nf_mc[i][j],nuf_mc[i][j]);       
   
       delta_eff_err[i][j] = err( eff_data[i][j],eff_data_err[i][j],eff_mc[i][j],eff_mc_err[i][j] );    
       
       out<<eff_data[i][j]<<"\t"<<eff_data_err[i][j]<<"\t"<<eff_mc[i][j]<<"\t"<<eff_mc_err[i][j]<<"\t"<<delta_eff[i][j]<<"\t"<<delta_eff_err[i][j]<<"\t"<<diff[i][j]<<"\n";
       }
   }
   
}

    double err_single(double e_a,double e_b)
    {
      Double_t d = sqrt(e_a*e_b/(e_a+e_b))/(e_a+e_b);
      return d;
    }

    double err(double m_1,double e_1,double m_2,double e_2)
   {
    double c=sqrt(e_1**2/m_2**2+m_1**2/m_2**4*e_2**2);
    return c;
   }

