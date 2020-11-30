#include <iomanip>
void trk()
{
      double   trk_mc[11] = {0},   trk_mc_err[11] = {0},
               trk_data[11]={0},   trk_data_err[11]={0},
               trk_corr[11]={0},   trk_corr_err[11]={0};
      //input 
      //ifstream mc("../dat/k_trk_mc.dat");
      //ifstream data("../dat/k_trk_data.dat");
       
      //ifstream mc("../dat/pip_trk_mc.dat");
      //ifstream data("../dat/pip_trk_data.dat");

      ifstream mc("../dat/pin_trk_mc.dat");
      ifstream data("../dat/pin_trk_data.dat");

      //K
      //TString P_cut[11]={"ptkK>0&&ptkK<0.1","ptkK>=0.1&&ptkK<0.2","ptkK>=0.2&&ptkK<0.3","ptkK>=0.3&&ptkK<0.4","ptkK>=0.4&&ptkK<0.5","ptkK>=0.5&&ptkK<0.6","ptkK>=0.6&&ptkK<0.7","ptkK>=0.7&&ptkK<0.8","ptkK>=0.8&&ptkK<0.9","ptkK>=0.9&&ptkK<1.0","ptkK>=1.0"};
      //pi1
      //TString P_cut[11]={"ptkPi1<0.1","ptkPi1>=0.1&&ptkPi1<0.2","ptkPi1>=0.2&&ptkPi1<0.3","ptkPi1>=0.3&&ptkPi1<0.4","ptkPi1>=0.4&&ptkPi1<0.5","ptkPi1>=0.5&&ptkPi1<0.6","ptkPi1>=0.6&&ptkPi1<0.7","ptkPi1>=0.7&&ptkPi1<0.8","ptkPi1>=0.8&&ptkPi1<0.9","ptkPi1>=0.9&&ptkPi1<1.0","ptkPi1>=1.0"};
      //pi2
      TString P_cut[11]={"ptkPi<0.1","ptkPi>=0.1&&ptkPi<0.2","ptkPi>=0.2&&ptkPi<0.3","ptkPi>=0.3&&ptkPi<0.4","ptkPi>=0.4&&ptkPi<0.5","ptkPi>=0.5&&ptkPi<0.6","ptkPi>=0.6&&ptkPi<0.7","ptkPi>=0.7&&ptkPi<0.8","ptkPi>=0.8&&ptkPi<0.9","ptkPi>=0.9&&ptkPi<1.0","ptkPi>=1.0"};

      TChain tsig("truth");
      tsig.Add("../truth.root");
      int n1=tsig.GetEntries();
      double n[11]={0},w[11]={0},w_err[11]={0};
      double trk_mc_sum=0,trk_data_sum=0;
      for (int i=0;i<11;i++){
    
          mc>>trk_mc[i]>>trk_mc_err[i];       
          data>>trk_data[i]>>trk_data_err[i];
          
          n[i]=tsig.GetEntries(P_cut[i]);
          w[i]=n[i]/n1;      //weight 
          if(w[i]==0)  continue;
          w_err[i]=sqrt(w[i]*(1-w[i])/n1);      //err of weight     
          
          trk_mc_sum+=trk_mc[i]*w[i];           // mc*w
          trk_data_sum+=trk_data[i]*w[i];       //data*w
          // wrong err
          //trk_data_err_sum += err_new(trk_data[i],trk_data_err[i],w[i],w_err[i]); // (data+err_data)*(w+w_err)        
          //trk_mc_err_sum += err_new(trk_mc[i],trk_mc_err[i],w[i],w_err[i]);    
    }

      //cout<<"track K:"<<trk_data_sum/trk_mc_sum<<endl;
      //cout<<"track pi1:"<<trk_data_sum/trk_mc_sum<<endl;
      cout<<"track pi2:"<<trk_data_sum/trk_mc_sum<<endl;
    
}

double err(double m_1,double e_1,double m_2,double e_2)
{
      double c=sqrt(e_1**2/m_2**2+m_1**2/m_2**4*e_2**2);//error transform fomular
      return c;
}

double err_new(double a1, double e1, double a2, double e2)
{
    double d = sqrt(a2**2*(e1)**2+a1**2*(e2*0.4)**2);
    return d;
}

