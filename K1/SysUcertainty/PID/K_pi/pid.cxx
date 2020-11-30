#include <iomanip>
void pid()
{
      double   pid_mc[11] = {0},   pid_mc_err[11] = {0},
               pid_data[11]={0},   pid_data_err[11]={0},
               pid_corr[11]={0},   pid_corr_err[11]={0};
 
       ifstream mc("./dat/pi_pid_mc.dat");
       ifstream data("./dat/pi_pid_data.dat");
       //ifstream mc("./dat/k_pid_mc.dat");
       //ifstream data("./dat/k_pid_data.dat");

       for (int i=0;i<11;i++){
          mc>>pid_mc[i]>>pid_mc_err[i];
          data>>pid_data[i]>>pid_data_err[i];

      }

////////////////////////////////////////////////////truth weight
      //pi2
      TString P_cut[11]={"p_pion2p3<0.1","p_pion2p3>=0.1&&p_pion2p3<0.2","p_pion2p3>=0.2&&p_pion2p3<0.3","p_pion2p3>=0.3&&p_pion2p3<0.4","p_pion2p3>=0.4&&p_pion2p3<0.5","p_pion2p3>=0.5&&p_pion2p3<0.6","p_pion2p3>=0.6&&p_pion2p3<0.7","p_pion2p3>=0.7&&p_pion2p3<0.8","p_pion2p3>=0.8&&p_pion2p3<0.9","p_pion2p3>=0.9&&p_pion2p3<1.0","p_pion2p3>=1.0"};
      //pi1
      //TString P_cut[11]={"p_pionp3<0.1","p_pionp3>=0.1&&p_pionp3<0.2","p_pionp3>=0.2&&p_pionp3<0.3","p_pionp3>=0.3&&p_pionp3<0.4","p_pionp3>=0.4&&p_pionp3<0.5","p_pionp3>=0.5&&p_pionp3<0.6","p_pionp3>=0.6&&p_pionp3<0.7","p_pionp3>=0.7&&p_pionp3<0.8","p_pionp3>=0.8&&p_pionp3<0.9","p_pionp3>=0.9&&p_pionp3<1.0","p_pionp3>=1.0"};
     //k
     //TString P_cut[11]={"p_kaonp3<0.1","p_kaonp3>=0.1&&p_kaonp3<0.2","p_kaonp3>=0.2&&p_kaonp3<0.3","p_kaonp3>=0.3&&p_kaonp3<0.4","p_kaonp3>=0.4&&p_kaonp3<0.5","p_kaonp3>=0.5&&p_kaonp3<0.6","p_kaonp3>=0.6&&p_kaonp3<0.7","p_kaonp3>=0.7&&p_kaonp3<0.8","p_kaonp3>=0.8&&p_kaonp3<0.9","p_kaonp3>=0.9&&p_kaonp3<1.0","p_kaonp3>=1.0"};
      


      TChain tsig("tagD");
      tsig.Add("tagD_pid.root");
      int n1=tsig.GetEntries();
      double n[11]={0},w[11]={0};
      double pid=0;
      double pid_sum=0;
      double pid_sum_err=0;
      double pid_sum_err1=0;
          for (int i=0;i<11;i++){

               n[i]=tsig.GetEntries(P_cut[i]);
               w[i]=n[i]/n1;
               pid = (pid_data[i]/pid_mc[i])*w[i];
               pid_sum += pid;               
               
               double pid_err = err(pid_data[i],pid_data_err[i],pid_mc[i],pid_mc_err[i])*w[i]; 
               pid_sum_err1 += (pid_err[i])**2; 
               }

               pid_sum_err = sqrt(pid_sum_err1);
cout<<"pid pi2:"<<pid_sum<<"\t"<<"err:"<<pid_sum_err<<endl;
//cout<<"pid pi1:"<<pid_sum<<"\t"<<"err:"<<pid_sum_err<<endl;
//cout<<"pid k:"<<pid_sum<<"\t"<<"err:"<<pid_sum_err<<endl;

}

double err(double m_1,double e_1,double m_2,double e_2)
{
   double c=sqrt(e_1**2/m_2**2+m_1**2/m_2**4*e_2**2);
   return c;
}


