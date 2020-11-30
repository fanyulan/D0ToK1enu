#include <iomanip>
void e_trk()
{
      double   trk_mc[6][11] = {0},  trk_mc_err[6][11] = {0},
               trk_data[6][11]={0},  trk_data_err[6][11]={0},   
               trk_weight[6][11]={0}, w_err[6][11]={0};
      //input 
      ifstream mc("./eff_mc.log");
      ifstream data("./eff_data.log");
      ifstream mc_err("./err_mc.log");
      ifstream data_err("./err_data.log");
      ifstream weight("./weight.log");
      ifstream weight_err("./weight_err.log");

      //read
      double trk_mc_sum = 0;
      double trk_data_sum = 0;
      double trk_mc_sum_err = 0.,trk_mc_sum_err1=0.,trk_mc_err_sum=0. ;
      double trk_data_sum_err = 0., trk_data_sum_err1=0.,trk_data_err_sum=0.;
      for (int i=0;i<6;i++){
           for(int j=0;j<10;j++){
           mc>>trk_mc[i][j];       
           data>>trk_data[i][j];
           mc_err>>trk_mc_err[i][j];
           data_err>>trk_data_err[i][j];
           weight>>trk_weight[i][j];
           weight_err>>w_err[i][j];
   
           if(trk_weight[i][j]==0) continue;
          trk_mc_sum += trk_mc[i][j]*trk_weight[i][j];
          trk_data_sum += trk_data[i][j]*trk_weight[i][j];
         
          trk_data_err_sum += err_new(trk_data[i][j],trk_data_err[i][j],trk_weight[i][j],w_err[i][j]);
          trk_mc_err_sum += err_new(trk_mc[i][j],trk_mc_err[i][j],trk_weight[i][j],w_err[i][j]);
          
          } 
      }        
         cout<<"ratio: "<<trk_data_sum/trk_mc_sum<<endl;
         //cout<<"ratio: "<<trk_data_sum/trk_mc_sum<<"\t"<<"err:"<<err(trk_data_sum,trk_data_sum_err1,trk_mc_sum,trk_mc_sum_err1)<<"\t"<<"new err:"<<err(trk_data_sum,trk_data_err_sum, trk_mc_sum,trk_mc_err_sum)<<endl;
}   
//x/y
double err(double m_1,double e_1,double m_2,double e_2)
{
      double c=sqrt(e_1*e_1*m_2*m_2+m_1*m_1*e_2*e_2)/(m_2*m_2);//error transform fomular
      return c;
}

//x*y
double err_new(double a1, double e1, double a2, double e2)
{
   double d = sqrt(a2**2*(e1*0.25)**2+a1**2*(e2*0.25)**2);
   return d;
}


