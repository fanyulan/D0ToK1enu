#include <iomanip>
void err()
{
      double   mc[6][11] = {0},   data[6][11]={0}, 
               eff[6][11]={0},    eff_err[6][11]={0}, 
               weight[6][11]={0}, w_err[6][11]={0};
      //input 
      ifstream fmc("./dat/mc_eff.log");
      ifstream fdata("./dat/data_eff.log");
      ifstream fratio("./dat/ratio_eff.log");
      ifstream fratio_err("./dat/ratio_eff_err.log");
      ifstream fweight("./dat/weight.log");
      ifstream fweight_err("./dat/weight_err.log");

      double sum=0, sum_err=0;
      for (int i=0;i<6;i++){
           for(int j=0;j<10;j++){
           fmc>>mc[i][j];       
           fdata>>data[i][j];
           fratio>>eff[i][j];
           fratio_err>>eff_err[i][j];
           fweight>>weight[i][j];
           fweight_err>>w_err[i][j];
   
           sum += eff[i][j]*weight[i][j];
           sum_err += (weight[i][j]*eff_err[i][j])**2+(2*w_err[i][j]*eff[i][j])**2;

         
          
          } 
      }        
     cout<<" err:"<<sqrt(sum_err)<<endl;
}


