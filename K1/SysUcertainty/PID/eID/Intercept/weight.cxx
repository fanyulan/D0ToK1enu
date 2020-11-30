#include <iomanip>

void weight(){
  
    Double_t pid_differ[19][10],   pid_differ_err[19][10];
    Double_t weight[19][10];
    Double_t temp1[19][10], temp2[19][10], temp3[19][10], temp4[19][10], temp5[19][10]; 
    //ifstream differ("./dat/minus_delta_eff_1.txt");
    ifstream differ("./dat/plus_delta_eff1.txt");
    
    TCut cut_cos[10] = { "cosE<-0.8","cosE>=-0.8&&cosE<-0.6","cosE>=-0.6&&cosE<-0.4","cosE>=-0.4&&cosE<-0.2","cosE>=-0.2&&cosE<0.0","cosE>=0.0&&cosE<0.2","cosE>=0.2&&cosE<0.4","cosE>=0.4&&cosE<0.6","cosE>=0.6&&cosE<0.8","cosE>=0.8" };
    TCut cut_p[19] = {"elecP>=0.1&&elecP<0.2","elecP>=0.2&&elecP<0.3","elecP>=0.3&&elecP<0.4","elecP>=0.4&&elecP<0.5","elecP>=0.5&&elecP<0.6","elecP>=0.6&&elecP<0.7","elecP>=0.7&&elecP<0.8","elecP>=0.8&&elecP<0.9","elecP>=0.9&&elecP<1.0","elecP>=1.0&&elecP<1.1","elecP>=1.1&&elecP<1.2","elecP>=1.2&&elecP<1.3","elecP>=1.3&&elecP<1.4","elecP>=1.4&&elecP<1.5","elecP>=1.5&&elecP<1.6","elecP>=1.6&&elecP<1.7","elecP>=1.7&&elecP<1.8","elecP>=1.8&&elecP<1.9","elecP>=1.9&&elecP<2.0"};

    TChain tsig("tagD");
    tsig.Add("../../tagD_pid.root");
    int n1=tsig.GetEntries();   
    Double_t n[19][10];
    Double_t pid=0, pid_sum=0, pid_err=0, pid_err_sum=0; 
    
    ofstream OUT("figure/res.dat", ios::app);
    for(int i=0;i<19;i++){
       for(int j=0;j<10;j++){

        differ >> temp1[i][j] >> temp2[i][j] >> temp3[i][j] >> temp4[i][j] >> pid_differ[i][j] >> pid_differ_err[i][j] >> temp5[i][j];
//        differ >> pid_differ[i][j] >> pid_differ_err[i][j];
        n[i][j] = tsig.GetEntries(cut_cos[j]&&cut_p[i]);
        weight[i][j] = n[i][j]/n1; 

        pid = pid_differ[i][j]*weight[i][j];
        pid_sum += pid;

        pid_err = (pid_differ_err[i][j] * weight[i][j])**2;
        pid_err_sum += pid_err;

        cout<<i<<"\t"<<j<<" weight "<<weight[i][j]<<", pid "<<pid<<", err "<<pid_err<<endl;
      }
   }

    double err = sqrt(pid_err_sum); 
    cout <<"diff  "<< pid_sum << "\t" << err << endl; 
    //OUT<<"minus: "<<pid_sum<<"\t"<<err;
    OUT<<"plus: "<<pid_sum<<"\t"<<err;
}
