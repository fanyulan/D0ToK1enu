#include<cstdio>
#include<string> 
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <iomanip>
    
#include "RooAbsReal.h"
#include "RooFit.h"
#include "RooRealVar.h"
#include "TMath.h"

using namespace RooFit ;
using std::vector; 

void correction(){
    double before[4];
    double err_before[4];
    double after[4];
    double err_after[4];
    ifstream eff_before("./dat/eff.txt");
    eff_before>>before[0]>>before[1]>>before[2]>>before[3];
    double K_pi_e[8] = {1.027, 1.000, 0.999, 0.999, 1.000, 0.999, 1.010, 0.986};//k track, k pid, pi1 track, pi1 pid, pi2 track, pi2 pid, e pid
    ifstream err("./dat/err.txt"); 
    err>>err_before[0]>>err_before[1]>>err_before[2]>>err_before[3];
    ofstream correction("./dat/correction.txt"); 
    ofstream corr("./dat/correction.dat"); 
    double sum=1.;  
    double eff_0[2], eff_1[2], eff_3[2], ave[2];
    for(int j=0;j<8;j++)
              {sum = sum*K_pi_e[j];} 
    for (int i=0; i<4; i++) {
              after[i] = before[i]*sum;
              correction<<"eff: "<<after[i]<<"\n";
              err_after[i] = err_before[i]*sum;
              correction<<"err: "<<err_after[i]<<"\n";
              corr<<after[i]<<"\n";
              
              std::cout<<sum<<"\t"<<after[i]<<"\t"<<err_after[i]<<std::endl;
}
}
