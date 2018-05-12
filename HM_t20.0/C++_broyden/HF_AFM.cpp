//c++ HF_AFM.cpp -std=c++1y ////this code giving the result that is matching with the self consistent calculation result where we calculated n_{\alpha,\sigma} self consistently
// this for AFM phase of HM with nearest neighbour hopping only
#include<iostream>
//#include <array>
#include <math.h>       /* cos */
#define PI 3.14159265
#include <vector>
#include <cstdlib>     /* abs */
#include<fstream>// header file for input and output
#include <sstream> 
#include"roots_multidim.h"

using namespace std;
class HF {
  long double t,t2,delta,A,U;
  const static int dim = 2,GRID = 500, bin = 500;
  long double cos_[2*GRID] = {};
  long double lambda_p[2*GRID][2*GRID] = {{}};
  long double lambda_n[2*GRID][2*GRID] = {{}};
  long double lambda_upp[2*GRID][2*GRID] = {{}};
  long double lambda_dnp[2*GRID][2*GRID] = {{}};
  long double lambda_upn[2*GRID][2*GRID] = {{}};
  long double lambda_dnn[2*GRID][2*GRID] = {{}};
  long double L_upp[2*GRID][2*GRID] = {{}};
  long double L_dnp[2*GRID][2*GRID] = {{}};
  long double L_upn[2*GRID][2*GRID] = {{}};
  long double L_dnn[2*GRID][2*GRID] = {{}};
  long double CA_upp[2*GRID][2*GRID] = {{}};
  long double CA_dnp[2*GRID][2*GRID] = {{}};
  long double CA_upn[2*GRID][2*GRID] = {{}};
  long double CA_dnn[2*GRID][2*GRID] = {{}};
  long double CB_upp[2*GRID][2*GRID] = {{}};
  long double CB_dnp[2*GRID][2*GRID] = {{}};
  long double CB_upn[2*GRID][2*GRID] = {{}};
  long double CB_dnn[2*GRID][2*GRID] = {{}};
  long double dos_upp[bin] = {};
  long double dos_upn[bin] = {};
  long double dos_dnp[bin] = {};
  long double dos_dnn[bin] = {};

  public:
  long double nA_upn,nA_dnn,nB_upn,nB_dnn;
  int cnt ;
  HF(long double t_, long double t2_, long double delta_, int dim_){
  	t = t_;
	t2 = t2_;
	delta = delta_;
  
  for(int i = 0; i<2*GRID; i++){
	cos_[i] = cos( (i-GRID)*PI/(GRID));
	}
	} // HF

  long double division(long double a, long double b)
	{
	   if( b == 0 )
	   {
	      throw "Division by zero condition!";
	   }
	   return (a/b);
	}  

  VecDoub solve(VecDoub_I x){
  auto ms = x[0];
  auto mu = U/2.0;
  //auto ms = x[1];
  auto dn = 0.0;
  Doub total_no_part = 0.0;
  Doub free_energy_con = 0.0;
  auto n = x.size();
  VecDoub fvec(n);
  int N_p, N_n;
  auto nA_up_ =  new (nothrow) long double[2*GRID]();
  auto nB_up_ =  new (nothrow) long double[2*GRID]();
  auto nA_dn_ =  new (nothrow) long double[2*GRID]();
  auto nB_dn_ =  new (nothrow) long double[2*GRID]();
  cnt = 0;
  auto nA_up =  0.5*(1.0-dn+ms);
  auto nB_up =  0.5*(1.0+dn-ms);
  auto nA_dn =  0.5*(1.0-dn-ms);
  auto nB_dn =  0.5*(1.0+dn+ms);
  N_p = 0;
  N_n = 0;
  auto weight = 0.0;
  // cout<<"input"<<nA_up<<nA_dn<<nB_up<<nB_dn<<endl;
   for(int i = 0; i<GRID; i++){
	for(int j = 0; j<(i+1); j++){
		if (j!=i){weight = 8.0;}
		else{weight = 4.0;}
		auto energy       =  -2*t*(cos_[i] + cos_[j]) ; 
		auto energy_prime = -2*dim*t2*cos_[i]*cos_[j]; 
		auto g = U*ms/2.0;
		auto C = sqrt( energy*energy + g*g );
		lambda_p[i][j] =  energy_prime + C ;
		lambda_n[i][j] =  energy_prime - C ;
		if (energy!=0.0){
		if (lambda_p[i][j]<0.0){
			N_p = N_p + weight*1.0;
			//total_no_part = total_no_part + 1.0;
			 try {
	     		 auto r  = division(1,C);
	   		}catch (const char* msg) {
		                 cout<<"C coming zer0:"<<C<<endl;
				}
			free_energy_con  = free_energy_con - weight/C;
			}
		if (lambda_n[i][j]<0.0){ 
			N_n = N_n + weight;
			total_no_part = total_no_part + weight;
			free_energy_con  = free_energy_con + weight/C;
			}
		}//energy!=0.0
		
	} // for(int i = 0; i<2*GRID; i++)
  } // for(int j = 0; j<2*GRID; j++)

 
 

 //fvec[0] = total_no_part/(4*GRID*GRID) - 1.0;
 fvec[0] = (free_energy_con*U)/(4*GRID*GRID) -2.0;
 cout<<"N_p,N_n:"<<N_p/(4*GRID*GRID)<<" "<<N_n/(4*GRID*GRID)<<endl;
 return fvec;
 }	//solve

  long double get_GRID(){ return 2*t; };
  auto get_lambda(){return lambda_upp[0][0];}
  void set_U(double U_){U=U_;}


};//class

long double abs(long double a){
	if (a>0.0){
		return a;
	}
	else{
		return -1.0*a;
	}
	}


 VecDoub vecfunc(VecDoub_I x){
  long double t,t2,delta,U,mu,conv,ms,dn;
  int dim,GRID,loop;
  
  VecDoub fnvec;
  t = 0.5;
  t2 = 0.0;
  delta = 0.0;
  dim = 2;
  //GRID = 200;
  
  ifstream infile; 
  infile.open("U_file.in"); 
  //if(infile.is_open()){std::cout<<"file occupation_in.dat is open"<<endl;}  
  std::string line;
  std::getline(infile, line);
  std::istringstream ss(line);
  ss >> U ;
  
  HF hf(t,t2,delta,dim);
 
  hf.set_U(U);
  fnvec = hf.solve(x);
  cout << "U:"<<U<<"  "<<fnvec[0]<<endl;
return fnvec;
}
int main(){
  bool check = false;
  VecDoub_IO  x(1);
  
  x[0] = 0.1;
  //vector<long double> U_list ={1.6,1.5,1.4,1.3,1.2,1.1,1.0,0.9,0.8,0.7};
  vector<long double> U_list ={0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.5,3.0,3.5,4.0};
  ofstream soumen("broyden_AFM_t20.dat");
  soumen <<"# U,mu,ms,ntotal-2,ms-ms1"<<endl;
  soumen <<"# only positive solution was kept"<<endl;
  for(long double U : U_list) {
	  check = false;
          //x[0] = U/2.0;
	  ofstream U_file("U_file.in");	  
	  U_file <<U<<endl;
	  U_file.close();

  	broydn(x, check, *vecfunc );
 	
  	auto m = vecfunc( x);
  	auto mu = U/2.0;
        auto ms = x[0];
	
  	// write n values
         cout <<"# U,mu,ms,ms-ms1"<<endl;
  	 cout  << U<<"  "<<mu<<"   "<<ms<<"	"<<m[0]<<endl;
  	 soumen << U<<"  "<<mu<<"   "<<ms<<"	"<<m[0]<<endl;
	//if(x[1]<0.01) x[1] = 0.1;
	//if(x[0]<0.01) x[0] = 0.1;
  }
  soumen.close();
  

  return 0;

} //main
//TODO does it have dispersion for upspin and downspin.
