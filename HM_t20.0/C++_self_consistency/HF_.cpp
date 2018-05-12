//c++ HF_.cpp -std=c++1y
#include<iostream>
//#include <array>
#include <math.h>       /* cos */
#define PI 3.14159265
#include <vector>
#include <cstdlib>     /* abs */
#include<fstream>// header file for input and output
#include <sstream> 

using namespace std;
class HF {
  long double t,t2,delta,A;
  const static int dim = 2,GRID = 300, bin = 500;
  long double cos_[2*GRID] = {};
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

  void solve(long double U, long double mu, long double nA_up, long double nA_dn, long double nB_up, long double nB_dn){
  auto nA_up_ =  new (nothrow) long double[2*GRID]();
  auto nB_up_ =  new (nothrow) long double[2*GRID]();
  auto nA_dn_ =  new (nothrow) long double[2*GRID]();
  auto nB_dn_ =  new (nothrow) long double[2*GRID]();
  cnt = 0;
   //cout<<"solver class is created-------------------"<<endl;
   for(int i = 0; i<2*GRID; i++){
	for(int j = 0; j<2*GRID; j++){
		auto energy = 2*t*(cos_[i] + cos_[j]);
		auto energyA = 2*dim*t2*cos_[i]*cos_[j];
		//cout<<energy<<endl;
		auto A_up = delta -mu + U*nA_up - energyA;
		auto A_dn = delta -mu + U*nA_dn - energyA;
		auto B_up = -delta -mu + U*nB_up - energyA;
		auto B_dn = -delta -mu + U*nB_dn - energyA;
		lambda_upp[i][j] = 0.5*( A_dn + B_dn + sqrt( (A_dn - B_dn)*(A_dn - B_dn) + 4*energy*energy) );
		lambda_dnp[i][j] = 0.5*( A_up + B_up + sqrt( (A_up - B_up)*(A_up - B_up) + 4*energy*energy) );
		lambda_upn[i][j] = 0.5*( A_dn + B_dn - sqrt( (A_dn - B_dn)*(A_dn - B_dn) + 4*energy*energy) );
		lambda_dnn[i][j] = 0.5*( A_up + B_up - sqrt( (A_up - B_up)*(A_up - B_up) + 4*energy*energy) );	
		L_upp[i][j] = (A_up - lambda_dnp[i][j]);
                L_dnp[i][j] = (A_dn - lambda_upp[i][j]);
		L_upn[i][j] = (A_up - lambda_dnn[i][j]);
		L_dnn[i][j] = (A_dn - lambda_upn[i][j]);
		if(energy != 0.0 ){
		        try {
	     			CA_upp[i][j] = division(energy*energy,energy*energy + L_upp[i][j]*L_upp[i][j]);
	   		}catch (const char* msg) {
		                 cnt = cnt + 1;
	     			//cout<<"L_upp"<<A_up<< "  "<<B_up<<"  "<<i<<"   "<<j<<":energy:   "<<energy<<"  " <<energyA<<"   "<<(A_up - lambda_dnp[i][j])<<"  "<<cnt<< endl;
	  		}
		try {
	     			CA_upn[i][j] = division(energy*energy,energy*energy + L_upn[i][j]*L_upn[i][j]);
	   		}catch (const char* msg) {
		                 cnt = cnt + 1;
	     		//	cout<<"L_upp"<<A_up<< "  "<<B_up<<"  "<<i<<"   "<<j<<"   "<<energy<<"  " <<energyA<<"   "<<(A_up - lambda_dnp[i][j])<< endl;
			}
		try {
	     			CA_dnp[i][j] = division(energy*energy,energy*energy + L_dnp[i][j]*L_dnp[i][j]);
	   		}catch (const char* msg) {
		                 cnt = cnt + 1;
	     		//	cout<<"L_upp"<<A_up<< "  "<<B_up<<"  "<<i<<"   "<<j<<"   "<<energy<<"  " <<energyA<<"   "<<(A_up - lambda_dnp[i][j])<< endl;
			}
		try {
	     			CA_dnn[i][j] = division(energy*energy,energy*energy + L_dnn[i][j]*L_dnn[i][j]);
	   		}catch (const char* msg) {
		                 cnt = cnt + 1;
	     		//	cout<<"L_upp"<<A_up<< "  "<<B_up<<"  "<<i<<"   "<<j<<"   "<<energy<<"  " <<energyA<<"   "<<(A_up - lambda_dnp[i][j])<< endl;
			}

		// now B sublattice calculation start here.
		 try {
     				CB_upp[i][j] = division(L_upp[i][j]*L_upp[i][j],energy*energy + L_upp[i][j]*L_upp[i][j]);	
   		}catch (const char* msg) {
			cnt = cnt + 1;
			//cout<<"L_dnp"<<A_dn<< "  "<<B_dn<<"  "<<i<<"   "<<j<<"   "<<energy<<"  " <<energyA<<"   "<<(A_dn - lambda_upp[i][j])<< endl;
  		}
                 try {
     				CB_dnp[i][j] = division(L_dnp[i][j]*L_dnp[i][j],energy*energy + L_dnp[i][j]*L_dnp[i][j]);	
   		}catch (const char* msg) {
			cnt = cnt + 1;
			//cout<<"L_dnp"<<A_dn<< "  "<<B_dn<<"  "<<i<<"   "<<j<<"   "<<energy<<"  " <<energyA<<"   "<<(A_dn - lambda_upp[i][j])<< endl;
  		}
                try {
     				CB_upn[i][j] = division(L_upn[i][j]*L_upn[i][j],energy*energy + L_upn[i][j]*L_upn[i][j]);	
   		}catch (const char* msg) {
			cnt = cnt + 1;
			//cout<<"L_dnp"<<A_dn<< "  "<<B_dn<<"  "<<i<<"   "<<j<<"   "<<energy<<"  " <<energyA<<"   "<<(A_dn - lambda_upp[i][j])<< endl;
  		}
                 try {
     				CB_dnn[i][j] = division(L_dnn[i][j]*L_dnn[i][j],energy*energy + L_dnn[i][j]*L_dnn[i][j]);	
   		}catch (const char* msg) {
			cnt = cnt + 1;
			//cout<<"L_dnp"<<A_dn<< "  "<<B_dn<<"  "<<i<<"   "<<j<<"   "<<energy<<"  " <<energyA<<"   "<<(A_dn - lambda_upp[i][j])<< endl;
  		}
		}//end of if(energyA != 0.0)	

		//CA_upp[i][j] = (energy*energy)/(energy*energy + L_upp[i][j]*L_upp[i][j] );
		//CA_dnp[i][j] = (energy*energy)/( energy*energy + L_dnp[i][j]*L_dnp[i][j]);
		//CB_upp[i][j] = (L_upp[i][j]*L_upp[i][j])/( energy*energy + L_upp[i][j]*L_upp[i][j]);
		//CB_dnp[i][j] = (L_dnp[i][j]*L_dnp[i][j])/( energy*energy + L_dnp[i][j]*L_dnp[i][j]);
		//CA_upn[i][j] = (energy*energy)/(energy*energy + L_upn[i][j]*L_upn[i][j] );
		//CA_dnn[i][j] = (energy*energy)/( energy*energy + L_dnn[i][j]*L_dnn[i][j]);
		//CB_upn[i][j] = (L_upn[i][j]*L_upn[i][j])/( energy*energy + L_upn[i][j]*L_upn[i][j]);
		//CB_dnn[i][j] = (L_dnn[i][j]*L_dnn[i][j])/( energy*energy + L_dnn[i][j]*L_dnn[i][j]);
	
	} // for(int i = 0; i<2*GRID; i++)
  } // for(int j = 0; j<2*GRID; j++)

 
  nA_upn = 0.0;
  nA_dnn = 0.0;
  nB_upn = 0.0;
  nB_dnn = 0.0;
  for(int j = 0; j<2*GRID; j++){
	for(int i = 0; i<2*GRID; i++){
		
			if (lambda_upp[i][j]<0.0){
				nA_up_[j] = nA_up_[j] + CA_dnp[i][j]; 
			        nB_up_[j] = nB_up_[j] + CB_dnp[i][j]; 
				}
			if (lambda_dnp[i][j]<0.0){ 
				nA_dn_[j] = nA_dn_[j] + CA_upp[i][j];
				nB_dn_[j] = nB_dn_[j] + CB_upp[i][j];
				}
			if (lambda_upn[i][j]<0.0){ 
				nA_up_[j] = nA_up_[j] + CA_dnn[i][j]; 
				nB_up_[j] = nB_up_[j] + CB_dnn[i][j];
				}
			if (lambda_dnn[i][j]<0.0){
				nA_dn_[j] = nA_dn_[j] + CA_upn[i][j];
				nB_dn_[j] = nB_dn_[j] + CB_upn[i][j];
				}
			
	
	} //for(int i = 0; i<2*GRID; i++)

    
		nA_upn = nA_upn + nA_up_[j];
		nA_dnn = nA_dnn + nA_dn_[j];
		nB_upn = nB_upn + nB_up_[j];
		nB_dnn = nB_dnn + nB_dn_[j];
		
  } //for(int j = 0; j<2*GRID; j++)
 delete[] nA_up_;
 delete[] nB_up_;
 delete[] nA_dn_;
 delete[] nB_dn_;
 nA_upn = nA_upn/(4*GRID*GRID);
 nA_dnn = nA_dnn/(4*GRID*GRID);
 nB_upn = nB_upn/(4*GRID*GRID);
 nB_dnn = nB_dnn/(4*GRID*GRID);
 }	//solve

  long double get_GRID(){ return 2*t; };
  auto get_lambda(){return lambda_upp[0][0];}

  auto dos_(long double U, long double mu, long double nA_up, long double nA_dn, long double nB_up, long double nB_dn, long double min_energy, long double max_energy){
	auto en_step = (max_energy - min_energy)/ bin;
	int l;		
	for(int j = 0; j<2*GRID; j++){
		for(int i = 0; i<2*GRID; i++){
			auto energy = 2*t*(cos_[i] + cos_[j]);
			auto energyA = 2*dim*t2*cos_[i]*cos_[j];
			auto A_up = delta -mu + U*nA_up - energyA;
			auto A_dn = delta -mu + U*nA_dn - energyA;
			auto B_up = -delta -mu + U*nB_up - energyA;
			auto B_dn = -delta -mu + U*nB_dn - energyA;
			lambda_upp[i][j] = 0.5*( A_dn + B_dn + sqrt( (A_dn - B_dn)*(A_dn - B_dn) + 4*energy*energy) ) - min_energy;
			l = int(lambda_upp[i][j]/en_step);
			dos_upp[l] = dos_upp[l] + 1; 

			lambda_dnp[i][j] = 0.5*( A_up + B_up + sqrt( (A_up - B_up)*(A_up - B_up) + 4*energy*energy) )- min_energy;
			l = int(lambda_dnp[i][j]/en_step);
			dos_dnp[l] = dos_dnp[l] + 1;

			lambda_upn[i][j] = 0.5*( A_dn + B_dn - sqrt( (A_dn - B_dn)*(A_dn - B_dn) + 4*energy*energy) )- min_energy;
			l = int(lambda_upn[i][j]/en_step);
			dos_upn[l] = dos_upn[l] + 1;

			lambda_dnn[i][j] = 0.5*( A_up + B_up - sqrt( (A_up - B_up)*(A_up - B_up) + 4*energy*energy) )- min_energy;	
			l = int(lambda_dnn[i][j]/en_step);
			dos_dnn[l] = dos_dnn[l] + 1;
			}
		}
	 stringstream dos_file;
	 dos_file<<"dos"<<"delta"<<delta<<"U"<<U <<"t2"<<t2<<"h.txt";
         ofstream dos(dos_file.str());
         dos<< "# U,mu,nA_up, nA_dn, nB_up, nB_dn, ms, mf, ntotal,mu,cint"<<endl;
         dos << "#energy    dos_upn   dos_upp   dos_dnn  dos_dnp GRID="<<GRID<<" bin="<<bin<<endl;
         for(int p = 0; p<bin; p++){
	   dos_upp[p] = dos_upp[p]/(4*GRID*GRID);
	   dos_upn[p] = dos_upn[p]/(4*GRID*GRID);
	   dos_dnp[p] = dos_dnp[p]/(4*GRID*GRID);
	   dos_dnn[p] = dos_dnn[p]/(4*GRID*GRID);
	   dos << (min_energy + p*en_step) <<"     "<<dos_upn[p]<<"      "<<dos_upp[p]<<"      "<<dos_dnn[p]<<"     "<<dos_dnp[p]<<endl;
	}//for(p=)	
  }//dos

};//class

long double abs(long double a){
	if (a>0.0){
		return a;
	}
	else{
		return -1.0*a;
	}
	}

int main(){
  long double t,t2,delta,U,mu,nA_up,nA_dn,nB_up,nB_dn,conv;
  int dim,GRID,loop;
  t = 0.5;
  t2 = 0.0;
  delta = 0.0;
  dim = 2;
  GRID = 200;
  U = 0.1;
  mu = 0.1;
  conv = 0.000000001;

  HF hf(t,t2,delta,dim);

  cout<<hf.get_GRID()<<endl;
  auto lambda = hf.get_lambda();
  //cout<< GRID <<hf.GRID<< hf.t<<endl;
  //hf.lambda_upp[0][0] = 100.0;
  cout<<hf.get_lambda()<<endl;
  
  stringstream data_file;
  data_file<<"delta"<<delta<<"t2"<<t2<<"_.txt";
  ofstream soumen(data_file.str());
  soumen <<"# U,mu,nA_up, nA_dn, nB_up, nB_dn, ms, mf, ntotal,mu,cint,dnA_up, dnA_dn, dnB_up, dnB_dn,loop"<<endl;
  //0.471591	0.466898	0.527911	0.53287
  //0.963728   0.342376	0.275014	0.657618	0.724984
  nA_up = 0.5;
  nA_dn = 0.5;
  nB_up = 0.5;
  nB_dn = 0.5;
  vector<long double> U_list ={0.2,0.22,0.24,0.26,0.28};
  //vector<long double> U_list ={1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9};//,0.1,0.2,0.4,0.43,0.46,0.49,0.6};
  //vector<long double> U_list ={0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.56};
  for(long double U : U_list) {
   cout<<"mu.........."<<mu<<endl;
   //------------------------------
   mu = U/2.0;
   auto nloops= 2000;
   auto mix = 0.7;
   auto coeffi = 0.3;
   int l =0;
    nA_up = nA_up + 0.001;
    nA_dn = nA_dn - 0.001;
    nB_up = nB_up - 0.001;
    nB_dn = nB_dn + 0.001;
   for (loop = 0; loop<nloops;loop++){
      if (loop > 500) {mix=0.3; coeffi=0.2;}
      if (loop > 1000) {mix=0.2; coeffi=0.1;}
      hf.solve(U, mu, nA_up, nA_dn, nB_up,nB_dn);
      auto ntotal = hf.nA_upn + hf.nA_dnn + hf.nB_upn + hf.nB_dnn;
      if( abs(ntotal-2)< 0.0001  and abs(hf.nA_upn-nA_up)<conv and abs(hf.nA_dnn-nA_dn)<conv and abs(hf.nB_upn-nB_up)<conv and abs(hf.nB_dnn-nB_dn)<conv ){break;}
      
      nA_up = mix*hf.nA_upn + (1-mix)*nA_up;
      nA_dn = mix*hf.nA_dnn + (1-mix)*nA_dn;
      nB_up = mix*hf.nB_upn + (1-mix)*nB_up;
      nB_dn = mix*hf.nB_dnn + (1-mix)*nB_dn; 
   
   if (((ntotal-2)> 0.0001 ||(ntotal-2)< 0.0001) && l> 10 && loop >30 && loop < (nloops-20)){
		mu = mu - coeffi*(ntotal-2.0);
		l = 0;
		cout<<loop<<""<<U<<"  "<<hf.nA_upn<<"	"<<hf.nA_dnn<<"	"<<hf.nB_upn<<"	"<<hf.nB_dnn<<"	"<< ntotal<<"   "<<hf.cnt<<"  "<<mu<<" "<<GRID<<endl;
		}//if
   l = l+1;
    
   //cout<<loop<<""<<U<<"  "<<hf.nA_upn<<"	"<<hf.nA_dnn<<"	"<<hf.nB_upn<<"	"<<hf.nB_dnn<<"	"<< ntotal<<"   "<<hf.cnt<<"  "<<mu<<" "<<GRID<<endl;
   
   }//loop
  
 auto ntotal = nA_up + nA_dn + nB_up + nB_dn;
 auto mz= 0.5*(nA_up - nA_dn - nB_up + nB_dn);
 auto mf = 0.5*(nA_up - nA_dn + nB_up - nB_dn);
 
 soumen << U<<"  "<<mu<<"   "<<nA_up<<"	"<<nA_dn<<"	"<<nB_up<<"	"<<nB_dn<<"	"<< mz<<"	"<< mf<<"	"<< ntotal<<"   "<<mu <<"	"<<hf.cnt<<"	"<<abs(hf.nA_upn-nA_up)<<"	"<<abs(hf.nA_dnn-nA_dn)<<"	"<<abs(hf.nB_upn-nB_up)<<"	"<<abs(hf.nB_dnn-nB_dn)<<"	"<<loop<<endl;
 cout<<"  "<<U<<"  "<<mu<<"  "<<nA_up<<"	"<<nA_dn<<"	"<<nB_up<<"	"<<nB_dn<<"	"<< mz<<"	"<< mf<<"	"<< ntotal<<"   "<<hf.cnt<<"  "<<mu<<" "<<GRID<<endl;
//dos_(long double U, long double mu, long double nA_up, long double nA_dn, long double nB_up, long double nB_dn, long double min_energy, long double max_energy)
 hf.dos_(U, mu, nA_up, nA_dn, nB_up,nB_dn,-3.0,3.0);
 }//U
  
  return 0;

} //main
