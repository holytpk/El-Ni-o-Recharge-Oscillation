//Euler Approximation for El Nino Southern Oscillation Model with Recharge Indecies and Random Factors. 
// L. He 06/28/2016

using namespace std;
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <cstdlib> 
#define USE_MATH_DEFINES

double noise(double);
double GD(double x);

int main(int argc,char *argv[])

{
   	int i, j, k, tau, N, Nt, seed, kmonth, ensemble, Nensemble = 10, tmax = atoi(argv[1]);
   	double t, dt, t0, *Tt, Tt0, Tti, Ttii, Tsum, *Taverage, *ht, ht0, hti, htii, hsum, *haverage, *r, w, x, *zt, zsum, *zaverage, B, G;; //declare time series variables
	double u, v, sumT, sumstdTe, *ensemblestd; //declare correlation coefficient variables
		
	//Time Rearrangement
	N = tmax*1000;
	
	r = new double [300000]();

	Tt = new double [3000000]();   //temperature as function of time
	ht = new double [3000000]();   //thermocline as fuction of time
	zt = new double [3000000]();   //noise as function of time

	Taverage = new double [3000000](); //monthly average
	haverage = new double [3000000]();
	zaverage = new double [3000000]();

	ensemblestd = new double [3000000]();

	ofstream outfile;  
	outfile.open("osc2010.dat");
	outfile.precision(5); 

	for(ensemble=1;ensemble<=Nensemble;ensemble++){

		seed = 7296616 + ensemble;
		srand48(seed); 
	
		Tt0 = 1.0;                //initial temperature, Celsius
		ht0 = 0.0;                //initial thermocline, meters
		t0 = 0.0;                 //initial time, month
		dt = 0.001;		  //time interva1, month
		Nt = (int)(1.0/dt);	   	  //number of time interval in a month

		r[k] = 1.0/6.0;		   //recharge index, month^-1	
		w = 2.0*M_PI/48.0;	   //annual frequency, month^-1
		B = 1.2;
		G = 1.0; 		   //ENSO SST anomalies
	
		u = 0.0; 		   //Average Temperature Throughout the Ensembles
		
		Tti = Tt0;         //initial conditions
		hti = ht0;
		Tt[1]= 0.1;
		ht[1]= 0.0;
		t = t0;

		//Iteration-dependent Frame
		for(k=1; k<=tmax; k++){
		
			r[k] = 1.0/6.0;
			
			for(i=1; i<=Nt; i++){			  	     
				Tt[i+(k-1)*Nt+1] = Tt[i+(k-1)*Nt] + dt*(-r[k]*Tt[i+(k-1)*Nt]+w*ht[i+(k-1)*Nt]+G*zt[i+(k-1)*Nt]/sqrt(24.0));        
				ht[i+(k-1)*Nt+1] = ht[i+(k-1)*Nt] + dt*(-w*Tt[i+(k-1)*Nt]);	                               	
				zt[i+(k-1)*Nt+1] = zt[i+(k-1)*Nt] + dt*(-zt[i+(k-1)*Nt]/1.5+noise(x));
				outfile << ensemble <<" "<< k <<" "<< i+(k-1)*Nt <<" "<< Tt[i+(k-1)*Nt] <<" "<< ht[i+(k-1)*Nt] <<" "<< zt[i+(k-1)*Nt] <<" "<< noise(x) << endl;
			}
			//cout << k <<" "<< r[k] << endl;
			}
			
		//Monthly Average
		for(i=1; i<=tmax; i++){
			Tsum = 0.0;
			hsum = 0.0;
			zsum = 0.0;
			for(j=1; j<=Nt; j++){
				Tsum += Tt[j+(i-1)*Nt];		
				hsum += ht[j+(i-1)*Nt];	
				zsum += zt[j+(i-1)*Nt];	
				//outfile << j <<" "<< Tsum <<" "<< hsum << endl;	
			}

			Taverage[i+(ensemble-1)*tmax] = Tsum/Nt; 	 	       		//average T
			haverage[i+(ensemble-1)*tmax] = hsum/Nt; 		       		//average h
			zaverage[i+(ensemble-1)*tmax] = zsum/Nt;  			 		//average z	
			//cout << i <<" "<< Taverage[i] << endl;
		}
	} //End of Thousands of Ensembles   
	
	for(kmonth=1; kmonth<=tmax; kmonth++){
		
		u = 0.0;		//reset STD terms
		sumT = 0.0;
		sumstdTe = 0.0;	
		
		for(i=kmonth; i<=Nensemble*tmax; i+=tmax){
			sumT += Taverage[i];
		}	
			
		u = sumT/(Nensemble*tmax);
  
		for(i=kmonth;i<=Nensemble*tmax;i+=tmax){
	
			sumstdTe += pow(Taverage[i]-u, 2.0);

		}

		//ensemblestd[kmonth] = pow(sumstdTe/(Nensemble-1), 0.5);
		//cout << kmonth <<" "<< ensemblestd[kmonth] << endl; 
		
	    }

	delete [] r;
  	delete [] Tt;	
   	delete [] ht;	
	delete [] zt;
	delete [] Taverage;
	delete [] haverage;
	delete [] zaverage;
	delete [] ensemblestd;
   	return 0;	
}

//VonNeumann Method 
double noise(double)  //noise, month^-1
{
	double x, W, Px, Pxmax, xmax;
	Pxmax = 1.0;  
	xmax = 10.0; //maximum value of noise
	W = Pxmax;
	Px = -99.0;
	while(W>Px){   
		W = Pxmax*(drand48());  
		x = xmax*2.0*(drand48()-0.5);  
		Px = GD(x);
	}
	return(0);
}

//Gaussian Distribution
double GD(double x)
{
	double xprime, mean, SD;
	SD = 30.0;
	mean = 0.0; 	
	xprime = x - mean; 
	return(exp(-xprime*xprime/(2.0*SD*SD))/(SD*sqrt(2.0*M_PI))); //Probability Distribution
}
