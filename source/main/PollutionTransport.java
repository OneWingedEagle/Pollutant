package main;
import math.*;
import static java.lang.Math.*;
import io.Console;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Scanner;

import javax.swing.JOptionPane;

import org.apache.commons.math3.special.Erf;

import java.io.FileDescriptor;
import java.io.FileOutputStream;

import components.GUI;


	public class PollutionTransport {
		// 1 D 
		
		boolean steady=false;
		
		double x0=0, xn=10000;
		double t0=0, tn;
		int xDiv=200;
		double dx, dt;
		int nT1, nTsub;
		double[] knownValues;
		int[] unIndex;
		int[] un_kn_map;
		double V, D;


		public static void main(String[] args){
			
			new PollutionTransport();
		}

		public PollutionTransport(){
		//	x0=0;
		//	xn=10000;
			xDiv=200;
			 nT1=10;
			 nTsub=10;
			if(steady){
				 nT1=1;
				 nTsub=1;
			}


			 V=4e-3;
			 D=2e-1;
		}
		
		public void Run()
		{		
		
			 
		
		int nT=nT1*nTsub;

		 dt=(tn-t0)/nT;

		 dx=(xn-x0)/xDiv;

		 int N=xDiv;
		
		knownValues=new double[N];
		un_kn_map=new int[N];
		unIndex=new int[N];
		
		for(int i=0;i<N;i++){
			un_kn_map[i]=-1;
			unIndex[i]=-1;
		}

		int nunx=0;
		
		for(int i=0;i<N;i++){

			if(i>=0 && i<1){
				knownValues[i]=1;
			}else
				if(i>=N-1 && i<-N){
					knownValues[i]=.5;
				}else{
				unIndex[i]=nunx;
				un_kn_map[nunx]=i;
				nunx++;
			}
			
		}
		
		Vect qC0=new Vect(N);
		
		for(int i=0;i<N;i++){

		//	if(i>=95 && i<=105)
		///	qC0.el[100]=-1.6e-4;
					
		}
		//qC0.el[N-1]=-1.2e-3;
		
		//qC0.el[50]=-1.0e-4;
		

		double Cr =V*dt/dx; // Cre should be less than 1.
		double Pe =V*dx/D; // Pe should be less than 2.
		util.pr("Cr = "+Cr+ "     *****  Cre = V*dt/dx  should be <= 1.");
		util.pr("Pe = "+Pe+ "     *****  Pe = V*dx/D  should be <= 2.");
		
		double C0=1;
		
			double dx2=dx*dx;
		double at=1/dt;
		double am=-(V/(2*dx)+D/dx2);
		double a=2*D/dx2;
		double an=(V/(2*dx)-D/dx2);

		Mat H=new Mat(N, N);
		
		for(int i=0;i<N;i++){

			if(i>0)
				H.el[i][i-1]=am;
			
			H.el[i][i]=a;
			
		if(i<N-1) H.el[i][i+1]=an;

			if(nT>1)
				H.el[i][i]+=at;

		}
		

		int nun=nunx;
	//	H.show();

		Mat H2=new Mat(nun, nun);

		for(int k=0;k<N;k++){
			if(unIndex[k]>=0) {
			Vect v=new Vect(nun);

			for(int j=0;j<N;j++){
				if(unIndex[j]>=0){
				v.el[unIndex[j]]=H.el[k][j];
				H.el[k][j]=0;
				}
			}
			//v.hshow();
			H2.setRow(v, unIndex[k]);
			}
			}
//H2.show();
//util.hshow(un_kn_map);
//util.hshow(unIndex);
		Vect Cx=new Vect(N);
		for(int i=0;i<N;i++){
			if(unIndex[i]<0)
				Cx.el[i]=knownValues[i];
		}
	
//		Cx.show();
//H.show();
		Vect bx=H.mul(Cx);
		
		
//	bx.show();
		Vect b0=new Vect(nun);
		Vect qC=new Vect(nun);
		for(int i=0;i<N;i++){
			if(unIndex[i]>=0){
				b0.el[unIndex[i]]=-bx.el[i];
				qC.el[unIndex[i]]=-qC0.el[i];
			}
		}

		Mat M=H2;
		
		Vect[] C=new Vect[nT1];
		
		
		C[0]= new Vect(N);
		for(int k=0;k<N;k++)
			C[0].el[k]=Cx.el[k];
		
		Vect sol=new Vect(nun);
			
		Vect b=new Vect(nun);
		Mat Minv=M.inv();

		int m=0;
		for(int m1=0;m1<nT1;m1++)
			for(int m2=0;m2<nTsub;m2++){
			m++;
			double t=m*dt;
			double f=1;//exp(-1*t*.000002);//+.9*cos(t*m/N*2*4*PI);
			//if(t>5*86400) f=0;
		//	if(t>src_duration)
			//	b=sol.times(at);
		//	else
				b=b0.add(qC.times(f)).add(sol.times(1*at));	
			
			sol=Minv.mul(b);

	//		sol.hshow();

			if(m2==0 && (steady ||m1>0)){
			 C[m1]= new Vect(N);
			// if(t<=src_duration)

				for(int k=0;k<N;k++){
					if(unIndex[k]>=0)
						C[m1].el[k]=sol.el[unIndex[k]];
					else
						C[m1].el[k]=knownValues[k];
				}
			 
		//	 C[m1].hshow();
			}
			}
//
		//analytical
		
		boolean analyt=false;

	if(analyt){	
		Mat y=new Mat(nT1,N);
	
		for(int m1=0;m1<nT1;m1++){
			double t=m1*nTsub*dt;
			if(t==0) t=1e-6;
		for(int k=0;k<N;k++){
			double x=k*dx;
			double tt=C0*.5*(Erf.erfc((x-V*t)/(2*sqrt(D*t)))+Math.exp(V*x/D)*Erf.erfc((x+V*t)/(2*sqrt(D*t))));;
	
			y.el[m1][k]=tt;
		}
		}
		
	Mat AN=new Mat(N,1+nT1);

	for(int k=0;k<N;k++){
		AN.el[k][0]=k*dx;
		for(int m1=0;m1<nT1;m1++){
			AN.el[k][m1+1]=y.el[m1][k];
		}
		

	}
	String[] titles1=new String[AN.nCol];
	for(int k=0;k<titles1.length;k++)
		titles1[k]="AN, t= "+(k)*nTsub*dt;
	
	//util.plotBunch(AN.el,titles1);
	
	}
	Mat NUM=new Mat(N,nT1+1);
	
	for(int k=0;k<N;k++){
		NUM.el[k][0]=k*dx;
		
		for(int m1=0;m1<nT1;m1++){
			NUM.el[k][m1+1]=C[m1].el[k];
		}
	//	X.el[k][1+nT]=y.el[k];
	}
	
	
	NUM.show("%15.5e");
	

	
	
	String[] titles2=new String[NUM.nCol];
	for(int k=0;k<titles2.length;k++)
		titles2[k]="t= "+(k)+" days";


	util.plotBunch(NUM.el,titles2);

		}	
		
	}

