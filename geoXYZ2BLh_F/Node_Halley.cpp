#include "stdafx.h"
#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#define _USE_MATH_DEFINES // for C++  
#include <cmath> 
#include "geoXYZ2BLh_F.h"
#include "Node_Halley.h"
using namespace std;


Node_Halley::Node_Halley(void)
{
	s0 			= 0.0; 
	s02			= 0.0; 
	s03	   		= 0.0;
	pn 			= 0.0;
	zc 			= 0.0;
	c0 			= 0.0; 
	c02			= 0.0;
	c03			= 0.0; 
	a0 			= 0.0; 
	a02			= 0.0; 
	a03			= 0.0; 
	d0 			= 0.0; 
	f0 			= 0.0 ;	 	
	b0 			= 0.0; 
	s1 			= 0.0; 
	cc 			= 0.0; 
	fi_Halley_DD	= 0.0; 
	fi_Halley_rad	= 0.0;
	h_Halley_err_nm	= 0.0; 
	h_Halley_err_mm	= 0.0;
	fi_Halley_err_ss = 0.0;
	s12			= 0.0 ;
	cc2			= 0.0 ;
	h_Halley	= 0.0 ;
	//Delta_Halley_err = 0.0;
}


Node_Halley::~Node_Halley(void)
{
}

void Node_Halley::set_ID(int i)
{
	ID = i+1;
	return;
}

void Node_Halley::compute_s0(long double a) 	
{
	// absz=abs(z)
	// ainv=1.d0/a
	// s0=absz*ainv

	s0 = Zp/a ;		
	return;
}
void Node_Halley::compute_s02()
{
	s02 = s0 * s0 ;
	return;
}
void Node_Halley::compute_s03()
{
	s03 = s02 * s0 ;
	return;
}
void Node_Halley::compute_pn(long double a)
{
	// pn=p*ainv
	pn = rp/a ;
	return;
}
void Node_Halley::compute_zc(long double a, long double b)
{		
	//zc=ec*s0 
	long double ec = b/a; 	 
	zc= s0 * ec ;
	return;
}
void Node_Halley::compute_c0(long double a, long double b)
{
	// c0 = (b * rp)/(a*a)	 ;
	long double ec = b/a;
	c0 = pn * ec ;
	return;
}
void Node_Halley::compute_c02()
{
	c02 = c0 * c0 ;
	return;
}
void Node_Halley::compute_c03()
{
	c03 = c02 * c0 ;
	return;
}
void Node_Halley::compute_a02()
{
	a02 = c02 + s02 ;
	return;
}
void Node_Halley::compute_a0()
{
	a0 = sqrt(a02) ;
	return;
}
void Node_Halley::compute_a03()
{
	a03 = a02 * a0 ;
	return;
}
void Node_Halley::compute_d0(long double ee)
{
	// d0 =zc*a03+e2*s03
	d0 = zc * a03 + ee * s03 ;
	return;
}
void Node_Halley::compute_f0(long double ee)
{
	// f0 =pn*a03-e2*c03
	f0 = pn *a03 - ee * c03 ;
	return;
}
void Node_Halley::compute_b0(long double a, long double b, long double ee)
{		
	// b0=e4T*s02*c02*pn*(a0-ec)
	long double ec = b/a;
	long double e4T = (ee * ee * 1.5)  ;
	long double kk = (a0-ec);
	b0 = e4T * s02 * c02 * pn * kk;
	return;
}
void Node_Halley::compute_s1()
{
	// s1=d0*f0-b0*s0  	
	s1 = (d0 * f0) - (b0 * s0) ;
	return;
}
void Node_Halley::compute_cc(long double a, long double b)		
{
	// cc=ec*(f0*f0-b0*c0)
	long double ec = b/a;
	long double kk = (f0 * f0) - (b0 * c0);
	cc = ec * kk ;
	return;
}
void Node_Halley::compute_fi_Halley()
{
	fi_Halley_rad = atan(s1/cc) ;
	fi_Halley_DD = 180 * fi_Halley_rad / PI	 ;
	return;
}
void Node_Halley::compute_s12()
{
	s12 = s1 * s1 ;
	return;
}
void Node_Halley::compute_cc2()
{
	cc2 = cc * cc ;
	return;
}
void Node_Halley::compute_h_Halley(long double a, long double ee)
{
	// h=(p*cc+absz*s1-a*sqrt(ec2*s12+cc2))/sqrt(s12+cc2)
	long double ec2 = 1 - ee ;
	long double k1 = rp * cc ;
	long double k2 = Zp * s1 ; 
	long double k3 = sqrt(ec2 * s12 + cc2)  ;
	long double k4 = a * k3;  	
	long double k5 = sqrt(s12 + cc2);
	
	h_Halley = (k1 + k2 -k4) / k5;
	return;
}

void Node_Halley::solve_Halley(long double a, long double b,  long double ee )
{
	compute_s0(a) 		; 
	compute_s02()		; 
	compute_s03()		;
	compute_pn(a) 		;
	compute_zc(a, b) 	;
	compute_c0(a, b) 	; 
	compute_c02()   	;
	compute_c03()		; 	 
	compute_a02()		;
	compute_a0()		;
	compute_a03()		; 
	compute_d0(ee)		; 
	compute_f0(ee)		;	 
	compute_b0(a, b, ee); 
	compute_s1()		; 
	compute_cc(a, b)	; 
	compute_fi_Halley() ; 
	compute_s12()		;
	compute_cc2()		;
	compute_h_Halley(a, ee) ;
	
	return;
}

void Node_Halley::compute_h_Halley_err(long double h)
{
	h_Halley_err_mm = 1000 * (h_Halley - h) ;
	h_Halley_err_nm = 1000000000 * (h_Halley - h) ;
	return;
}

void Node_Halley::compute_fi_Halley_err(long double fi)
{
	fi_Halley_err_ss = 3600 * (fi_Halley_DD - fi) ;
	return;
}

void Node_Halley::write_Halley(char ff[20])
{
	ofstream fn;
	fn.open(ff,ios::app); // open file in append mode
	fn << std::fixed << std::setprecision(18) ;
		
	fn << ID <<"," << fi_Halley_DD <<"," << h_Halley <<"," << (fi_Halley_err_ss * 1e6) <<"," ;
	fn << (h_Halley_err_mm * 1e6) ;
	fn << endl	;
	fn.close();
	return;
}

void Node_Halley::write_Halley_Trace(char ff[20])
{
	ofstream fn;
	fn.open(ff,ios::app); // open file in append mode
	fn << std::fixed << std::setprecision(18) ;
		
	//fn << ID <<"," << rp <<"," << Zp <<"," << s0 <<"," << a0;
	//fn << endl	;	
	fn << "* * * * * * * * * * * * * *    Halley Correction Factors for Node ID =" << ID  << endl ; 
    fn << "rp =" << rp  << endl ;
    fn << "Zp =" << Zp  << endl ;
    fn << "s0 =" << s0  << endl ;
    fn << "s02=" << s02 << endl ;
    fn << "s03=" << s03 << endl ;
    fn << "pn =" << pn  << endl ;
    fn << "c0 =" << c0  << endl ;
    fn << "c02=" << c02 << endl ;
    fn << "c03=" << c03 << endl ;
    fn << "zc =" << zc  << endl ;
    fn << "a0 =" << a0  << endl ;
    fn << "a02=" << a02 << endl ;
    fn << "a03=" << a03 << endl ;
    fn << "d0 =" << d0  << endl ;
    fn << "f0 =" << f0  << endl ;
    fn << "b0 =" << b0  << endl ;
    fn << "s1 =" << s1  << endl ;
    fn << "cc =" << cc  << endl ; 
    fn << "s12=" << s12 << endl ;
    fn << "cc2=" << cc2 << endl ;
    fn << "h_Halley=" << h_Halley  << endl ;     
    fn << "h_Halley_err_mm=" << h_Halley_err_mm << endl ; 
	fn << "h_Halley_err_nm=" << h_Halley_err_nm << endl ;
    fn << "fi_Halley_err_ss=" << fi_Halley_err_ss <<endl ;
    fn << "fi_Halley_rad="    << fi_Halley_rad <<endl ;
    fn << "fi_Halley_DD="     << fi_Halley_DD <<endl ;
	//fn << "Delta_Halley_err="     << Delta_Halley_err <<endl ;
	fn.close();
	return;
}

//void Node_Halley::compute_Delta_Halley_err(long double a)
//{
//	long double qq, vv;
//	qq =0.0 ;
//	vv = 0.0;
//
//  qq = abs(fi_Halley_err_ss) ;	  // 
//  vv = 0.001 * abs(h_Halley_err_mm);	   // in meters
//  Delta_Halley_err = qq + (vv/(h_Halley+a)) ;
//  return;
//}