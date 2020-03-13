// This is the C++ implementation for the ForTran Program of Fukushima
#pragma once
class Node_Halley
{
public:	

	int ID   ;
	long double rp   ;	// sqrt(Xp*Xp + Yp*Yp)
	long double Zp   ;
	// Newton Facotrs
	long double s0   ;  // z/a
	long double s02  ;  // s0_square
	long double s03  ;  // s0_cube

    long double pn   ;  // r/a

	long double c0   ;  // b*r/aa
	long double c02  ;  // c0_square
	long double c03  ;  // c0_cube

    long double zc   ;  // s0*b/a

	long double a0   ;  // sqrt(a02)
	long double a02  ;  // c02+s02
	long double a03  ;  // a0_cube	   

    long double d0   ;  // d0=zc*a03+ee*s03
    long double f0   ;  // pn*a03-ee*c03

	// Prepare Halley Correction Factor

	long double b0   ; // e4T*s02*c02*pn*(a0-ec)
	long double s1   ; // d0*f0-b0*s0
	long double cc   ; // ec*(f0*f0-b0*c0)
    
	//	Evaluate phi & h

	long double fi_Halley_rad; // atan(s1/cc)  - in radians
	long double fi_Halley_DD;  // decimal degree
	long double s12      ; // s1*s1
	long double cc2      ; // cc*cc
	long double h_Halley ; // (p*cc+absz*s1-a*sqrt(ec2*s12+cc2))/sqrt(s12+cc2)
	long double h_Halley_err_nm ; // h_Halley - h	in nm
	long double h_Halley_err_mm ; // h_Halley - h	in mm
	long double fi_Halley_err_ss ; // fi_Halley_DD - fi	   - in sconds
	// long double Delta_Halley_err ; // abs(fi_Halley_err_ss * 1000) + abs(h_Halley_err_mm/1000)/(h+a) 

public:
	Node_Halley(void);
	~Node_Halley(void);
	void set_ID(int i)     ;
	void compute_s0(long double a) 		; 
	void compute_s02()		; 
	void compute_s03()		;
	void compute_pn(long double a) 		;
	void compute_zc(long double a, long double b) 	;
	void compute_c0(long double a, long double b) 		; 
	void compute_c02()   	;
	void compute_c03()		; 	 
	void compute_a02()		;
	void compute_a0()		;
	void compute_a03()		; 
	void compute_d0(long double ee)		; 
	void compute_f0(long double ee)		;	 
	void compute_b0(long double a, long double b, long double ee)		; 
	void compute_s1()		; 
	void compute_cc(long double a, long double b)		; 
	void compute_fi_Halley(); 
	void compute_s12()		;
	void compute_cc2()		;
	void compute_h_Halley(long double a, long double ee) ;
	void solve_Halley(long double a, long double b, long double ee);
	void compute_h_Halley_err(long double h)	;
	void compute_fi_Halley_err(long double fi);
	//void compute_Delta_Halley_err(long double a)		; 
	void write_Halley(char ff[20]) ;
	void write_Halley_Trace(char ff[20]) ;

};

