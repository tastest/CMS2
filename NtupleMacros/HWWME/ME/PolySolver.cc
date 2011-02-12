#include "PolySolver.hpp"
#include "TMath.h"
/*
Using Sturm Sequences to Bracket Real Roots of Polynomial Equations
by D.G. Hook and P.R. McAree
from "Graphics Gems", Academic Press, 1990
*/

#include <stdio.h>
#include <math.h>
#include <iostream>
using std::cout;
using std::endl;
/*
 * a driver program for a root solver.
 */
/*
int main(){

 
 double epX= -0.635972;double  epY= 33.060679;double  epZ= 28.590870;double  epE= 43.713279;
 double emX=  7.550971;double  emY=-15.732424;double  emZ=105.253369;double  emE=106.690196;
 double nuX= 41.132373;double  nuY=-33.804731;double  nuZ= 67.206673;double  nuE= 85.740124;
 double nbX=-48.047372;double  nbY= 16.476476;double  nbZ= 40.928679;double  nbE= 65.231749;
 
 double wpX=epX+nuX;	double wmX=emX+nbX;
 double wpY=epY+nuY;	double wmY=emY+nbY;
 double wpZ=epZ+nuZ;	double wmZ=emZ+nbZ;
 double wpE=epE+nuE;	double wmE=emE+nbE;
 
 double Mw1=sqrt(wpE*wpE-wpX*wpX-wpY*wpY-wpZ*wpZ);
 double Mw2=sqrt(wmE*wmE-wmX*wmX-wmY*wmY-wmZ*wmZ);
  
 
 double A=Mw1*Mw1+epX*epX+epY*epY+epZ*epZ-epE*epE+2*epZ*nuZ;
 double B=2*epX;
 double C=2*epY;

 double p1=B*B-4*epE*epE;
 double p2=2*B*C;
 double p3=C*C-4*epE*epE;
 double p4=2*A*B;
 double p5=2*A*C;
 double p6=A*A-4*epE*epE*nuZ*nuZ;
 
 double Qx=-epX-emX;
 double Qy=-epY-emY;
 
 double D=Mw2*Mw2-emE*emE
         +epX*epX-Qx*Qx
	 +epY*epY-Qy*Qy
	 +emZ*emZ+2*emZ*nbZ;
	 
 double E=-2*emX;
 double F=-2*emY;

 double q1=E*E-4*emE*emE;
 double q2=2*E*F;
 double q3=F*F-4*emE*emE;
 double q4=2*D*E+8*Qx*emE*emE;
 double q5=2*D*F+8*Qy*emE*emE;
 double q6=D*D-4*emE*emE*nbZ*nbZ-4*emE*emE*(Qx*Qx+Qy*Qy);
 
// p1=1;p2=0;p3=1;p4=-2;p5=-2;p6=2;
// q1=1;q2=0;q3=1;q4=-2;q5=-2;q6=2;

 double M= q3*p1-p3*q1;
 double N= q3*p4-p3*q4;
 double O= q3*p6-p3*q6;
 double P=-q3*p2+p3*q2;
 double Q=-q3*p5+p3*q5;
 
 double R=     P*p2+M*p3;
 double S=P*p5+Q*p2+N*p3;
 double T=Q*p5     +O*p3;
 
 cout << M <<" " << N << " " 
      << O <<" " << P << " " 
      << Q <<" " << R << " " 
      << S <<" " << T << " " 
      <<endl;
      
 double k[5];
 k[4]=  	       P*P*p1+M*R;
 k[3]=        2*P*Q*p1+P*P*p4+M*S+N*R;
 k[2]= Q*Q*p1+2*P*Q*p4+P*P*p6+M*T+N*S+O*R;
 k[1]= Q*Q*p4+2*P*Q*p6  	 +N*T+O*S;
 k[0]= Q*Q*p6			     +O*T;
 
 
// double coeff[5]={-300*300*300*300,0,0,0,1};
 double limits[2]={-1000,1000};
 double roots[10];
 PolySolver::FindRoots(4,k,limits,roots);

 cout << "Input value: Mw1=" << Mw1 <<" Mw2=" <<Mw2 <<endl; 

 for(int i=0;i<4;i++){
   double x=roots[i];
   double y=(M*x*x+N*x+O)/(P*x+Q);
   double eqn=k[0]
             +k[1]*x
	     +k[2]*x*x
	     +k[3]*x*x*x
	     +k[4]*x*x*x*x;
   
   double eqn1=p1*x*x
              +p2*x*y
	      +p3*y*y
	      +p4*x
	      +p5*y
	      +p6;

   double eqn2=q1*x*x
              +q2*x*y
	      +q3*y*y
	      +q4*x
	      +q5*y
	      +q6;
	      
   cout <<"chek sol :" << x <<" y= " <<y 
        <<" StrumEqn= " <<eqn 
	<<" eqn1= " << eqn1
	<<" eqn2= " << eqn2
	<<endl;
   
   nuX=x;              nuY=y;             nuE=sqrt(nuX*nuX+nuY*nuY+nuZ*nuZ);
   nbX=-(epX+emX+nuX); nbY=-(epY+emY+nuY);nbE=sqrt(nbX*nbX+nbY*nbY+nbZ*nbZ); 
   
  wpX=epX+nuX;    wmX=emX+nbX;
  wpY=epY+nuY;    wmY=emY+nbY;
  wpZ=epZ+nuZ;    wmZ=emZ+nbZ;
  wpE=epE+nuE;    wmE=emE+nbE;
 
  Mw1=sqrt(wpE*wpE-wpX*wpX-wpY*wpY-wpZ*wpZ);   
  Mw2=sqrt(wmE*wmE-wmX*wmX-wmY*wmY-wmZ*wmZ);

cout <<"nuX= " << nuX <<" nuY= " << nuY <<" Wm1= " << Mw1 <<endl;
cout <<"nbX= " << nbX <<" nbY= " << nbY <<" Wm2= " << Mw2 <<endl;

 cout <<endl;
 }


 return 0;
}


*/

int PolySolver::FindRoots(int order, double coeff[], double limits[2],double roots[])
{

	poly	sseq[MAX_ORDER];
	double 	min, max;
	int		i, nroots, nchanges, np, atmin, atmax;
	
	//cout << "Coeff: ";
	for (i = order; i >= 0; i--) 
	  {	    
	    sseq[0].coef[i] = coeff[i];	    
	    //  cout << " " << sseq[0].coef[i];
	  }
	//	cout << endl;

	/*
	 * build the Sturm sequence
	 */
	np = buildsturm(order, sseq);

	/*
	printf("Sturm sequence for:\n");

	for (i = order; i >= 0; i--)
			printf("%lf ", sseq[0].coef[i]);

	printf("\n\n");

	for (i = 0; i <= np; i++) {
			for (int j = sseq[i].ord; j >= 0; j--)
				printf("%lf ", sseq[i].coef[j]);
			printf("\n");
	}

	printf("\n");
	*/
	/* 
	 * get the number of real roots
	 */
	nroots = numroots(np, sseq, &atmin, &atmax);

	//cout << "[] Number of roots: " << nroots << endl;
	if (nroots == 0) {
	  //	printf("solve: no real roots\n");
			return 0;
	}

	/*
	 * calculate the bracket that the roots live in
	 */

	/*
	min = -1.0;
	nchanges = numchanges(np, sseq, min);
	for (i = 0; nchanges != atmin && i != MAXPOW; i++) { 
			min *= 10.0;
			nchanges = numchanges(np, sseq, min);
			}
			*/
	// new from dw
	min = limits[0];
	nchanges = numchanges(np, sseq, min);

	if (nchanges != atmin) {
	  //	printf("solve: unable to bracket all negative roots\n");
			atmin = nchanges;
	}

	/*
	max = 1.0;
	nchanges = numchanges(np, sseq, max);
	for (i = 0; nchanges != atmax && i != MAXPOW; i++) { 
			max *= 10.0;
			nchanges = numchanges(np, sseq, max);
	}
	*/

	max = limits[1];
	nchanges = numchanges(np, sseq, max);

	if (nchanges != atmax) {
	  //	printf("solve: unable to bracket all positive roots\n");
			atmax = nchanges;
	}

	nroots = atmin - atmax;
	
	//cout << "Roots within limits: " << nroots << endl;
	/*
	 * perform the bisection.
	 */
	if (nroots) sbisect(np, sseq, min, max, atmin, atmax, roots);

	return nroots;

}

////////////////////////////////////////////// sturm.c

/*
 * modp
 *
 *	calculates the modulus of u(x) / v(x) leaving it in r, it
 *  returns 0 if r(x) is a constant.
 *  note: this function assumes the leading coefficient of v 
 *	is 1 or -1
 */
static int modp(poly *u, poly *v, poly *r)
{
	int		k, j;
	double	*nr, *end, *uc;

	nr = r->coef;
	end = &u->coef[u->ord];

	uc = u->coef;
	while (uc <= end)
			*nr++ = *uc++;

	if (v->coef[v->ord] < 0.0) {


			for (k = u->ord - v->ord - 1; k >= 0; k -= 2)
				r->coef[k] = -r->coef[k];

			for (k = u->ord - v->ord; k >= 0; k--)
				for (j = v->ord + k - 1; j >= k; j--)
					r->coef[j] = -r->coef[j] - r->coef[v->ord + k]
					* v->coef[j - k];
	} else {
			for (k = u->ord - v->ord; k >= 0; k--)
				for (j = v->ord + k - 1; j >= k; j--)
				r->coef[j] -= r->coef[v->ord + k] * v->coef[j - k];
	}

	k = v->ord - 1;
	while (k >= 0 && fabs(r->coef[k]) < SMALL_ENOUGH) {
		r->coef[k] = 0.0;
		k--;
	}

	r->ord = (k < 0) ? 0 : k;

	return(r->ord);
}

/*
 * buildsturm
 *
 *	build up a sturm sequence for a polynomial in smat, returning
 * the number of polynomials in the sequence
 */
int buildsturm(int ord, poly *sseq)
{
	int		i;
	double	f, *fp, *fc;
	poly	*sp;

	sseq[0].ord = ord;
	sseq[1].ord = ord - 1;


	/*
	 * calculate the derivative and normalise the leading
	 * coefficient.
	 */
	f = fabs(sseq[0].coef[ord] * ord);
	fp = sseq[1].coef;
	fc = sseq[0].coef + 1;
	for (i = 1; i <= ord; i++)
			*fp++ = *fc++ * i / f;

	/*
	 * construct the rest of the Sturm sequence
	 */
	for (sp = sseq + 2; modp(sp - 2, sp - 1, sp); sp++) {

		/*
		 * reverse the sign and normalise
		 */
		f = -fabs(sp->coef[sp->ord]);
		for (fp = &sp->coef[sp->ord]; fp >= sp->coef; fp--)
				*fp /= f;
	}

	sp->coef[0] = -sp->coef[0];	/* reverse the sign */

	return(sp - sseq);
}

/*
 * numroots
 *
 *	return the number of distinct real roots of the polynomial
 * described in sseq.
 */
int
numroots(int np, poly *sseq, int *atneg, int *atpos)
{
		int		atposinf, atneginf;
		poly	*s;
		double	f, lf;

		atposinf = atneginf = 0;


	/*
	 * changes at positive infinity
	 */
	lf = sseq[0].coef[sseq[0].ord];

	for (s = sseq + 1; s <= sseq + np; s++) {
			f = s->coef[s->ord];
			if (lf == 0.0 || lf * f < 0)
				atposinf++;
		lf = f;
	}

	/*
	 * changes at negative infinity
	 */
	if (sseq[0].ord & 1)
			lf = -sseq[0].coef[sseq[0].ord];
	else
			lf = sseq[0].coef[sseq[0].ord];

	for (s = sseq + 1; s <= sseq + np; s++) {
			if (s->ord & 1)
				f = -s->coef[s->ord];
			else
				f = s->coef[s->ord];
			if (lf == 0.0 || lf * f < 0)
				atneginf++;
			lf = f;
	}

	*atneg = atneginf;
	*atpos = atposinf;

	return(atneginf - atposinf);
}

/*
 * numchanges
 *
 *	return the number of sign changes in the Sturm sequence in
 * sseq at the value a.
 */
int
numchanges(int np, poly *sseq, double a)
{
	int		changes;
	double	f, lf;
	poly	*s;

	changes = 0;

	lf = evalpoly(sseq[0].ord, sseq[0].coef, a);

	for (s = sseq + 1; s <= sseq + np; s++) {
			f = evalpoly(s->ord, s->coef, a);
			if (lf == 0.0 || lf * f < 0)
				changes++;
			lf = f;
	}

	return(changes);
}

/*
 * sbisect
 *
 *	uses a bisection based on the sturm sequence for the polynomial
 * described in sseq to isolate intervals in which roots occur,
 * the roots are returned in the roots array in order of magnitude.
 */
void sbisect(int np, poly *sseq, double min, double max, int atmin, int atmax, double *roots)
{
	double	mid;
	int		n1 = 0, n2 = 0, its, atmid, nroot;

	if ((nroot = atmin - atmax) == 1) {

		/*
		 * first try a less expensive technique.
		 */
	  //if (modrf(sseq->ord, sseq->coef, min, max, &roots[0]))
	  //	return;


		/*
		 * if we get here we have to evaluate the root the hard
		 * way by using the Sturm sequence.
		 */
		for (its = 0; its < MAXIT; its++) {
				mid = (min + max) / 2;

				atmid = numchanges(np, sseq, mid);

				if (fabs(mid) > RELERROR) {
					if (fabs((max - min) / mid) < RELERROR) {
						roots[0] = mid;
						return;
					}
				} else if (fabs(max - min) < RELERROR) {
					roots[0] = mid;
					return;
				}

				if ((atmin - atmid) == 0)
					min = mid;
				else
					max = mid;
			}

		if (its == MAXIT) {
				fprintf(stderr, "sbisect: overflow min %f max %f\
					diff %e nroot %d n1 %d n2 %d\n",
					min, max, max - min, nroot, n1, n2);
			roots[0] = mid;
		}

		return;
	}

	/*
	 * more than one root in the interval, we have to bisect...
	 */
	for (its = 0; its < MAXIT; its++) {

			mid = (min + max) / 2;

			atmid = numchanges(np, sseq, mid);

			n1 = atmin - atmid;
			n2 = atmid - atmax;


			if (n1 != 0 && n2 != 0) {
				sbisect(np, sseq, min, mid, atmin, atmid, roots);
				sbisect(np, sseq, mid, max, atmid, atmax, &roots[n1]);
				break;
			}

			if (n1 == 0)
				min = mid;
			else
				max = mid;
	}

	if (its == MAXIT) {
			fprintf(stderr, "sbisect: roots too close together\n");
			fprintf(stderr, "sbisect: overflow min %f max %f diff %e\
				nroot %d n1 %d n2 %d\n",
				min, max, max - min, nroot, n1, n2);
			for (n1 = atmax; n1 < atmin; n1++)
			roots[n1 - atmax] = mid;
	}
}

////////////////////////////////////////////////////// util.c

/*
 * evalpoly
 *
 *	evaluate polynomial defined in coef returning its value.
 */
double evalpoly (int ord, double *coef, double x)
{
	double	*fp, f;

	fp = &coef[ord];
	f = *fp;

	for (fp--; fp >= coef; fp--)
	f = x * f + *fp;

	return(f);
}


/*
 * modrf
 *
 *	uses the modified regula-falsi method to evaluate the root
 * in interval [a,b] of the polynomial described in coef. The
 * root is returned is returned in *val. The routine returns zero
 * if it can't converge.
 */
int modrf(int ord, double *coef,double a, double b, double *val)
{
	int		its;
	double	fa, fb, x, fx, lfx;
	double	*fp, *scoef, *ecoef;

	scoef = coef;
	ecoef = &coef[ord];

	fb = fa = *ecoef;
	for (fp = ecoef - 1; fp >= scoef; fp--) {
		fa = a * fa + *fp;
		fb = b * fb + *fp;
	}

	/*
	 * if there is no sign difference the method won't work
	 */
	if (fa * fb > 0.0)
		return(0);

	if (fabs(fa) < RELERROR) {
		*val = a;
		return(1);
	}

	if (fabs(fb) < RELERROR) {
		*val = b;
		return(1);
	}

	lfx = fa;


	for (its = 0; its < MAXIT; its++) {

		x = (fb * a - fa * b) / (fb - fa);

		fx = *ecoef;
		for (fp = ecoef - 1; fp >= scoef; fp--)
				fx = x * fx + *fp;

		if (fabs(x) > RELERROR) {
				if (fabs(fx / x) < RELERROR) {
					*val = x;
					return(1);
				}
		} else if (fabs(fx) < RELERROR) {
				*val = x;
				return(1);
		}

		if ((fa * fx) < 0) {
				b = x;
				fb = fx;
				if ((lfx * fx) > 0)
					fa /= 2;
		} else {
				a = x;
				fa = fx;
				if ((lfx * fx) > 0)
					fb /= 2;
		}

		lfx = fx;
	}

	fprintf(stderr, "modrf overflow %f %f %f\n", a, b, fx);

	return(0);
}
	


int PolySolver::FindRoots4thOrder(double coeff[], double roots[],double errFlag[])
{
 
 double minVal=-1E-5;
 int Flag=0;
 double a0=(double)coeff[0]/coeff[4];
 double a1=(double)coeff[1]/coeff[4];
 double a2=(double)coeff[2]/coeff[4];
 double a3=(double)coeff[3]/coeff[4];

// cout << a0 << a1 << a2 << a3<<endl;
 double B1=12*a0+a2*a2-3*a1*a3;
 double B2=27*a1*a1-72*a0*a2+2*a2*a2*a2-9*a1*a2*a3+27*a0*a3*a3*a3;
 double B3=-4*B1*B1*B1+B2*B2;
 if(B3>0) B3=TMath::Sqrt(B3); else if(B3>minVal){ B3=0;} else { B3=0; Flag=1;}
 double G=B2+B3;

// cout << " B " << B2 <<" " << B3 <<" " <<G <<endl;
 if(G <=0) {
     Flag=2;	 
     errFlag[0]=errFlag[1]=errFlag[2]=errFlag[3]=Flag;
     roots[0]=roots[1]=roots[2]=roots[3]=0;
cout <<"Flag: " << Flag << " G=B2+B3=" <<B2<<"+"<<B3<<"="<< G<<endl;
cout <<"coef: " << coeff[0]/coeff[4] 
          <<" " << coeff[1]/coeff[4] 
          <<" " << coeff[2]/coeff[4] 
          <<" " << coeff[3]/coeff[4] 
	  <<endl;
     return Flag;
	 
 }
 G=TMath::Power(G ,1./3); 

// cout << G <<endl; 
 //real solution for cubic equations
 double y1=a2/3+ (4*TMath::Power(2,1./3)*B1/G +TMath::Power(2,5./3)*G)/12;

 
 double R=a3*a3/4-a2+y1;
 if(R >0) R=sqrt(R); else if(R>minVal){ R=0;} else { 
  Flag=3;
  errFlag[0]=errFlag[1]=errFlag[2]=errFlag[3]=Flag;
  roots[0]=roots[1]=roots[2]=roots[3]=0;

  cout <<"Flag: " << Flag << " R=" <<R<<endl;

  return Flag;
 } 
 
 double D=0; double E=0;
 if(R ==0){

   double D1=y1*y1-4*a0;
   if(D1>0) D1=sqrt(D1); else if(D1>minVal){ D1=0;} else { 
	    errFlag[0]=errFlag[1]=errFlag[2]=errFlag[3]=4;
	    roots[0]=roots[1]=roots[2]=roots[3]=0;

	    cout <<"Flag: " << Flag << " D1=" <<D1<<endl;

	    return 4;
   }
   D=3*a3*a3/4-2*a2+2*D1;
   E=3*a3*a3/4-2*a2-2*D1;
   
 }else{
   D=3*a3*a3/4-R*R-2*a2+(4*a2*a3-8*a1-a3*a3*a3)/4/R;
   E=3*a3*a3/4-R*R-2*a2-(4*a2*a3-8*a1-a3*a3*a3)/4/R;
   
 }
 
 
 if(D>0) D=TMath::Sqrt(D); else if(D/fabs(coeff[4])>minVal){ D=0;} else {
cout <<"Flag 5 D=" <<D<<endl;
cout <<"coef: " << coeff[0]/coeff[4] 
          <<" " << coeff[1]/coeff[4] 
          <<" " << coeff[2]/coeff[4] 
          <<" " << coeff[3]/coeff[4] 
	  <<endl;
  D=0; errFlag[0]=errFlag[1]=5; 
  
  }
 if(E>0) E=TMath::Sqrt(E); else if(E/fabs(coeff[4]>minVal)){ E=0;} else {

cout <<"Flag 5 E=" <<E<<endl;
cout <<"coef: " << coeff[0]/coeff[4] 
          <<" " << coeff[1]/coeff[4] 
          <<" " << coeff[2]/coeff[4] 
          <<" " << coeff[3]/coeff[4] 
	  <<endl;
  E=0; errFlag[2]=errFlag[3]=5;   
  }
 
 roots[0]=(-a3/2+R+D)/2;
 roots[1]=(-a3/2+R-D)/2;
 roots[2]=(-a3/2-R+E)/2;
 roots[3]=(-a3/2-R-E)/2;

 return Flag;
 
}



