#include <stdio.h>
#include <math.h>

// Declaration of the Cohesion5 and Cohesion6 functions
double Cohesion5(double *s1, double *s2, double phi, int dim, int lg);
double Cohesion6(double *s1, double *s2, double phi, int dim, int lg);

// Implementation of the Cohesion5 function
double Cohesion5(double *s1, double *s2, double phi, int dim, int lg) {
    int ii;
    double cent1 = 0, cent2 = 0, sdist = 0, out;
    for (ii = 0; ii < dim; ii++) {
        cent1 += s1[ii] / (double) dim;
        cent2 += s2[ii] / (double) dim;
    }

    for (ii = 0; ii < dim; ii++) {
        sdist += sqrt((s1[ii] - cent1) * (s1[ii] - cent1) +
                      (s2[ii] - cent2) * (s2[ii] - cent2));
    }

    out = (-phi * sdist);
    if (!lg) out = exp(out);
    return out;
}

// Implementation of the Cohesion6 function
double Cohesion6(double *s1, double *s2, double phi, int dim, int lg) {
    int ii;
    double cent1 = 0, cent2 = 0, sdist = 0, out;
    for (ii = 0; ii < dim; ii++) {
        cent1 += s1[ii] / (double) dim;
        cent2 += s2[ii] / (double) dim;
    }

    for (ii = 0; ii < dim; ii++) {
        sdist += sqrt((s1[ii] - cent1) * (s1[ii] - cent1) +
                      (s2[ii] - cent2) * (s2[ii] - cent2));
    }

    out = -phi * log(sdist);
    if (dim == 1) out = 0.0;
    if (!lg) out = exp(out);
    return out;
}
#include <stdio.h>
#include <math.h>

// Function declarations
double Cohesion1(double *s1, double *s2, double epsilon, int dim, int lg);
double Cohesion2(double *s1, double *s2, double a, int dim, int lg);
double Cohesion5(double *s1, double *s2, double phi, int dim, int lg);
double Cohesion6(double *s1, double *s2, double phi, int dim, int lg);
double Cohesion3_4(double *s1, double *s2, double *mu0, double k0, double v0, double *L0, int dim, int Cohesion, int lg);
double G2a(double a, int lg);

// Cohesion1 function implementation
double Cohesion1(double *s1, double *s2, double epsilon, int dim, int lg) {
    int ii;
    double cent1 = 0, cent2 = 0, sdist = 0, maxdist = 0, dist, out;
    for (ii = 0; ii < dim; ii++) {
        cent1 += s1[ii] / (double) dim;
        cent2 += s2[ii] / (double) dim;
    }
	printf("cent1 %f cent2 %f\n",cent1,cent2);

    maxdist = sqrt((s1[0] - cent1) * (s1[0] - cent1) + (s2[0] - cent2) * (s2[0] - cent2));
    for (ii = 0; ii < dim; ii++) {
        dist = sqrt((s1[ii] - cent1) * (s1[ii] - cent1) + (s2[ii] - cent2) * (s2[ii] - cent2));
        sdist += dist;
        if (maxdist < dist) maxdist = dist;
    }
	printf("sdist %f\n",sdist);

    if (sdist >= 1) {
        out = -lgamma(sdist * epsilon);
		printf("Case if\n");
    } else if (sdist != 0) {
        out = -log(sdist);
		printf("Case elseif\n");
    } else {
        out = log(1);
    }
	printf("out %f\n",out);
    if (!lg) out = exp(out);
    return out;
}

// Cohesion2 function implementation
double Cohesion2(double *s1, double *s2, double a, int dim, int lg) {
    int i, ii;
    double dist, out = 1.0;
    for (i = 0; i < dim; i++) {
        for (ii = 0; ii < dim; ii++) {
            dist = sqrt((s1[i] - s1[ii]) * (s1[i] - s1[ii]) + (s2[i] - s2[ii]) * (s2[i] - s2[ii]));
            if (dist >= a) {
                out = 0.0;
                break;
            }
        }
        if (out == 0.0) break;
    }

    if (lg) out = log(out);
    return out;
}

// Cohesion3_4 function implementation

double Cohesion3_4(double *s1, double *s2, double *mu0, double k0, double v0, double *L0,
					int dim, int Cohesion, int lg){

	int ii;
	double kn,knn,vn,vnn,out,sbar1,sbar2,dL0,dLn,dLnn;
	double s_sbar1,s_sbar2;
	double sbar_mu01, sbar_mu02, sbar_mun1, sbar_mun2, mun1, mun2; 
	// double munn1, munn2;
	double Vs1, Vs2, Vs3, Vs4;
	double Vsbarmu01, Vsbarmu02, Vsbarmu03, Vsbarmu04;
	double Vsbarmun1, Vsbarmun2, Vsbarmun3, Vsbarmun4;
	double Ln1, Ln2, Ln3, Ln4, Lnn1, Lnn2, Lnn3, Lnn4;

	sbar1=0.0, sbar2=0.0;
	for(ii=0; ii<dim; ii++){

		sbar1 = sbar1 + s1[ii]/(double) dim;
		sbar2 = sbar2 + s2[ii]/(double) dim;

	}


	Vs1=0.0, Vs2=0.0, Vs3=0.0, Vs4=0.0;
	for(ii=0; ii<dim; ii++){

		s_sbar1 = s1[ii] - sbar1;
		s_sbar2 = s2[ii] - sbar2;


		Vs1 = Vs1 + s_sbar1*s_sbar1;
		Vs2 = Vs2 + s_sbar1*s_sbar2;
		Vs3 = Vs3 + s_sbar2*s_sbar1;
		Vs4 = Vs4 + s_sbar2*s_sbar2;

	}


	kn = k0 + dim; vn = v0 + dim;
	knn = kn + dim; vnn = vn + dim;

//	Rprintf("kn = %f\n", kn);
//	Rprintf("vn = %f\n", vn);

//	Rprintf("knn = %f\n", knn);
//	Rprintf("vnn = %f\n", vnn);

	mun1 = k0/(k0+dim)*mu0[0] + dim/(k0+dim)*sbar1;
	mun2 = k0/(k0+dim)*mu0[1] + dim/(k0+dim)*sbar2;
	//munn1 = kn/(kn+dim)*mun1 + dim/(kn+dim)*sbar1;
	//munn2 = kn/(kn+dim)*mun2 + dim/(kn+dim)*sbar2;


	sbar_mu01 = sbar1 - mu0[0];
	sbar_mu02 = sbar2 - mu0[1];
	sbar_mun1 = sbar1 - mun1;
	sbar_mun2 = sbar2 - mun2;


	Vsbarmu01 = sbar_mu01*sbar_mu01, Vsbarmun1 = sbar_mun1*sbar_mun1;
	Vsbarmu02 = sbar_mu01*sbar_mu02, Vsbarmun2 = sbar_mun1*sbar_mun2;
	Vsbarmu03 = sbar_mu02*sbar_mu01, Vsbarmun3 = sbar_mun2*sbar_mun1;
	Vsbarmu04 = sbar_mu02*sbar_mu02, Vsbarmun4 = sbar_mun2*sbar_mun2;

	Ln1 = L0[0] + Vs1 + k0*(dim)/(k0+dim)*Vsbarmu01;
	Ln2 = L0[1] + Vs2 + k0*(dim)/(k0+dim)*Vsbarmu02;
	Ln3 = L0[2] + Vs3 + k0*(dim)/(k0+dim)*Vsbarmu03;
	Ln4 = L0[3] + Vs4 + k0*(dim)/(k0+dim)*Vsbarmu04;

	Lnn1 = Ln1 + Vs1 + kn*(dim)/(kn+dim)*Vsbarmun1;
	Lnn2 = Ln2 + Vs2 + kn*(dim)/(kn+dim)*Vsbarmun2;
	Lnn3 = Ln3 + Vs3 + kn*(dim)/(kn+dim)*Vsbarmun3;
	Lnn4 = Ln4 + Vs4 + kn*(dim)/(kn+dim)*Vsbarmun4;


	dL0 = L0[0]*L0[3] - L0[1]*L0[2];
	dLn = Ln1*Ln4 - Ln2*Ln3;
	dLnn = Lnn1*Lnn4 - Lnn2*Lnn3;

//	RprintVecAsMat("L0", L0, 2, 2);
//	RprintVecAsMat("Ln", Ln, 2, 2);
//	RprintVecAsMat("Lnn", Lnn, 2, 2);



	if(Cohesion==3){
		out = -dim*log(M_PI) +
		      (G2a(0.5*vn, 1) - G2a(0.5*v0, 1)) +
		      (0.5*v0*log(dL0) - 0.5*vn*log(dLn)) +
		      (log(k0) - log(kn));
	}
	if(Cohesion==4){
//		Rprintf("G2a(0.5*vn) = %f\n", G2a(0.5*vn,1));
//		Rprintf("G2a(0.5*vnn) = %f\n", G2a(0.5*vnn,1));
		out = -dim*log(M_PI) +
		      (G2a(0.5*vnn,1) - G2a(0.5*vn,1)) +
		      (0.5*vn*log(dLn) - 0.5*vnn*log(dLnn)) +
		      (log(kn) - log(knn));
	}
	printf("out %f\n",out);
//	Rprintf("out = %f\n", out);
	if(!lg) out = exp(out);
	return(out);

}

// G2a function implementation
double G2a(double a, int lg) {
    double out;
    out = log(M_PI) + lgamma(a) + lgamma(a - 0.5);
    if (!lg) out = exp(out);
    return out;
}

int main() {
    // Test data
    // double s1[] = {1.0, 2.0, 3.0};
    // double s2[] = {4.0, 5.0, 6.0};
	// double s1[] = {1.0, 1.0, 1.4};
    // double s2[] = {1.0, 1.1, 1.4};
	double s1[] = {1.0,2};
    double s2[] = {1.1,2};
	
    int dim = 2;
    double phi = 0.5;
    double epsilon = 0.1;
    double a = 5.0;
    double mu0[] = {2.0, 1.0};
    double k0 = 0.5, v0 = 0.1;
    double L0[] = {1.0,0.5, 0.5, 1.0};
    int lg = 1;

    // Test Cohesion1
    double result1 = Cohesion1(s1, s2, epsilon, dim, lg);
    printf("Cohesion1 result: %f\n", result1);

    // Test Cohesion2
    double result2f = Cohesion2(s1, s2, a, dim, 0);
    printf("Cohesion2 result: %f\n", result2f); 
	double result2t = Cohesion2(s1, s2, a, dim, 1);
    printf("Cohesion2 result: %f\n", result2t);


    // Test Cohesion3_4 for Cohesion3
    int Cohesion = 3;
    double result3 = Cohesion3_4(s1, s2, mu0, k0, v0, L0, dim, Cohesion, lg);
    printf("Cohesion3 result: %f\n", result3);

    // Test Cohesion3_4 for Cohesion4
    Cohesion = 4;
    double result4 = Cohesion3_4(s1, s2, mu0, k0, v0, L0, dim, Cohesion, lg);
    printf("Cohesion4 result: %f\n", result4);

	    // Test Cohesion5
    double result5 = Cohesion5(s1, s2, phi, dim, lg);
    printf("Cohesion5 result: %f\n", result5);

    // Test Cohesion6
    double result6 = Cohesion6(s1, s2, phi, dim, lg);
    printf("Cohesion6 result: %f\n", result6);


    return 0;
}
