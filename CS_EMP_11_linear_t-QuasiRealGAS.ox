/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Estimate all linear realized t-quasi GAS model parameters (gamma, omega, alpha and beta)
**
**  Date:
**    	10/01/2015
**
**  Author:
**	  	Tamer Dilaver, Koen de Man & Sina Zolnoor
**
**	Supervisor:
**		L.F. Hoogerheide & S.J. Koopman
**
*/

#include <oxstd.h>
#include <oxdraw.h>
#include <oxprob.h>
#include <maximize.h>
#import <modelbase>
#import <simula>
#include <oxfloat.h>

static decl iB;	 					//Repeats
static decl iSIZE;					//Size of time series
static decl iSTEPS;					//#Steps to divide the size
static decl iSIMS;					//# of Zt ~ N(0,1)
static decl dALPHA;			 		//actually gamma in notes
static decl dBETA;
static decl dOMEGA;
static decl dGAMMA;					//:=h_1
static decl dTAU_1;
static decl dTAU_2;
static decl dXI;
static decl dSIGMA2_U;
static decl dPHI;					//denoted as \varphi in notes
static decl s_vY; 					//Simulated returns
static decl s_vX; 					//Simulated realized measure
static decl s_vDate;
static decl dOMEGA_START;
static decl dBETA_START;
static decl dGAMMA_START;
static decl dXI_START;
static decl dPHI_START;
static decl dTAU_1_START;
static decl dTAU_2_START;
static decl dSIGMA2_U_START;
static decl dLAMBDA_Z_START;
static decl dLAMBDA_U_START;
static decl dRATIO;

/*
**  Function:	Transform (start) parameters
**
**  Input: 		vTheta [parametervalues]
**
**  Output: 	vThetaStar
*/

fTransform(const avThetaStar, const vTheta){
	avThetaStar[0]=		vTheta;
	
	avThetaStar[0][0] = log(vTheta[0]);
	avThetaStar[0][1] = log((vTheta[2]*vTheta[3])/(1-vTheta[2]*vTheta[3] - vTheta[1]));	//aanpassen later		 
	avThetaStar[0][2] = log((vTheta[1])/(1-vTheta[2]*vTheta[3] - vTheta[1]));			//aanpassen later	 
//	avThetaStar[0][3] = log(vTheta[3]);
	avThetaStar[0][3] = log(vTheta[3]);
	avThetaStar[0][4] = vTheta[4];
	avThetaStar[0][5] = vTheta[5];
	avThetaStar[0][6] = log(vTheta[6]);		  
//	avThetaStar[0][8] = log(vTheta[8]-4)-log(100-vTheta[8]);
//	avThetaStar[0][9] = log(vTheta[9]-4)-log(100-vTheta[9]);
	return 1;
}

/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adOmega, adBeta, adGamma, adXi, adPhi, adTau_1, adTau_2, adSigma2_u,  vTheta){
**				0		1	     2	      3	     4	    5	    6	     7		  
**  Output: 	1 
*/
fGetPars(const adOmega, const adBeta, const adGamma, const adPhi, const adTau_1, const adTau_2, const adSigma2_u, const vTheta){

	adOmega[0] 	= exp(vTheta[0]);
	adBeta[0] 	= exp(vTheta[2])/(exp(vTheta[1])+exp(vTheta[2])+1);
	adGamma[0] 	= exp(vTheta[1]-vTheta[3])/(exp(vTheta[1])+exp(vTheta[2])+1);
//	adXi[0] 	= exp(vTheta[3]);
	adPhi[0] 	= exp(vTheta[3]);
	adTau_1[0] 	= vTheta[4];
	adTau_2[0] 	= vTheta[5];
	adSigma2_u[0] = exp(vTheta[6]);
//	adLambda_z[0] = 4+(100-4)*exp(vTheta[8])/(1+exp(vTheta[8]));
//	adLambda_u[0] = 4+(100-4)*exp(vTheta[9])/(1+exp(vTheta[9]));
	return 1;
}


/*
**  Function:	Calculates average value loglikelihood for linear realized t-quasi GAS model given parameter values
**
**  Input: 		vTheta [parameter values], adFunc [adress function value], avScore [the score], amHessian [hessian matrix]
**
**  Output:		1
**
*/
fLogLike_Real_t_GAS(const vTheta, const adFunc, const avScore, const amHessian){
	decl dOmega, dBeta, dGamma, dXi, dPhi, dTau_1, dTau_2, dSigma2_u, dLambda_z, dLambda_u;
	dLambda_u = 4;
	dLambda_z = 4;
	dXi = 0;
	fGetPars(&dOmega, &dBeta, &dGamma, &dPhi, &dTau_1, &dTau_2, &dSigma2_u, vTheta);

	decl dH = (dOmega +dGamma*dXi)/(1-dBeta-dPhi*dGamma); //initialise with unconditional expectation of log conditional variance					 											//initial condition by definition
	decl vEta = zeros(sizerc(s_vY), 1);

	for(decl i = 0; i < sizerc(s_vY); ++i){
			//likelihood contribution
			vEta[i] = -log(M_PI)-1/2*log(dLambda_z-2) -1/2*log(dH) -log(gammafact(dLambda_z/2))+log(gammafact((dLambda_z+1)/2)) -(dLambda_z+1)/2*log(1+ s_vY[i]^2 / ((dLambda_z-2)*dH)) -1/2*log(dLambda_u-2) -1/2*log(dSigma2_u) -log(gammafact(dLambda_u/2))+log(gammafact((dLambda_u+1)/2)) -(dLambda_u+1)/2*log(1+ (s_vX[i] - dXi - dPhi*dH - dTau_1*s_vY[i]/sqrt(dH) - dTau_2*(s_vY[i]^2/dH-1))^2 / ((dLambda_u-2)*dSigma2_u));		//Standardized Student's t

			//recursion
			  dH = dOmega + dGamma*(dLambda_u+3)/dLambda_u*((dLambda_u+1)/(dLambda_u-2)*(1+s_vX[i]/((dLambda_u-2)*  dSigma2_u))^(-1)*s_vX[i] -   dSigma2_u) + dBeta*dH;
	}
	adFunc[0] = sumc(vEta)/sizerc(s_vY); 									 	//Average
	return 1;
}

/*
**  Function:	Transform parameters back
**
**  Input: 		vThetaStar
**
**  Output: 	vTheta [parametervalues]
*/
fTransformBack(const avTheta, const vThetaStar){
	avTheta[0]=		vThetaStar;
	
	avTheta[0][0] = exp(vThetaStar[0]);
	avTheta[0][1] = exp(vThetaStar[2])/(exp(vThetaStar[1])+exp(vThetaStar[2])+1);
	avTheta[0][2] = exp(vThetaStar[1]-vThetaStar[3])/(exp(vThetaStar[1])+exp(vThetaStar[2])+1);
//	avTheta[0][3] = exp(vThetaStar[3]);
	avTheta[0][3] = exp(vThetaStar[3]);
	avTheta[0][4] = vThetaStar[4];
	avTheta[0][5] = vThetaStar[5];
	avTheta[0][6] = exp(vThetaStar[6]);
//	avTheta[0][8] = 4+(100-4)*exp(vThetaStar[8])/(1+exp(vThetaStar[8]));
//	avTheta[0][9] = 4+(100-4)*exp(vThetaStar[9])/(1+exp(vThetaStar[9]));
	return 1;
}

/*
**  Function:	calculate standard errors
**
**  Input: 		vThetaStar
**
**  Output: 	vStdErrors
*/
fSigmaStdError(const vThetaStar){

 		decl iN, mHessian, mHess, mJacobian, vStdErrors, vP;

		iN 			= sizerc(s_vY);
		Num2Derivative(fLogLike_Real_t_GAS, vThetaStar, &mHessian);
		NumJacobian(fTransformBack, vThetaStar, &mJacobian);	  //numerical Jacobian
		mHessian 	= mJacobian*invert(-iN*mHessian)*mJacobian';
		vStdErrors 	= sqrt(diagonal(mHessian)');

		return 	vStdErrors;
}

/*
**  Function:	calculate variance of model
**
**  Input: 		vTheta
**
**  Output: 	vH [vector with variances]
*/

fVariance(const vTheta){
	decl dOmega, dBeta, dGamma, dXi, dPhi, dTau_1, dTau_2, dSigma2_u, dLambda_z, dLambda_u, vH;
	dLambda_u = 4;
	dLambda_z = 4;
	dXi = 0;
	fGetPars(&dOmega, &dBeta, &dGamma, &dPhi, &dTau_1, &dTau_2, &dSigma2_u, vTheta);

	
	vH = zeros(sizerc(s_vY),1);
	vH[0] = (dOmega +dGamma*dXi)/(1-dBeta-dPhi*dGamma);	
	for(decl i = 1; i < sizerc(s_vY); i++){													   //mixed 	
		 	vH[i] = dOmega + dGamma*(dLambda_u+3)/dLambda_u*((dLambda_u+1)/(dLambda_u-2)*(1+s_vX[i-1]/((dLambda_u-2)*dSigma2_u))^(-1)*s_vX[i-1] - dSigma2_u) + dBeta*vH[i-1];
	}		   
	
	return 	vH;


}

/*
**  Function:	Estimate linear realized t-quasi GAS model parameters
**
**  Input: 		 vReturns, vRealMeasure, adOmega_hat, adBeta_hat, adGamma_hat, adPhi_hat, adTau_1_hat, adTau_2_hat, adSigma2_u_hat, avVariance
**
**  Output: 	vTheta [estimated parametervalues]
*/

fEstimateReal_t_GAS(const vReturns, const vRealMeasure, const adOmega_hat, const adBeta_hat, const adGamma_hat, const adPhi_hat, const adTau_1_hat, const adTau_2_hat, const adSigma2_u_hat, const avVariance){

	//initialise parameter values
	decl vTheta = zeros(7,1);
	vTheta[0] = dOMEGA_START;
	vTheta[1] = dBETA_START;
	vTheta[2] = dGAMMA_START;
//	vTheta[3] = dXI_START;
	vTheta[3] = dPHI_START;
	vTheta[4] = dTAU_1_START;
	vTheta[5] = dTAU_2_START;
	vTheta[6] = dSIGMA2_U_START;
//	vTheta[8] = dLAMBDA_Z_START;
//	vTheta[9] = dLAMBDA_U_START;
//	vTheta = <-2.7129 ; 0.70111 ; 0.29256 ; 9.0774 ; 0.99598 ; -0.0066165 ; 0.099432 ; 0.21021 ; 7 ; 7>; //note we changed starting value because no improvement in linesearch
																										 //we picked optimal values Gaussian version
	decl vThetaStart = vTheta;

	//globalalize returns and vectorize true pars
	s_vY = vReturns;
	s_vX = vRealMeasure;

	//transform parameters
	decl vThetaStar; 
	fTransform(&vThetaStar, vTheta);

	//Maximize the LL
	decl dFunc;
	decl iA;
	iA=MaxBFGS(fLogLike_Real_t_GAS, &vThetaStar, &dFunc, 0, TRUE);

	//Transform thetasStar back
  	fTransformBack(&vTheta, vThetaStar);

	//return pars
	adOmega_hat[0] = vTheta[0];
	adBeta_hat[0] = vTheta[1];	
	adGamma_hat[0] = vTheta[2];
//	adXi_hat[0] = vTheta[3];
	adPhi_hat[0] = vTheta[3];
	adTau_1_hat[0] = vTheta[4];
	adTau_2_hat[0] = vTheta[5];
	adSigma2_u_hat[0] = vTheta[6];
//	adLambda_z_hat[0] = vTheta[8];
//	adLambda_u_hat[0] = vTheta[9];

	decl vSigmaStdError = fSigmaStdError(vThetaStar);
	decl vVariance = fVariance(vThetaStar);
	avVariance[0] = vVariance;
	
	print("\n",MaxConvergenceMsg(iA));
	println("\nFunctiewaarde likelihood eindwaardes:", dFunc);
	print("\nOptimale parameters met standaarderrors \n",
          	"%r", { "omega",  "beta",  "gamma",  "xi=0", "phi","tau_1","tau_2","sigma2_u","lambda_z=4","lambda_u =4"},
          	"%c", {"thetaStart","theta","std.error"}, (vThetaStart[0:2]'~0~vThetaStart[3:6]'~4~4)'~(vThetaStart[0:2]'~0~vThetaStart[3:6]'~4~4)'~(vSigmaStdError[0:2]'~0~vSigmaStdError[3:6]'~0~0)');
			
	return 1;
}

/*
**  Function:	Determine Forecast
**
**  Input: 		vTheta
**
**  Output: 	vH [vector of forecasts]
*/

fForecast(const vTheta){
	decl dOmega, dBeta, dGamma, dXi, dPhi, dTau_1, dTau_2, dSigma2_u, dLambda_z, dLambda_u, vH;
	fGetPars(&dOmega, &dBeta, &dGamma, &dPhi, &dTau_1, &dTau_2, &dSigma2_u, vTheta);
	dLambda_u =4;
	dLambda_u =4;
	dXi = 0;
	
	vH = zeros((sizerc(s_vY)+1),1);
	vH[0] = (dOmega +dGamma*dXi)/(1-dBeta-dPhi*dGamma);	
	
	for(decl i = 1; i < sizerc(s_vY)+1; i++){								
		vH[i] = dOmega + dGamma*(dLambda_u+3)/dLambda_u*((dLambda_u+1)/(dLambda_u-2)*(1+s_vX[i-1]/((dLambda_u-2)*dSigma2_u))^(-1)*s_vX[i-1] - dSigma2_u) + dBeta*vH[i-1];
	}		   
	
	return 	vH[sizerc(s_vY)];
}

/*
**  Function:	Compute MAE	and MSE
**
**  Input: 		adMAE_OC [adress of MAE], adMSE_OC [adress of MSE] vReturns_1 [return series], vBenchmark [Benchmark], dC [ratio]
**
**  Output: 	1
*/

fMAE(const adMAE, const adMSE, const vReturns, const vRV, const vBenchmark, const dC){

	decl iWindow = 250;
	decl iT = sizerc(vReturns);
//	decl vTemp_returns = vReturns_1;
	decl vH_forecast = zeros(iWindow, 1);
	decl vError = zeros(iWindow, 1);

	dOMEGA_START = 0.098956;
	dBETA_START = 0.64342;
	dGAMMA_START = 0.45721;
//	dXI_START = 0.12208;
	dPHI_START = 0.72717;
	dTAU_1_START = -0.034277;
	dTAU_2_START = 0.28478 ;
	dSIGMA2_U_START = 0.2122;
//	dLAMBDA_Z_START = 10;
//	dLAMBDA_U_START = 4
	
	for(decl j = 0; j<iWindow; j++){
		s_vY = 	vReturns[j:(iT - iWindow +j)];
		s_vX = 	vRV[j:(iT - iWindow +j)];

		//initialise parametervalues
		decl vTheta = zeros(7,1);
		vTheta[0] = dOMEGA_START;
		vTheta[1] = dBETA_START;
		vTheta[2] = dGAMMA_START;
	//	vTheta[3] = dXI_START;
		vTheta[3] = dPHI_START;
		vTheta[4] = dTAU_1_START;
		vTheta[5] = dTAU_2_START;
		vTheta[6] = dSIGMA2_U_START;
	//	vTheta[8] = dLAMBDA_Z_START;
	//	vTheta[9] = dLAMBDA_U_START;
	
		//transform parameters
		decl vThetaStar; 
		fTransform(&vThetaStar, vTheta);
	
		//Maximize the LL
		decl dFunc;
		decl iA;
		iA=MaxBFGS(fLogLike_Real_t_GAS, &vThetaStar, &dFunc, 0, TRUE);
	
		//Transform thetasStar back
	  	fTransformBack(&vTheta, vThetaStar);

		dOMEGA_START = vTheta[0];
		dBETA_START = vTheta[1];
		dGAMMA_START = vTheta[2];
	//	dXI_START = vTheta[3];
		dPHI_START = vTheta[3];
		dTAU_1_START = vTheta[4];
		dTAU_2_START = vTheta[5];
		dSIGMA2_U_START = vTheta[6];
	//	dLAMBDA_Z_START = vTheta[8];
	//	dLAMBDA_U_START = 99.9999;
		vH_forecast[j] = fForecast(vThetaStar);
		vError[j] = (dC*vH_forecast[j] - dRATIO*vBenchmark[(iT - iWindow +j)]);

	}

	if(dC==1){
		savemat("vAE_CC_RK_LINEAR_RK_t-QuasiREALGAS.xls", fabs(vError));

		DrawTMatrix(0, (dC*vH_forecast~dRATIO*vBenchmark[iT-iWindow:])', {"Close-to-close"}, s_vDate[iT-iWindow:]');
//		ShowDrawWindow();
	}else{
		savemat("vAE_OC_RK_LINEAR_RK_t-QuasiREALGAS.xls", fabs(vError));
//		SetDrawWindow("CS_EMP_10_FORECASTS");
		DrawTMatrix(1, (dC*vH_forecast~dRATIO*vBenchmark[iT-iWindow:])', {"Open-to-close"}, s_vDate[iT-iWindow:]');
//		ShowDrawWindow();
	}
	
	adMAE[0] = meanc(fabs(vError));
	adMSE[0] = meanc(sqr(vError));
	return 1;

}

/*
**				MAIN PROGRAM
**
**  Purpose:	Estimate linear realized t-quasi GAS model parameters 
**
**  Output: 	Figures
*/

main(){
	//laad SBUX returns
	decl mData_1 = loadmat("ReturnsOpenToClose.csv");
	decl mData_2 = loadmat("ReturnsCloseToClose.csv"); 
	decl vReturns_1 = 100*mData_1[:][1];
	decl vReturns_2 = 100*mData_2[:][1];

	decl vRV = loadmat("RV.csv");
	decl vBV = loadmat("BV.csv");
	decl vRK = loadmat("RK.csv");

	dRATIO = (varc(vReturns_1) +varc(vReturns_2))/varc(vReturns_1);
	
	s_vX = vRK;   		//pick vRV, vBV or vKV

	//laad Dates SBUX returns
	decl vTemp_Date = mData_2[][0];
	decl vYear 		= floor(vTemp_Date/10000);							
	decl vMonth 	= floor((vTemp_Date-floor(vTemp_Date/10000)*10000)/100);	
	decl vDay 		= vTemp_Date-floor(vTemp_Date/100)*100;
	s_vDate 		= dayofcalendar(vYear, vMonth, vDay);

	dOMEGA_START = 2.1514;
	dBETA_START = 0.50939;
	dGAMMA_START = 0.38231;
	dPHI_START = 0.56632;
	dTAU_1_START = -0.10168;
	dTAU_2_START = 0.34765 ;
	dSIGMA2_U_START = 3.9044;


//                thetaStart        theta    std.error
//omega             0.098956       2.1514         .NaN
//beta               0.64342      0.50939         .NaN
//gamma              0.45721      0.38231         .NaN
//xi                 0.12208  2.1649e-225         .NaN
//phi                0.72717      0.56632         .NaN
//tau_1            -0.034277     -0.10168         .NaN
//tau_2              0.28478      0.34765         .NaN
//sigma2_u           0.21220       3.9044         .NaN
//lambda_z            4.0000       4.0000       0.0000
//lambda_u =4         4.0000       4.0000       0.0000

	decl dOmega_hat, dBeta_hat, dGamma_hat, dXi_hat, dPhi_hat, dTau_1_hat, dTau_2_hat, dSigma2_u_hat, dLambda_z;
	decl vVariance_1, vVariance_2;
	print("\nO-C");
	fEstimateReal_t_GAS(vReturns_1, s_vX, &dOmega_hat, &dBeta_hat, &dGamma_hat, &dPhi_hat, &dTau_1_hat, &dTau_2_hat, &dSigma2_u_hat, &vVariance_1);

	print("\nC-C");
	fEstimateReal_t_GAS(vReturns_2, s_vX, &dOmega_hat, &dBeta_hat, &dGamma_hat, &dPhi_hat, &dTau_1_hat, &dTau_2_hat, &dSigma2_u_hat, &vVariance_2);

	//graphs
	SetDrawWindow("CS_EMP_11_linear_t-QuasiRealGAS(1,1)");
	DrawTMatrix(0, (vReturns_1~sqrt(vVariance_1))', {"Open-to-close"}, s_vDate');
	DrawTMatrix(1, (vReturns_2~sqrt(vVariance_2))', {"Close-to-close"}, s_vDate');
	ShowDrawWindow();

	//forecasts MAE	 and MSE
	decl vBenchmark = vRV;
	decl dMAE_OC;
	decl dMSE_OC;
	SetDrawWindow("CS_EMP_10_FORECASTS");
	fMAE(&dMAE_OC, &dMSE_OC, vReturns_1, vRV, vBenchmark, dRATIO);
	print("\n dMAE_OC = ",dMAE_OC);
	print("\n dMSE_OC = ",dMSE_OC);
	
	decl dMAE_CC;
	decl dMSE_CC;
	fMAE(&dMAE_CC, &dMSE_CC, vReturns_2, vRV, vBenchmark, 1);
	print("\n dMAE_CC = ",dMAE_CC);
	print("\n dMSE_CC = ",dMSE_CC);
	ShowDrawWindow();
}