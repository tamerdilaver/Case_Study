/*
**	Case Study Financial Econometrics 4.3 		
**
**  Purpose:
**  	Estimate all Gaussian log-linear RealGARCH(1,1) model parameters (alpha, gamma, omega, xi, phi, tau_1, tau_2, sigma2_u and beta)
**
**  Date:
**    	17/01/2015
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
static decl	iSTD_ERROR;				//0 or 1 binary

/*
**  Function:	Simulate Gaussian log-linear RealGARCH(1,1) returns for given parameters
**
**  Input:		dOmega, dBeta, dGamma, dXi, dPhi, dTau_1, dTau_2, dSigma2_u, avReturns, avRealMeasure [this is the realized measure]]
**
**  Output:		1
**
*/

fSimLogRealGARCH(const dOmega, const dBeta, const dGamma, const dXi, const dPhi, const dTau_1, const dTau_2, const dSigma2_u, const avReturns, const avRealMeasure){
	decl vTemp, vLogH, vLogX, dZ_Temp, dU_Temp;
	vTemp = vLogH = vLogX =  zeros(iSIZE+1, 1);

	vLogH[0] = (dOmega +dGamma*dXi)/(1-dBeta-dPhi*dGamma);		//by definition
	
	for(decl i = 0; i < iSIZE; i++){
		dZ_Temp = rann(1,1);
		dU_Temp = rann(1,1);
		vTemp[i] = sqrt(exp(vLogH[i]))*dZ_Temp;
		vLogX[i] = dXi + dPhi*vLogH[i] + dTau_1*dZ_Temp + dTau_2*(dZ_Temp^2-1)+ sqrt(dSigma2_u)*dU_Temp;
		vLogH[i+1] = dOmega  + dBeta*vLogH[i]+ dGamma*vLogX[i];
	}

	vTemp = dropr(vTemp,iSIZE);
	vLogX = dropr(vLogX,iSIZE);
	
	avReturns[0] = vTemp;
	avRealMeasure[0] = exp(vLogX);
	return 1;
}

/*
**  Function:	Transform (start)parameters
**
**  Input: 		vTheta [parametervalues]
**
**  Output: 	vThetaStar
*/

fTransform(const avThetaStar, const vTheta){
	avThetaStar[0]=		vTheta;
	
	avThetaStar[0][0] = vTheta[0];
	avThetaStar[0][1] = log((vTheta[2]*vTheta[4])/(1-vTheta[2]*vTheta[4] - vTheta[1]));		//aanpassen later	 
	avThetaStar[0][2] = log((vTheta[1])/(1-vTheta[2]*vTheta[4] - vTheta[1]));				//aanpassen later	 	 
	avThetaStar[0][3] = vTheta[3];
	avThetaStar[0][4] = log(vTheta[4]);
	avThetaStar[0][5] = vTheta[5];
	avThetaStar[0][6] = vTheta[6];
	avThetaStar[0][7] = log(vTheta[7]);
	return 1;
}

/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adOmega, adBeta, adGamma, adXi, adPhi, adTau_1, adTau_2, adSigma2_u,  vTheta
**				0		1	     2	      3	     4	    5	    6	     7		  
**  Output: 	1 
*/

fGetPars(const adOmega, const adBeta, const adGamma, const adXi, const adPhi, const adTau_1, const adTau_2, const adSigma2_u, const vTheta){

	adOmega[0] 	= vTheta[0];
	adBeta[0] 	= exp(vTheta[2])/(exp(vTheta[1])+exp(vTheta[2])+1);
	adGamma[0] 	= exp(vTheta[1]-vTheta[4])/(exp(vTheta[1])+exp(vTheta[2])+1);
	adXi[0] 	= vTheta[3];
	adPhi[0] 	= exp(vTheta[4]);
	adTau_1[0] 	= vTheta[5];
	adTau_2[0] 	= vTheta[6];
	adSigma2_u[0] = exp(vTheta[7]);
	return 1;
}

/*
**  Function:	Calculates average value loglikelihood for log-linear Realized GARCH given parameter values
**
**  Input: 		vTheta [parametervalues], adFunc [adres functievalue], avScore [the score], amHessian [hessianmatrix]
**
**  Output:		1
**
*/

fLogLike_LogRealGARCH(const vTheta, const adFunc, const avScore, const amHessian){
	decl dOmega, dBeta, dGamma, dXi, dPhi, dTau_1, dTau_2, dSigma2_u;
	fGetPars(&dOmega, &dBeta, &dGamma, &dXi, &dPhi, &dTau_1, &dTau_2, &dSigma2_u, vTheta);

	decl dlogH = (dOmega +dGamma*dXi)/(1-dBeta-dPhi*dGamma); //initialise with unconditional expectation of log conditional variance					 											//initial condition by definition

	decl vlogEta = zeros(sizerc(s_vY), 1);

	for(decl i = 0; i < sizerc(s_vY); ++i){
			//likelihood contribution
			vlogEta[i] = 2*log(M_2PI) +dlogH + s_vY[i]^2 / exp(dlogH) + log(dSigma2_u) + (log(s_vX[i]) - dXi - dPhi*dlogH - dTau_1*s_vY[i]/sqrt(exp(dlogH)) - dTau_2*(s_vY[i]^2/exp(dlogH)-1))^2/dSigma2_u;		//Gaussian

			//log-linear Realized GARCH recursion
			dlogH = dOmega + dBeta * dlogH + dGamma * log(s_vX[i]);
	}
	
	adFunc[0] = sumc(vlogEta)/(-2*sizerc(s_vY)); //Average
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
	
	avTheta[0][0] = vThetaStar[0];
	avTheta[0][1] = exp(vThetaStar[2])/(exp(vThetaStar[1])+exp(vThetaStar[2])+1);
	avTheta[0][2] = exp(vThetaStar[1]-vThetaStar[4])/(exp(vThetaStar[1])+exp(vThetaStar[2])+1);
	avTheta[0][3] = vThetaStar[3];
	avTheta[0][4] = exp(vThetaStar[4]);
	avTheta[0][5] = vThetaStar[5];
	avTheta[0][6] = vThetaStar[6];
	avTheta[0][7] = exp(vThetaStar[7]);
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
		Num2Derivative(fLogLike_LogRealGARCH, vThetaStar, &mHessian);
		NumJacobian(fTransformBack, vThetaStar, &mJacobian);	  //numerical Jacobian
		mHessian 	= mJacobian*invert(-iN*mHessian)*mJacobian';
		vStdErrors 	= sqrt(diagonal(mHessian)');

		return 	vStdErrors;
}

/*
**  Function:	Estimate log-linear Realized GARCH parameters
**
**  Input: 		vReturns, adAlpha_hat, adBeta_hat, adOmega_hat, adGamma_hat, adXi_hat, adPhi_hat, adTau_1_hat, adTau_2_hat, adSigma2_u_hat
**
**  Output: 	vTheta [estimated parametervalues]
*/

fEstimateLogRealGARCH(const vReturns, const vRealMeasure, const adOmega_hat, const adBeta_hat, const adGamma_hat, const adXi_hat, const adPhi_hat, const adTau_1_hat, const adTau_2_hat, const adSigma2_u_hat){

	//initialise parametervalues
	decl vTheta = zeros(8,1);
	vTheta = <0.2 ; 0.6 ; 0.4 ; -0.3 ; 0.8 ; -0.05 ; 0.10 ; 0.2>;
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
	iA=MaxBFGS(fLogLike_LogRealGARCH, &vThetaStar, &dFunc, 0, TRUE);

	//Transform thetasStar back
  	fTransformBack(&vTheta, vThetaStar);

	//return pars
	adOmega_hat[0] = vTheta[0];
	adBeta_hat[0] = vTheta[1];	
	adGamma_hat[0] = vTheta[2];
	adXi_hat[0] = vTheta[3];
	adPhi_hat[0] = vTheta[4];
	adTau_1_hat[0] = vTheta[5];
	adTau_2_hat[0] = vTheta[6];
	adSigma2_u_hat[0] = vTheta[7];

	if(iSTD_ERROR){		//only do this for fMonteCarlo2
		decl vSigmaStdError = fSigmaStdError(vThetaStar);
		return vSigmaStdError;
	}else{
		return 1;
	}
}

/*
**  Function:	Simulates and Estimates log-linear Realized data and parameters many times
**				to illustrate Asymptotic normality
**
**  Input: 		amMonteCarlo [matrix of many estimated parameters];
**
**  Output: 	1
*/

fMonteCarlo(const amMonteCarlo){
	decl mTemp;
	mTemp = zeros(iB,8);

	for(decl i = 0; i<iB ; i++){
		decl vReturns, vRealMeasure;
		fSimLogRealGARCH(dOMEGA, dBETA, dGAMMA, dXI, dPHI, dTAU_1, dTAU_2, dSIGMA2_U, &vReturns, &vRealMeasure);

		decl dOmega_hat, dBeta_hat, dGamma_hat, dXi_hat, dPhi_hat, dTau_1_hat, dTau_2_hat, dSigma2_u_hat;
		fEstimateLogRealGARCH(vReturns, vRealMeasure, &dOmega_hat, &dBeta_hat, &dGamma_hat, &dXi_hat, &dPhi_hat, &dTau_1_hat, &dTau_2_hat, &dSigma2_u_hat);	

		mTemp[i][0] =  dOmega_hat;			
		mTemp[i][1]	=  dBeta_hat;
		mTemp[i][2]	=  dGamma_hat;
		mTemp[i][3]	=  dXi_hat;
		mTemp[i][4]	=  dPhi_hat;
		mTemp[i][5]	=  dTau_1_hat;
		mTemp[i][6]	=  dTau_2_hat;
		mTemp[i][7]	=  dSigma2_u_hat;
	}
	amMonteCarlo[0] = mTemp;
	return 1;
}

/*
**  Function:	Simulated and Estimates log linear realized GARCH data and parameters many times
**				to illustrate consistency it returns minimum, mean and maximum values for the estimated parameters
**
**  Input: 		amOmega [matrix containing the min, max and mean of estimated omega],...
**
**  Output: 	1
*/

fMonteCarlo2(const amOmega, const amBeta, const amGamma, const amXi, const amPhi, const amTau_1, const amTau_2, const amSigma2_u, const amOmega2, const amBeta2, const amGamma2, const amXi2, const amPhi2, const amTau_12, const amTau_22, const amSigma2_u2){
	decl mTemp, mTempOmega, mTempBeta, mTempGamma, mTempXi, mTempPhi, mTempTau_1, mTempTau_2, mTempSigma2_u;
	decl mTemp2, mTempOmega2, mTempBeta2, mTempGamma2, mTempXi2, mTempPhi2, mTempTau_12, mTempTau_22, mTempSigma2_u2;
	mTempOmega = mTempBeta = mTempGamma = mTempXi = mTempPhi = mTempTau_1 = mTempTau_2 = mTempSigma2_u = zeros((iSIZE/iSTEPS),3);
	mTempOmega2 = mTempBeta2 = mTempGamma2 = mTempXi2 = mTempPhi2 = mTempTau_12 = mTempTau_22 = mTempSigma2_u2 = zeros((iSIZE/iSTEPS),3);
	mTemp = mTemp2 = zeros(iB,8);

	decl iSize = iSIZE;

	for(decl j = 0; j<(iSize/iSTEPS) ; j++){
		iSIZE = ((iSTEPS)*(j+1));
		for(decl i = 0; i<iB ; i++){
			decl vReturns, vRealMeasure;
			fSimLogRealGARCH(dOMEGA, dBETA, dGAMMA, dXI, dPHI, dTAU_1, dTAU_2, dSIGMA2_U, &vReturns, &vRealMeasure);
	
			decl dOmega_hat, dBeta_hat, dGamma_hat, dXi_hat, dPhi_hat, dTau_1_hat, dTau_2_hat, dSigma2_u_hat, vSE;
			vSE = fEstimateLogRealGARCH(vReturns, vRealMeasure, &dOmega_hat, &dBeta_hat, &dGamma_hat, &dXi_hat, &dPhi_hat, &dTau_1_hat, &dTau_2_hat, &dSigma2_u_hat);	

			mTemp[i][0] =  sqrt(iSIZE)*(dOmega_hat-dOMEGA);			
			mTemp[i][1]	=  sqrt(iSIZE)*(dBeta_hat-dBETA);
			mTemp[i][2]	=  sqrt(iSIZE)*(dGamma_hat-dGAMMA);
			mTemp[i][3]	=  sqrt(iSIZE)*(dXi_hat-dXI);
			mTemp[i][4]	=  sqrt(iSIZE)*(dPhi_hat-dPHI);
			mTemp[i][5]	=  sqrt(iSIZE)*(dTau_1_hat - dTAU_1);
			mTemp[i][6]	=  sqrt(iSIZE)*(dTau_2_hat-dTAU_2);
			mTemp[i][7]	=  sqrt(iSIZE)*(dSigma2_u_hat-dSIGMA2_U);

			mTemp2[i][0] 	=  (dOmega_hat-dOMEGA)/vSE[0];		
			mTemp2[i][1]	=  (dBeta_hat-dBETA)/vSE[1];
			mTemp2[i][2]	=  (dGamma_hat-dGAMMA)/vSE[2];
			mTemp2[i][3]	=  (dXi_hat-dXI)/vSE[3];
			mTemp2[i][4]	=  (dPhi_hat-dPHI)/vSE[4];
			mTemp2[i][5]	=  (dTau_1_hat - dTAU_1)/vSE[5];
			mTemp2[i][6]	=  (dTau_2_hat-dTAU_2)/vSE[6];
			mTemp2[i][7]	=  (dSigma2_u_hat-dSIGMA2_U)/vSE[7];
		} 

		// v0.025_quantile, vMean, v0.975_quantile;				We get 95%-intervals
		mTempOmega[j][0] = quantilec(mTemp[][],0.025)'[0];
		mTempOmega[j][1] = meanc(mTemp[][])'[0];
		mTempOmega[j][2] = quantilec(mTemp[][],0.975)'[0];
		
		mTempBeta[j][0] = quantilec(mTemp[][],0.025)'[1];
		mTempBeta[j][1] = meanc(mTemp[][])'[1];
		mTempBeta[j][2] = quantilec(mTemp[][],0.975)'[1];

		mTempGamma[j][0] = quantilec(mTemp[][],0.025)'[2];
		mTempGamma[j][1] = meanc(mTemp[][])'[2];
		mTempGamma[j][2] = quantilec(mTemp[][],0.975)'[2];

		mTempXi[j][0] = quantilec(mTemp[][],0.025)'[3];
		mTempXi[j][1] = meanc(mTemp[][])'[3];
		mTempXi[j][2] = quantilec(mTemp[][],0.975)'[3];
		
		mTempPhi[j][0] = quantilec(mTemp[][],0.025)'[4];
		mTempPhi[j][1] = meanc(mTemp[][])'[4];
		mTempPhi[j][2] = quantilec(mTemp[][],0.975)'[4];

		mTempTau_1[j][0] = quantilec(mTemp[][],0.025)'[5];
		mTempTau_1[j][1] = meanc(mTemp[][])'[5];
		mTempTau_1[j][2] = quantilec(mTemp[][],0.975)'[5];

		mTempTau_2[j][0] = quantilec(mTemp[][],0.025)'[6];
		mTempTau_2[j][1] = meanc(mTemp[][])'[6];
		mTempTau_2[j][2] = quantilec(mTemp[][],0.975)'[6];

		mTempSigma2_u[j][0] = quantilec(mTemp[][],0.025)'[7];
		mTempSigma2_u[j][1] = meanc(mTemp[][])'[7];
		mTempSigma2_u[j][2] = quantilec(mTemp[][],0.975)'[7];

		// v0.025_quantile, v0.5_quantile, v0.975_quantile;				We get 95%-intervals
		mTempOmega2[j][0] = quantilec(mTemp2[][],0.025)'[0];
		mTempOmega2[j][1] = quantilec(mTemp2[][],0.5)'[0];
		mTempOmega2[j][2] = quantilec(mTemp2[][],0.975)'[0];
		
		mTempBeta2[j][0] = quantilec(mTemp2[][],0.025)'[1];
		mTempBeta2[j][1] = quantilec(mTemp2[][],0.5)'[1];
		mTempBeta2[j][2] = quantilec(mTemp2[][],0.975)'[1];

		mTempGamma2[j][0] = quantilec(mTemp2[][],0.025)'[2];
		mTempGamma2[j][1] = quantilec(mTemp2[][],0.5)'[2];
		mTempGamma2[j][2] = quantilec(mTemp2[][],0.975)'[2];

		mTempXi2[j][0] = quantilec(mTemp2[][],0.025)'[3];
		mTempXi2[j][1] = quantilec(mTemp2[][],0.5)'[3];
		mTempXi2[j][2] = quantilec(mTemp2[][],0.975)'[3];
		
		mTempPhi2[j][0] = quantilec(mTemp2[][],0.025)'[4];
		mTempPhi2[j][1] = quantilec(mTemp2[][],0.5)'[4];
		mTempPhi2[j][2] = quantilec(mTemp2[][],0.975)'[4];

		mTempTau_12[j][0] = quantilec(mTemp2[][],0.025)'[5];
		mTempTau_12[j][1] = quantilec(mTemp2[][],0.5)'[5];;
		mTempTau_12[j][2] = quantilec(mTemp2[][],0.975)'[5];

		mTempTau_22[j][0] = quantilec(mTemp2[][],0.025)'[6];
		mTempTau_22[j][1] = quantilec(mTemp2[][],0.5)'[6];
		mTempTau_22[j][2] = quantilec(mTemp2[][],0.975)'[6];

		mTempSigma2_u2[j][0] = quantilec(mTemp2[][],0.025)'[7];
		mTempSigma2_u2[j][1] = quantilec(mTemp2[][],0.5)'[7];
		mTempSigma2_u2[j][2] = quantilec(mTemp2[][],0.975)'[7];
	}

	amOmega[0] = mTempOmega;
	amBeta[0] = mTempBeta;
	amGamma[0] = mTempGamma;
	amXi[0] = mTempXi;
	amPhi[0] = mTempPhi;
	amTau_1[0] = mTempTau_1;
	amTau_2[0] = mTempTau_2;
	amSigma2_u[0] = mTempSigma2_u;

	amOmega2[0] = mTempOmega2;
	amBeta2[0] = mTempBeta2;
	amGamma2[0] = mTempGamma2;
	amXi2[0] = mTempXi2;
	amPhi2[0] = mTempPhi2;
	amTau_12[0] = mTempTau_12;
	amTau_22[0] = mTempTau_22;
	amSigma2_u2[0] = mTempSigma2_u2;

	return 1;
}

/*
**				MAIN PROGRAM
**
**  Purpose:	Simulate Log-linear Realized GARCH returns 
**				Estimate Log-linear Realized GARCH parameters 
**
**  Input: 		dOMEGA, dBETA, dGAMMA, dXI, dPHI, dTAU_1, dTAU_2, dSIGMA2_U, iB, iSIZE, iSIMS, iSTEPS
**
**  Output: 	Figures
*/

main()
{
	//we don't simulate the startparameter here (too much computation)
	//we don't simulate the mean here (too much computation)
	//SET PARAMETERS
	dOMEGA = 0.2;
	dBETA = 0.6;
	dGAMMA = 0.4;
	dXI = -0.3;
	dPHI = 0.8;
	dTAU_1 = -0.05;
	dTAU_2 = 0.10;
	dSIGMA2_U = 0.2; 
	
/*
** ..................................................................................	
**	 		ASYMPTOTIC NORMALITY
**..................................................................................
*/

	//SET # OF SIMULATIONS 
	iB = 500; 			//max 5000
	iSIZE = 5000;		//max 5000
	iSIMS = iB*iSIZE;
	iSTD_ERROR = 0;

	//DO MANY SIMULATIONS AND ESITMATIONS	
	decl mMonteCarlo;
	fMonteCarlo(&mMonteCarlo);	  

	//DRAW GRAPHS
	SetDrawWindow("CS_SIM_3_asymptotic normality");
	DrawDensity(0, (mMonteCarlo[][0])', {"(i) Density omega"});
	DrawDensity(1, (mMonteCarlo[][1])', {"(ii) Density beta"});
	DrawDensity(2, (mMonteCarlo[][2])', {"(iii) Density gamma"});
	DrawDensity(3, (mMonteCarlo[][3])', {"(iv) Density xi"});
	DrawDensity(4, (mMonteCarlo[][4])', {"(v) Density phi"});
	DrawDensity(5, (mMonteCarlo[][5])', {"(vi) Density tau_1"});
	DrawDensity(6, (mMonteCarlo[][6])', {"(vii) Density tau_2"});
	DrawDensity(7, (mMonteCarlo[][7])', {"(viii) Density sigma2_u"});
	ShowDrawWindow();

	print("\nFirst Graph Finished at ",time(),"\n");
	
/*
** ..................................................................................	
**	 			CONSISTENCY
**	Check consistency for alpha and beta
** ..................................................................................
*/	

	//SET # OF SIMULATIONS 
	iB = 100;			//100
	iSIZE = 10000;		//10000
	iSIMS = iB*iSIZE;
	iSTD_ERROR = 1;
	
	//DO MANY SIMULATIONS AND ESITMATIONS
	decl mOmega, mBeta, mGamma, mXi, mPhi, mTau_1, mTau_2, mSigma2_u, mOmega2, mBeta2, mGamma2, mXi2, mPhi2, mTau_12, mTau_22, mSigma2_u2;
	iSTEPS = iSIZE/10;	//steps of iSIZE/100 takes a while (steps of iSIZE/10 is faster)
	fMonteCarlo2(&mOmega, &mBeta, &mGamma, &mXi, &mPhi, &mTau_1, &mTau_2, &mSigma2_u, &mOmega2, &mBeta2, &mGamma2, &mXi2, &mPhi2, &mTau_12, &mTau_22, &mSigma2_u2);

	//DRAW GRAPHS
	SetDrawWindow("CS_SIM_3_Consistency");
	Draw(0, (mOmega)',iSTEPS,iSTEPS);
	Draw(1, (mBeta)',iSTEPS,iSTEPS);
	Draw(2, (mGamma)',iSTEPS,iSTEPS);
	Draw(3, (mXi)',iSTEPS,iSTEPS);
	Draw(4, (mPhi)',iSTEPS,iSTEPS);
	Draw(5, (mTau_1)',iSTEPS,iSTEPS);
	Draw(6, (mTau_2)',iSTEPS,iSTEPS);
	Draw(7, (mSigma2_u)',iSTEPS,iSTEPS);
	DrawTitle(0,"(i) omega");	
	DrawTitle(1,"(ii) beta");
	DrawTitle(2,"(ii) gamma");	
	DrawTitle(3,"(iv) xi");
	DrawTitle(4,"(v) phi");	
	DrawTitle(5,"(vi) tau_1");
	DrawTitle(6,"(vii) tau_2");	
	DrawTitle(7,"(viii) sigma2_u");
	ShowDrawWindow();
	print("\nSecond Graph Finished at ",time(),"\n");

	SetDrawWindow("CS_SIM_3_Check_Normality_Consistency");
	Draw(0, mOmega2',iSTEPS,iSTEPS);
	Draw(1, mBeta2',iSTEPS,iSTEPS);
	Draw(2, mGamma2',iSTEPS,iSTEPS);
	Draw(3, mXi2',iSTEPS,iSTEPS);
	Draw(4, mPhi2',iSTEPS,iSTEPS);
	Draw(5, mTau_12',iSTEPS,iSTEPS);
	Draw(6, mTau_22',iSTEPS,iSTEPS);
	Draw(7, mSigma2_u2',iSTEPS,iSTEPS);
	DrawTitle(0,"(i) omega");	
	DrawTitle(1,"(ii) beta");
	DrawTitle(2,"(ii) gamma");	
	DrawTitle(3,"(iv) xi");
	DrawTitle(4,"(v) phi");	
	DrawTitle(5,"(vi) tau_1");
	DrawTitle(6,"(vii) tau_2");	
	DrawTitle(7,"(viii) sigma2_u");
	ShowDrawWindow();
	print("\nThird Graph Finished at ",time(),"\n");
}