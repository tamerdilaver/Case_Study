/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Estimate all t-GARCH model parameters 
**		with Maximum Likelikhood for SBUX return s.t. alpha + beta <1 (Since alpha>0 and beta>0) 
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
static decl dALPHA;
static decl dBETA;
static decl dOMEGA;
static decl dGAMMA;
static decl vSTD_NORM;				// Zt ~ N(0,1)
static decl s_vY; 					//Simulated returns
static decl s_vDate;
static decl dALPHA_START;
static decl dBETA_START;
static decl dOMEGA_START;
static decl dLAMBDA_START;
static decl dRATIO;

/*
**  Function:	Transform (start)parameters
**
**  Input: 		vTheta [parametervalues]
**
**  Output: 	vThetaStar
*/

fTransform(const avThetaStar, const vTheta){
	avThetaStar[0]=		vTheta;
	
	avThetaStar[0][0] = log(vTheta[0]/(1-vTheta[0]-vTheta[1]));
	avThetaStar[0][1] = log(vTheta[1]/(1-vTheta[0]-vTheta[1]));
	avThetaStar[0][2] = log(vTheta[2]);
	avThetaStar[0][3] = log(vTheta[3]-4)-log(100-vTheta[3]);
	return 1;
}

/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adAlpha, adBeta, aOmega, adGamma, adLambda, vTheta
**
**  Output: 	1 
*/

fGetPars(const adAlpha, const adBeta, const adOmega,const adLambda, const vTheta){

	adAlpha[0] = exp(vTheta[0])/(exp(vTheta[0])+exp(vTheta[1])+1);
	adBeta[0] = exp(vTheta[1])/(exp(vTheta[0])+exp(vTheta[1])+1);
	adOmega[0] = exp(vTheta[2]);
	adLambda[0] = 4+(100-4)*exp(vTheta[3])/(1+exp(vTheta[3]));
	return 1;
}


/*
**  Function:	Calculates average value loglikelihood for t-GARCH given parameter values
**
**  Input: 		vTheta [parametervalues], adFunc [adres functievalue], avScore [the score], amHessian [hessianmatrix]
**
**  Output:		1
**
*/
fLogLike_t_GARCH(const vTheta, const adFunc, const avScore, const amHessian){
	decl dAlpha, dBeta, dOmega, dLambda;
	fGetPars( &dAlpha,  &dBeta, &dOmega, &dLambda, vTheta);

	decl dS2 = dOmega/(1-dAlpha-dBeta);					 											//initial condition by definition
	decl vLogEta = zeros(sizerc(s_vY), 1);

	//could be made more efficient by deleting the logs
	for(decl i = 0; i < sizerc(s_vY); ++i){
			//likelihood contribution
			vLogEta[i] = -1/2*log(M_PI)-1/2*log(dLambda-2) -1/2*log(dS2) -log(gammafact(dLambda/2))+log(gammafact((dLambda+1)/2)) -(dLambda+1)/2*log(1+ s_vY[i]^2 / ((dLambda-2)*dS2));
			//GARCH recursion
			dS2 = dOmega + dBeta * dS2 +  dAlpha * s_vY[i]^2;
	}
	
	adFunc[0] = sumc(vLogEta)/sizerc(s_vY); 									 	//Average
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

	avTheta[0][0] = exp(vThetaStar[0])/(exp(vThetaStar[0])+exp(vThetaStar[1])+1);
	avTheta[0][1] = exp(vThetaStar[1])/(exp(vThetaStar[0])+exp(vThetaStar[1])+1);
	avTheta[0][2] = exp(vThetaStar[2]);
	avTheta[0][3] = 4+(100-4)*exp(vThetaStar[3])/(1+exp(vThetaStar[3]));
	//actually need to restrict dLambda_hat between (4,100), because we want existence of kurtosis of the model
	//otherwise there will be no convergence for small samples that occur Gaussian
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
		Num2Derivative(fLogLike_t_GARCH, vThetaStar, &mHessian);
		NumJacobian(fTransformBack, vThetaStar, &mJacobian); //numerical Jacobian
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
	decl dAlpha,  dBeta, dOmega, dLambda, vH;

	fGetPars(&dAlpha,  &dBeta, &dOmega, &dLambda, vTheta);
	
	vH = zeros(sizerc(s_vY),1);
	vH[0]= dOmega/(1-dAlpha-dBeta);	
	
	for(decl i = 1; i < sizerc(s_vY); i++){ //mixed 	
		vH[i] = dOmega + dAlpha*s_vY[i-1]^2 + dBeta*vH[i-1];
	}
	
	return 	vH;


}

/*
**  Function:	Estimate t-GARCH parameters
**
**  Input: 		vReturns, adAlpha_hat, adBeta_hat, adOmega_hat, adLambda_hat, avVariance 
**
**  Output: 	vTheta [estimated parametervalues]
*/

fEstimate_t_GARCH(const vReturns, const adAlpha_hat, const adBeta_hat, const adOmega_hat, const adLambda_hat, const avVariance){

	//initialise parameter values
	decl vTheta = zeros(4,1);
	vTheta[0] = dALPHA_START;
	vTheta[1] = dBETA_START;
	vTheta[2] = dOMEGA_START;
	vTheta[3] = dLAMBDA_START;
	
	decl vThetaStart = vTheta;

	//globalize returns and vectorize true pars
	s_vY = vReturns;

	//transform parameters
	decl vThetaStar; 
	fTransform(&vThetaStar, vTheta);

	//Maximize the LL
	decl dFunc;
	decl iA;
	iA=MaxBFGS(fLogLike_t_GARCH, &vThetaStar, &dFunc, 0, TRUE);

	//Transform thetasStar back
  	fTransformBack(&vTheta, vThetaStar);

	//return alpha, beta, omega, gamma and lambda
	adAlpha_hat[0] = vTheta[0];
	adBeta_hat[0] = vTheta[1];
	adOmega_hat[0] = vTheta[2];
	adLambda_hat[0] = vTheta[3];

	decl vSigmaStdError = fSigmaStdError(vThetaStar);
	decl vVariance = fVariance(vThetaStar);
	avVariance[0] = vVariance;
	
	print("\n",MaxConvergenceMsg(iA));
	println("\nFunctiewaarde likelihood eindwaardes:", dFunc);
	print("\nOptimale parameters met standaarderrors \n",
          	"%r", { "alpha",  "beta",  "omega", "lambda"},
          	"%c", {"thetaStart","theta","std.error"}, vThetaStart~vTheta~vSigmaStdError);
			
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
	decl dAlpha,  dBeta, dOmega, dLambda, vH;
	fGetPars(&dAlpha,  &dBeta, &dOmega, &dLambda, vTheta);
	
	vH = zeros((sizerc(s_vY)+1),1);
	vH[0]= dOmega/(1-dAlpha-dBeta);	
	
	for(decl i = 1; i < sizerc(s_vY)+1; i++){			
		vH[i] = dOmega + dAlpha*s_vY[i-1]^2 + dBeta*vH[i-1];
	}
	
	return 	vH[sizerc(s_vY)];
}

/*
**  Function:	Compute MAE
**
**  Input: 		adMAE_OC [adress of MAE], vReturns_1 [return series], vBenchmark [Benchmark], dC [ratio]
**
**  Output: 	1
*/

fMAE(const adMAE_OC, const vReturns_1, const vBenchmark, const dC){

	decl iWindow = 250;
	decl iT = sizerc(vReturns_1);
//	decl vTemp_returns = vReturns_1;
	decl vH_forecast = zeros(iWindow, 1);
	decl vSqrd_error = zeros(iWindow, 1);

	dALPHA_START = 0.1;
	dBETA_START = 0.85;
	dOMEGA_START = 0.01;
	dLAMBDA_START = 10;

	for(decl j = 0; j<iWindow; j++){
		s_vY = 	vReturns_1[j:(iT - iWindow +j)];

		//initialise parametervalues
		decl vTheta = zeros(4,1);
		vTheta[0] = dALPHA_START;
		vTheta[1] = dBETA_START;
		vTheta[2] = dOMEGA_START;
		vTheta[3] = dLAMBDA_START;
	
		//transform parameters
		decl vThetaStar; 
		fTransform(&vThetaStar, vTheta);
	
		//Maximize the LL
		decl dFunc;
		decl iA;
		iA=MaxBFGS(fLogLike_t_GARCH, &vThetaStar, &dFunc, 0, TRUE);
	
		//Transform thetasStar back
	  	fTransformBack(&vTheta, vThetaStar);

		dALPHA_START = vTheta[0];
		dBETA_START = vTheta[1];
		dOMEGA_START = vTheta[2];
		dLAMBDA_START = vTheta[3];

		vH_forecast[j] = fForecast(vThetaStar);
		vSqrd_error[j] = fabs(dC*vH_forecast[j] - dRATIO*vBenchmark[(iT - iWindow +j)]);

	}

	savemat("vAE_CC_RK_GARCH_t.xls", vSqrd_error);
	adMAE_OC[0] = meanc(vSqrd_error);

	return 1;

}

/*
**  Function:	Compute MSE
**
**  Input: 		adMSE_OC [adress of MAE], vReturns_1 [return series], vBenchmark [Benchmark], dC [ratio]
**
**  Output: 	1
*/

fMSE(const adMSE_OC, const vReturns_1, const vBenchmark, const dC){

	decl iWindow = 250;
	decl iT = sizerc(vReturns_1);
//	decl vTemp_returns = vReturns_1;
	decl vH_forecast = zeros(iWindow, 1);
	decl vSqrd_error = zeros(iWindow, 1);

	dALPHA_START = 0.1;
	dBETA_START = 0.85;
	dOMEGA_START = 0.01;

	for(decl j = 0; j<iWindow; j++){
		s_vY = 	vReturns_1[j:(iT - iWindow +j)];

		//initialise parametervalues
		decl vTheta = zeros(4,1);
		vTheta[0] = dALPHA_START;
		vTheta[1] = dBETA_START;
		vTheta[2] = dOMEGA_START;
		vTheta[3] = dLAMBDA_START;
	
		//transform parameters
		decl vThetaStar; 
		fTransform(&vThetaStar, vTheta);
	
		//Maximize the LL
		decl dFunc;
		decl iA;
		iA=MaxBFGS(fLogLike_t_GARCH, &vThetaStar, &dFunc, 0, TRUE);
	
		//Transform thetasStar back
	  	fTransformBack(&vTheta, vThetaStar);

		dALPHA_START = vTheta[0];
		dBETA_START = vTheta[1];
		dOMEGA_START = vTheta[2];

		vH_forecast[j] = fForecast(vThetaStar);
		vSqrd_error[j] = (dC*vH_forecast[j] - dRATIO*vBenchmark[(iT - iWindow +j)])^2;

	}

	adMSE_OC[0] = meanc(vSqrd_error);

	return 1;

}

/*
**				MAIN PROGRAM
**
**  Purpose:	Estimate t-GARCH parameters alpha, beta, omega and lambda.
**
**  Input: 		dALPHA, dBETA, dOMEGA, dLAMBDA, iB, iSIZE, iSIMS, iSTEPS
**
**  Output: 	Figures
*/
main(){
	//laad SBUX returns
	decl mData_1 = loadmat("ReturnsOpenToClose.csv");
	decl mData_2 = loadmat("ReturnsCloseToClose.csv"); 
	decl vReturns_1 = 100*mData_1[:][1];
	decl vReturns_2 = 100*mData_2[:][1];

	dRATIO = (varc(vReturns_1) +varc(vReturns_2))/varc(vReturns_1);
	
	decl vRV = loadmat("RV.csv");
	decl vBV = loadmat("BV.csv");
	decl vRK = loadmat("RK.csv");  

	//laad Dates SBUX returns
	decl vTemp_Date = mData_2[][0];
	decl vYear 		= floor(vTemp_Date/10000);							
	decl vMonth 	= floor((vTemp_Date-floor(vTemp_Date/10000)*10000)/100);	
	decl vDay 		= vTemp_Date-floor(vTemp_Date/100)*100;
	s_vDate 		= dayofcalendar(vYear, vMonth, vDay);

	dALPHA_START = 0.1;
	dBETA_START = 0.85;
	dOMEGA_START = 0.01;
	dLAMBDA_START = 10;
	
	decl dAlpha_hat, dBeta_hat, dOmega_hat, dLambda_hat;
	decl vVariance_1, vVariance_2;
	print("\nO-C");
	fEstimate_t_GARCH(vReturns_1, &dAlpha_hat, &dBeta_hat, &dOmega_hat, &dLambda_hat, &vVariance_1);

	print("\nC-C");
	fEstimate_t_GARCH(vReturns_2, &dAlpha_hat, &dBeta_hat, &dOmega_hat, &dLambda_hat, &vVariance_2);

	//graphs
	SetDrawWindow("CS_EMP_1_t-GARCH(1,1)");
	DrawTMatrix(0, (vReturns_1~sqrt(vVariance_1))', {"Open-to-close"}, s_vDate');
	DrawTMatrix(1, (vReturns_2~sqrt(vVariance_2))', {"Close-to-close"}, s_vDate');
	ShowDrawWindow();

	//forecasts MAE
	decl vBenchmark = vRK;
	decl dMAE_OC;
	fMAE(&dMAE_OC, vReturns_1, vBenchmark, dRATIO);
	print("\n dMAE_OC = ",dMAE_OC);
	
	decl dMAE_CC;
	fMAE(&dMAE_CC, vReturns_2, vBenchmark, 1);
	print("\n dMAE_CC = ",dMAE_CC);

	decl dMSE_OC;
	fMSE(&dMSE_OC, vReturns_1, vBenchmark, dRATIO);
	print("\n dMSE_OC = ",dMSE_OC);
	
	decl dMSE_CC;
	fMSE(&dMSE_CC, vReturns_2, vBenchmark, 1);
	print("\n dMSE_CC = ",dMSE_CC);
	
}