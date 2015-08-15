/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Estimate all Student's t-GARCH model parameters (gamma, nu, omega, A and B)
**
**  Date:
**    	16/01/2015
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
static decl dLAMBDA;				//Degrees of freedom
static decl dALPHA;					
static decl dBETA;
static decl dOMEGA;										 
static decl dGAMMA;					//Initialised Variance
static decl vSTD_STUDENT_T;			// Zt ~ TID(lambda)   
static decl s_vY; 					//Simulated returns
static decl	iSTD_ERROR;				//0 or 1 binary

/*
**  Function:	Simulate t-GARCH returns for given parameters
**
**  Input:		dAlpha, dBeta, dOmega, dGamma, dLambda, avReturns, iIteration [to get different Zt's]
**
**  Output:		1
**
*/

fSim_t_GARCH(const dAlpha, const dBeta, const dOmega, const dGamma, const dLambda, const avReturns, const iIteration){
	decl vTemp, vH;
	vTemp = vH =  zeros(iSIZE+1, 1);

	vH[0]= dGamma;		//by definition
	
	for(decl i = 0; i < iSIZE; i++){	
		vTemp[i] =  sqrt(vH[i])*vSTD_STUDENT_T[(i + (iIteration*iSIZE))];
		vH[i+1] = dOmega + dBeta*vH[i] + dAlpha*sqr(vTemp[i]) ;
	}

	vTemp = dropr(vTemp,iSIZE);
	vH = dropr(vH,iSIZE);
	
	avReturns[0] = vTemp;
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
	
	avThetaStar[0][0] = log(vTheta[0]/(1-vTheta[0]-vTheta[1]));
	avThetaStar[0][1] = log(vTheta[1]/(1-vTheta[0]-vTheta[1]));
	avThetaStar[0][2] = log(vTheta[2]);
	avThetaStar[0][3] = log(vTheta[3]);
	avThetaStar[0][4] = log(vTheta[4]-2)-log(100-vTheta[4]);
	return 1;
}

/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adAlpha, adBeta, aOmega, adGamma, adLambda, vTheta
**
**  Output: 	1 
*/

fGetPars(const adAlpha, const adBeta, const adOmega, const adGamma,const adLambda, const vTheta){

	adAlpha[0] = exp(vTheta[0])/(exp(vTheta[0])+exp(vTheta[1])+1);
	adBeta[0] = exp(vTheta[1])/(exp(vTheta[0])+exp(vTheta[1])+1);
	adOmega[0] = exp(vTheta[2]);
	adGamma[0] = exp(vTheta[3]);
	adLambda[0] = 2+(100-2)*exp(vTheta[4])/(1+exp(vTheta[4]));
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
	decl dAlpha, dBeta, dOmega, dGamma, dLambda;
	fGetPars( &dAlpha,  &dBeta, &dOmega,  &dGamma, &dLambda, vTheta);

	decl dS2 = dGamma; //initial condition by definition
	decl vLogEta = zeros(sizerc(s_vY), 1);

	//Could be made more efficient by deleting logs
	for(decl i = 0; i < sizerc(s_vY); ++i){
			//likelihood contribution
			vLogEta[i] = -1/2*log(M_PI)-1/2*log(dLambda-2) -1/2*log(dS2) -log(gammafact(dLambda/2))+log(gammafact((dLambda+1)/2)) -(dLambda+1)/2*log(1+ s_vY[i]^2 / ((dLambda-2)*dS2));
			//GARCH recursion
			dS2 = dOmega + dBeta * dS2 +  dAlpha * s_vY[i]^2;
	}
	
	adFunc[0] = sumc(vLogEta)/sizerc(s_vY); //Average
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
	avTheta[0] =	vThetaStar;

	avTheta[0][0] = exp(vThetaStar[0])/(exp(vThetaStar[0])+exp(vThetaStar[1])+1);
	avTheta[0][1] = exp(vThetaStar[1])/(exp(vThetaStar[0])+exp(vThetaStar[1])+1);
	avTheta[0][2] = exp(vThetaStar[2]);
	avTheta[0][3] = exp(vThetaStar[3]);
	avTheta[0][4] = 2+(100-2)*exp(vThetaStar[4])/(1+exp(vThetaStar[4]));
	//actually need to restrict dLambda_hat between (2,100)
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

		iN = sizerc(s_vY); //number of observations
		Num2Derivative(fLogLike_t_GARCH, vThetaStar, &mHessian);
		NumJacobian(fTransformBack, vThetaStar, &mJacobian);	  //numerical Jacobian
		mHessian 	= mJacobian*invert(-iN*mHessian)*mJacobian';
		vStdErrors 	= sqrt(diagonal(mHessian)');

		return 	vStdErrors;
}

/*
**  Function:	Estimate t-GARCH parameters
**
**  Input: 		vReturns, adAlpha_hat, adBeta_hat, adOmega_hat, adGamma_hat, adLambda_hat (dBeta_0 not necessary)
**
**  Output: 	vTheta [estimated parametervalues]
*/

fEstimate_t_GARCH(const vReturns, const adAlpha_hat, const adBeta_hat, const adOmega_hat, const adGamma_hat, const adLambda_hat){

	//initialise parametervalues
	decl vTheta = zeros(5,1);
	vTheta = <0.1 ; 0.89 ; 0.05 ; 5; 7>; // Alpha, Beta, Omega, Gamma, Lambda starting values
	decl vThetaStart = vTheta;

	//globalalize returns and vectorize true pars
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
	adGamma_hat[0] = vTheta[3];
	adLambda_hat[0] = vTheta[4];
	
	if(iSTD_ERROR){		//only do this for fMonteCarlo2
		decl vSigmaStdError = fSigmaStdError(vThetaStar);
		return vSigmaStdError;
	}else{
		return 1;  //otherwise end function by returning 1
	}
}

/*
**  Function:	Simulates and Estimates t-GARCH data and parameters many times
**				to illustrate Asymptotic normality
**
**  Input: 		amMonteCarlo [matrix of many estimated parameters],  dBeta_0;
**
**  Output: 	1
*/

fMonteCarlo(const amMonteCarlo){
	decl mTemp;
	mTemp = zeros(iB,5);

	for(decl i = 0; i<iB ; i++){
		decl vReturns;
		fSim_t_GARCH(dALPHA, dBETA, dOMEGA, dGAMMA, dLAMBDA, &vReturns, i);

		decl dAlpha_hat, dBeta_hat, dOmega_hat, dGamma_hat, dLambda_hat;
		fEstimate_t_GARCH(vReturns, &dAlpha_hat, &dBeta_hat, &dOmega_hat, &dGamma_hat, &dLambda_hat);	 //Omega, Gamma and lambda also estimated

		mTemp[i][0] =  dAlpha_hat;
		mTemp[i][1]	=  dBeta_hat;
		mTemp[i][2] =  dOmega_hat;
		mTemp[i][3]	=  dGamma_hat;
		mTemp[i][4]	=  dLambda_hat;
	}
	amMonteCarlo[0] = mTemp;
	return 1;
}

/*
**  Function:	Simulated and Estimates t-GARCH data and parameters many times
**				to illustrate consistency it returns minimum, mean and maximum values for the estimated parameters
**
**  Input: 		amAlpha [matrix containing the min, max and mean of estimated alpha],
**				amBeta [matrix containing the min, max and mean of estimated beta], 
**				amOmega [matrix containing the min, max and mean of estimated omega],
**				amGamma [matrix containing the min, max and mean of estimated gamma], 
**				amLambda [matrix containing the min, max and mean of estimated lambda],	dBETA
**
**  Output: 	1
*/

fMonteCarlo2(const amAlpha, const amBeta, const amOmega, const amGamma, const amLambda, const amAlpha2, const amBeta2, const amOmega2, const amGamma2, const amLambda2){
	decl mTemp, mTempAlpha, mTempBeta, mTempOmega, mTempGamma, mTempLambda;
	decl mTemp2, mTempAlpha2, mTempBeta2, mTempOmega2, mTempGamma2, mTempLambda2;
	mTempAlpha = mTempBeta = mTempOmega = mTempGamma = mTempLambda = zeros((iSIZE/iSTEPS),3);
	mTempAlpha2 = mTempBeta2 = mTempOmega2 = mTempGamma2 = mTempLambda2 = zeros((iSIZE/iSTEPS),3);
	mTemp = mTemp2 = zeros(iB,5);

	decl iSize = iSIZE;

	for(decl j = 0; j<(iSize/iSTEPS) ; j++){
		iSIZE = ((iSTEPS)*(j+1));
		for(decl i = 0; i<iB ; i++){
			decl vReturns;
			fSim_t_GARCH(dALPHA, dBETA, dOMEGA, dGAMMA, dLAMBDA, &vReturns, i);
	
			decl dAlpha_hat, dBeta_hat, dOmega_hat, dGamma_hat, dLambda_hat, vSE;
			vSE = fEstimate_t_GARCH(vReturns, &dAlpha_hat, &dBeta_hat, &dOmega_hat, &dGamma_hat, &dLambda_hat);	 //Omega, Gamma and Lambda also estimated
	
			mTemp[i][0] =  sqrt(iSIZE)*(dAlpha_hat - dALPHA); //SQRT(T)*(\hat_\alpha_T - \alpha_0) ~ N(0, \SIGMA)
			mTemp[i][1]	=  sqrt(iSIZE)*(dBeta_hat - dBETA);
			mTemp[i][2]	=  sqrt(iSIZE)*(dOmega_hat - dOMEGA);
			mTemp[i][3]	=  sqrt(iSIZE)*(dGamma_hat - dGAMMA);
			mTemp[i][4]	=  sqrt(iSIZE)*(dLambda_hat - dLAMBDA);

			//2nd part
			mTemp2[i][0] 	=  (dAlpha_hat-dALPHA)/vSE[0]; //(\hat_\alpha_T - \alpha_0)/SE(\hat_\alpha) ~ N(0, 1)
			mTemp2[i][1]	=  (dBeta_hat-dBETA)/vSE[1];
			mTemp2[i][2]	=  (dOmega_hat-dOMEGA)/vSE[2];
			mTemp2[i][3]	=  (dGamma_hat-dGAMMA)/vSE[3];
			mTemp2[i][4]	=  (dLambda_hat - dLAMBDA)/vSE[4];
			
		}
		// v0.025_quantile, vMean, v0.975_quantile --> we get 95%-intervals
		mTempAlpha[j][0] = quantilec(mTemp[][],0.025)'[0];
		mTempAlpha[j][1] = meanc(mTemp[][])'[0];
		mTempAlpha[j][2] = quantilec(mTemp[][],0.975)'[0];
	
		mTempBeta[j][0] = quantilec(mTemp[][],0.025)'[1];
		mTempBeta[j][1] = meanc(mTemp[][])'[1];
		mTempBeta[j][2] = quantilec(mTemp[][],0.975)'[1];

		mTempOmega[j][0] = quantilec(mTemp[][],0.025)'[2];
		mTempOmega[j][1] = meanc(mTemp[][])'[2];
		mTempOmega[j][2] = quantilec(mTemp[][],0.975)'[2];
	
		mTempGamma[j][0] = quantilec(mTemp[][],0.025)'[3];
		mTempGamma[j][1] = meanc(mTemp[][])'[3];
		mTempGamma[j][2] = quantilec(mTemp[][],0.975)'[3];

		mTempLambda[j][0] = quantilec(mTemp[][],0.025)'[4];
		mTempLambda[j][1] = meanc(mTemp[][])'[4];
		mTempLambda[j][2] = quantilec(mTemp[][],0.975)'[4];

		//2nd part: v0.025_quantile, v0.5_quantile, v0.975_quantile;	
		mTempAlpha2[j][0] = quantilec(mTemp2[][],0.025)'[0];
		mTempAlpha2[j][1] = quantilec(mTemp2[][],0.5)'[0];	  //deletec()
		mTempAlpha2[j][2] = quantilec(mTemp2[][],0.975)'[0];
	
		mTempBeta2[j][0] = quantilec(mTemp2[][],0.025)'[1];
		mTempBeta2[j][1] = quantilec(mTemp2[][],0.5)'[1];
		mTempBeta2[j][2] = quantilec(mTemp2[][],0.975)'[1];

		mTempOmega2[j][0] = quantilec(mTemp2[][],0.025)'[2];
		mTempOmega2[j][1] = quantilec(mTemp2[][],0.5)'[2];
		mTempOmega2[j][2] = quantilec(mTemp2[][],0.975)'[2];
	
		mTempGamma2[j][0] = quantilec(mTemp2[][],0.025)'[3];
		mTempGamma2[j][1] = quantilec(mTemp2[][],0.5)'[3];
		mTempGamma2[j][2] = quantilec(mTemp2[][],0.975)'[3];

		mTempLambda2[j][0] = quantilec(mTemp2[][],0.025)'[4];
		mTempLambda2[j][1] = quantilec(mTemp2[][],0.5)'[4];
		mTempLambda2[j][2] = quantilec(mTemp2[][],0.975)'[4];
	}

	amAlpha[0] 	= mTempAlpha;
	amBeta[0] 	= mTempBeta;
	amOmega[0] 	= mTempOmega;
	amGamma[0] 	= mTempGamma;
	amLambda[0]	= mTempLambda;

	amAlpha2[0] = mTempAlpha2;
	amBeta2[0] 	= mTempBeta2;
	amOmega2[0] = mTempOmega2;
	amGamma2[0] = mTempGamma2;
	amLambda2[0]= mTempLambda2;

	return 1;
}

/*
**				MAIN PROGRAM
**
**  Purpose:	For a given alpha, omega, gamma and lambda, get beta such that E log(alpha_0 z_t^2 + beta_0) = 0.
**				Simulate t-GARCH returns for alpha, omega, gamma, lambda and the found beta many time.
**				Estimate t-GARCH parameters alpha, beta, omega, gamma and lambda.
**
**  Input: 		dALPHA, dOMEGA, dGAMMA, dLAMBDA, iB, iSIZE, iSIMS, iSTEPS
**
**  Output: 	Figures
*/
main()
{
	//SET PARAMETERS
	dALPHA = 0.1;
	dBETA = 0.89;
	dOMEGA = 0.05;
	dGAMMA = 5;
	dLAMBDA = 7;


///*
//** ..................................................................................	
//**	 		ASYMPTOTIC NORMALITY
//**	Get distributions of alpha and beta (to check for asymptotic normality)
//**..................................................................................
//*/

	//SET # OF SIMULATIONS 
	iB = 500;			//5000 possible 
	iSIZE = 5000;		//5000 possible
	iSIMS = iB*iSIZE;
	vSTD_STUDENT_T = sqrt((dLAMBDA-2)/dLAMBDA)*rant(iSIMS,1, dLAMBDA);	 //draw standardized student's t distributed 
	iSTD_ERROR = 0;
	
	//DO MANY SIMULATIONS AND ESITMATIONS	
	decl mMonteCarlo;
	fMonteCarlo(&mMonteCarlo);	  

	//DRAW GRAPHS
	SetDrawWindow("CS_SIM_1_asymptotic_normality");
	DrawDensity(0, (mMonteCarlo[][0])', {"(i) Density alpha"});
	DrawDensity(1, (mMonteCarlo[][1])', {"(ii) Density beta"});
	DrawDensity(2, (mMonteCarlo[][2])', {"(iii) Density omega"});
	DrawDensity(3, (mMonteCarlo[][3])', {"(iv) Density gamma"});
	DrawDensity(4, (mMonteCarlo[][4])', {"(iv) Density lambda"});
	ShowDrawWindow();

	print("\nFirst Graph Finished at ",time(),"\n");
/*
** ..................................................................................	
**	 			CONSISTENCY
**	Check consistency for alpha and beta
** ..................................................................................
*/	

	//SET # OF SIMULATIONS 
	iB = 100;			  //100
	iSIZE = 10000;		  //10000
	iSIMS = iB*iSIZE;
	vSTD_STUDENT_T = sqrt((dLAMBDA-2)/dLAMBDA)*rant(iSIMS,1,dLAMBDA);	 //draw standardized student's t distributed 
	iSTD_ERROR = 1;
		
	//DO MANY SIMULATIONS AND ESITMATIONS
	decl mAlpha, mBeta, mOmega, mGamma, mLambda, mAlpha2, mBeta2, mOmega2, mGamma2, mLambda2;
	iSTEPS = iSIZE/10;				 	//steps of iSIZE/100 takes a while (steps of iSIZE/10 is faster)
	fMonteCarlo2(&mAlpha, &mBeta, &mOmega, &mGamma, &mLambda, &mAlpha2, &mBeta2, &mOmega2, &mGamma2, &mLambda2);

	//DRAW GRAPHS
	SetDrawWindow("CS_SIM_1_Consistency");
	Draw(0, mAlpha',iSTEPS,iSTEPS);
	Draw(1, mBeta',iSTEPS,iSTEPS);
	Draw(2, mOmega',iSTEPS,iSTEPS);
	Draw(3, mGamma',iSTEPS,iSTEPS);
	Draw(4, mLambda',iSTEPS,iSTEPS);
	DrawTitle(0,"(i) alpha");	
	DrawTitle(1,"(ii) beta");
	DrawTitle(2,"(ii) omega");	
	DrawTitle(3,"(iv) gamma");
	DrawTitle(4,"(v) lambda");
	ShowDrawWindow();
	print("\nSecond Graph Finished at ",time(),"\n");

	SetDrawWindow("CS_SIM_1_Check_Normality_Consistency");
	Draw(0, mAlpha2',iSTEPS,iSTEPS);
	Draw(1, mBeta2',iSTEPS,iSTEPS);
	Draw(2, mOmega2',iSTEPS,iSTEPS);
	Draw(3, mGamma2',iSTEPS,iSTEPS);
	Draw(4, mLambda2',iSTEPS,iSTEPS);
	DrawTitle(0,"(i) alpha");	
	DrawTitle(1,"(ii) beta");
	DrawTitle(2,"(ii) omega");	
	DrawTitle(3,"(iv) gamma");
	DrawTitle(4,"(iv) lambda");
	ShowDrawWindow();
	print("\nThird Graph Finished at ",time(),"\n");
}