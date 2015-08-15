#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import <maximize>

/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Estimate all Log linear Real t-GARCH model parameters using Bayesian estimation
**
**  Date:
**    	28/01/2015
**
**  Author:
**	  	Tamer Dilaver, Koen de Man & Sina Zolnoor
**
**	Supervisor:
**		L.F. Hoogerheide & S.J. Koopman
**
*/


static decl s_dMaxLogLikelihood, s_vR, s_vRV, s_dNumberOfParameters, s_mEstimatedCovMatrixMLE, s_vLikelihoodKernel, s_vLogL, s_PLL;
 
/*
**  Function: 	Extract the parameters from vTheta
**
**  Input: 		adOmega, adBeta, adGamma, adXi, adPhi, adTau_1, adTau_2, adSigma2_u,  vTheta
**				0		1	     2	      3	     4	    5	    6	     7		  
*/

fGetPars(const adOmega, const adBeta, const adGamma, const adXi, const adPhi, const adTau_1, const adTau_2, const adSigma2_u, const adLambda_z, const adLambda_u, const vTheta){

	adOmega[0] = vTheta[0];
	adBeta[0] = exp(vTheta[2])/(exp(vTheta[1])+exp(vTheta[2])+1);
	adGamma[0] = exp(vTheta[1]-vTheta[4])/(exp(vTheta[1])+exp(vTheta[2])+1);
	adXi[0]	= vTheta[3];
	adPhi[0]= exp(vTheta[4]);
	adTau_1[0] 	= vTheta[5];
	adTau_2[0] 	= vTheta[6];
	adSigma2_u[0] = exp(vTheta[7]);
	adLambda_z[0] = 4+(100-4)*exp(vTheta[8])/(1+exp(vTheta[8]));
	adLambda_u[0] = 4+(100-4)*exp(vTheta[9])/(1+exp(vTheta[9]));
}

/*
**
**	Purpose:
**	  Transform parameters to estimate with the restrictions:
**
**
**	Inputs:
**	  vTheta: 10 x 1 vector which contains parameters 
**
**
*/

fTransform(const avThetaStar, const vTheta){
	avThetaStar[0]=		vTheta;
	avThetaStar[0][0] = vTheta[0];
	avThetaStar[0][1] = log((vTheta[2]*vTheta[4])/(1-vTheta[2]*vTheta[4] - vTheta[1]));	
	avThetaStar[0][2] = log((vTheta[1])/(1-vTheta[2]*vTheta[4] - vTheta[1]));			
	avThetaStar[0][3] = vTheta[3];
	avThetaStar[0][4] = log(vTheta[4]);
	avThetaStar[0][5] = vTheta[5];
	avThetaStar[0][6] = vTheta[6];
	avThetaStar[0][7] = log(vTheta[7]);		  
	avThetaStar[0][8] = log(vTheta[8]-4)-log(100-vTheta[8]);
	avThetaStar[0][9] = log(vTheta[9]-4)-log(100-vTheta[9]);

	return 1;
}

/*
**
**	Purpose:
**	  Transform parameters back
**
**
**	Inputs:
**	  vThetaStar: 10 x 1 vector which contains parameters 
**
**
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
	avTheta[0][8] = 4+(100-4)*exp(vThetaStar[8])/(1+exp(vThetaStar[8]));
	avTheta[0][9] = 4+(100-4)*exp(vThetaStar[9])/(1+exp(vThetaStar[9]));

	return 1;
}

/*
**
**	Purpose:
**    Computes a vector of variances needed for the likelihoodfunction
**
**
**	Inputs:
**    vP: 10 x 1 vector of parameters
**
**	Returns: exp(vLogH): vector of asked variances
**
*/

fVariance(const vTheta){
	decl dOmega, dBeta, dGamma, dXi, dPhi, dTau_1, dTau_2, dSigma2_u, dLambda_z, dLambda_u;
	fGetPars(&dOmega, &dBeta, &dGamma, &dXi, &dPhi, &dTau_1, &dTau_2, &dSigma2_u, &dLambda_z, &dLambda_u, vTheta);

	decl iNumberOfObserations, vLogH;
	iNumberOfObserations = sizerc(s_vR);
	vLogH = zeros(iNumberOfObserations,1);
	vLogH[0] = (dOmega +dGamma*dXi)/(1-dBeta-dPhi*dGamma);	

	decl i;
	for(i=1; i<iNumberOfObserations; i++){													   
		vLogH[i] = dOmega  + dBeta*vLogH[i-1]+ dGamma*log(s_vRV[i-1]);
	}		   
	
	return exp(vLogH);
}

/*
**
**	Purpose: 		
**	  Computes loglikelihood for given parameters
**
**
**  Inputs:
**    vPtrans			10 x 1 vector of transformed parameters
**	  						
**  Return Value:
**    iBoolean			integer, boolean value is 1 if loglikelihood is succesfully computed, 0 if not
**	  	
*/

fLogLike_LogReal_t_GARCH(const vTheta, const adFunc, const avScore, const amHessian){
	decl dOmega, dBeta, dGamma, dXi, dPhi, dTau_1, dTau_2, dSigma2_u, dLambda_z, dLambda_u;
	fGetPars(&dOmega, &dBeta, &dGamma, &dXi, &dPhi, &dTau_1, &dTau_2, &dSigma2_u, &dLambda_z, &dLambda_u, vTheta);

	decl dlogH, vlogEta;
	dlogH = (dOmega +dGamma*dXi)/(1-dBeta-dPhi*dGamma); 
	vlogEta = zeros(sizerc(s_vR), 1);

	decl i;
	for(i = 0; i < sizerc(s_vR); ++i){
			vlogEta[i] = -log(M_PI)-1/2*log(dLambda_z-2) -1/2*dlogH -log(gammafact(dLambda_z/2))+log(gammafact((dLambda_z+1)/2)) -(dLambda_z+1)/2*log(1+ s_vR[i]^2 / ((dLambda_z-2)*exp(dlogH))) -1/2*log(dLambda_u-2) -1/2*log(dSigma2_u) -log(gammafact(dLambda_u/2))+log(gammafact((dLambda_u+1)/2)) -(dLambda_u+1)/2*log(1+ (log(s_vRV[i]) - dXi - dPhi*dlogH - dTau_1*s_vR[i]/sqrt(exp(dlogH)) - dTau_2*(s_vR[i]^2/exp(dlogH)-1))^2 / ((dLambda_u-2)*dSigma2_u));		
			dlogH = dOmega + dBeta * dlogH + dGamma * log(s_vRV[i]);
	}

	decl iN, vLogVar;
	iN = rows(s_vR);
	vLogVar = log(fVariance(vTheta));
	s_PLL = (iN*log(gammafact((dLambda_z+1)/2)) - iN*log(gammafact(dLambda_z/2)) - 1/2*iN*log(M_PI*(dLambda_z-2)) - 1/2*sumc(vLogVar) - ((dLambda_z+1)/2)*sumc(log(1+(s_vR.^2 ./((dLambda_z-2).*exp(vLogVar))))));

	adFunc[0] = sumc(vlogEta)/sizerc(s_vR); 									
	return 1;
}

/*	fMultivariateTDraws(const vMu, const dDF, const iNumberOfDraws)	
**
**	Purpose: 		
**	  Computes draws from the multivariate student t distribution
**
**
**  Inputs:
**    vMu				vector of expected parameter values from the draws
**	  dDF				degrees of freedom
**    iNumberOfDraws	number of draws
**
**  Output:
**    mStudent			matrix of candidate draws
**	  	
*/

fMultivariateTDraws(const vMu, const dDF, const iNumberOfDraws){															   	

	decl mR, mStdNormal, dAuxFactor, mStudent;
  	mR 	= choleski(s_mEstimatedCovMatrixMLE);														   	
	mStdNormal 					= rann(iNumberOfDraws, sizerc(vMu));									
  	dAuxFactor 					= ones(iNumberOfDraws, 1) ./ sqrt( ranchi(iNumberOfDraws, 1, dDF) / dDF );
  	mStudent 					= ones(iNumberOfDraws, 1) * vMu' + (dAuxFactor*ones(1,sizerc(vMu))) .* mStdNormal * mR';		

	return mStudent;
}

/*	fMultivariateT(const vMu, const dDF, const mCandidateDraws)
**
**	Purpose: 		
**	  Computes the candidate probability density function (pdf)
**
**
**  Inputs:
**    vMu				vector of expected parameter values from the draws
**	  dDF				degrees of freedom
**    mCandidateDraws	matrix of candidate draws
**
**  Output:		
**    vPdf				vector of the candidate probability density function (pdf) == Q(theta~)
**	  	
*/

fMultivariateT(const vMu, const dDF, const mCandidateDraws){

	decl mInvCapSigma, dDetermCapSigma, mP; 
	mInvCapSigma = invert(s_mEstimatedCovMatrixMLE);								
	dDetermCapSigma	= determinant(s_mEstimatedCovMatrixMLE);
	mP = choleski(mInvCapSigma);
	
	decl iTotalNumberOfDraws, vPdf, vIota, vX_minus_mu_InvCapSigma_X_minus_mu;
	iTotalNumberOfDraws = rows(mCandidateDraws);								
	vIota = ones(iTotalNumberOfDraws,1);																		
	vX_minus_mu_InvCapSigma_X_minus_mu = sumsqrr((mCandidateDraws-vIota*vMu' ) * mP);
	vPdf = (gammafact((dDF + s_dNumberOfParameters)/2) / (gammafact(dDF/2)*(M_PI*dDF)^(s_dNumberOfParameters/2))*dDetermCapSigma^(-0.5))	
								  * ( vIota+(1/dDF)*vX_minus_mu_InvCapSigma_X_minus_mu ).^(-0.5*(dDF+s_dNumberOfParameters));		 					

	return vPdf;
}

/*	
**
**	Purpose: 		
**	  Computes the one step transition probability to get p(theta | y)
**			
**  Input:
**    mCandidateDrawsvP				matrix of multivariate student t candidate draws
**
**  Output:		
**    vPosteriorKernelValue			p(theta | y)
**	  	
*/

fGarchTPosteriorKernelMatrix(const mCandidateDrawsvP){

	decl iNumberOfObservations, iNumberOfDraws; 
	iNumberOfObservations = sizerc(s_vR);
	iNumberOfDraws = rows(mCandidateDrawsvP);

	decl vOmega, vBeta, vGamma, vXi, vPhi, vTau1, vTau2, vSigmaSquared_u, vLambda_z, vLambda_u;   
	vOmega = mCandidateDrawsvP[][0];										
	vBeta = mCandidateDrawsvP[][1];										
 	vGamma = mCandidateDrawsvP[][2];										
 	vXi = mCandidateDrawsvP[][3];										
	vPhi = mCandidateDrawsvP[][4];										
	vTau1 = mCandidateDrawsvP[][5];										
	vTau2 = mCandidateDrawsvP[][6];										
	vSigmaSquared_u	= mCandidateDrawsvP[][7];										
	vLambda_z = mCandidateDrawsvP[][8];										
	vLambda_u = mCandidateDrawsvP[][9];										

	decl vIndicatorParamaterValues, vEasyValues, mCandidateDrawsvPReplaced;
	vIndicatorParamaterValues = (vBeta .>= 0) .* (vGamma .>= 0) .* (vPhi .>= 0) .* (vSigmaSquared_u .>= 0) .* (vLambda_z .>= 2) .* (vLambda_u .>= 0).* (vBeta + vGamma .* vPhi .< 1);
	vEasyValues = varc(s_vR) ~ 0 ~ 0 ~ 0 ~ 0 ~ 0 ~ 0 ~ 0.25 ~ 13 ~ 7;																										
	mCandidateDrawsvPReplaced = (vIndicatorParamaterValues * ones(1,s_dNumberOfParameters)) .* mCandidateDrawsvP + (vIndicatorParamaterValues .==0) * vEasyValues;					
							
	vOmega = mCandidateDrawsvPReplaced[][0];										
	vBeta = mCandidateDrawsvPReplaced[][1];										
 	vGamma = mCandidateDrawsvPReplaced[][2];										
 	vXi = mCandidateDrawsvPReplaced[][3];										
	vPhi = mCandidateDrawsvPReplaced[][4];										
	vTau1 = mCandidateDrawsvPReplaced[][5];										
	vTau2 = mCandidateDrawsvPReplaced[][6];										
	vSigmaSquared_u	= mCandidateDrawsvPReplaced[][7];										
	vLambda_z = mCandidateDrawsvPReplaced[][8];										
	vLambda_u = mCandidateDrawsvPReplaced[][9];							

	decl vIota, vVar, vCorrectedLambda_z, vCorrectedLambda_u;
	vIota = ones(iNumberOfDraws,1);										
	s_vLogL	= zeros(iNumberOfDraws,1);								
	vVar = varc(s_vR)*vIota; 											
	vCorrectedLambda_z = (vLambda_z - 2*vIota)./vLambda_z;										
	vCorrectedLambda_u = (vLambda_u - 2*vIota)./vLambda_u;
	
	decl t, vRViota, vResidualT, vLogLAdditionT;
	for(t=1; t<iNumberOfObservations; t++){																				   
	   vRViota = s_vRV[t-1] * vIota; 											
	   vResidualT = s_vR[t] * vIota;										
	   vVar	= vOmega + vBeta .*log(vVar) + vGamma .* log(vRViota);																														
	   vLogLAdditionT = -0.5*log(vVar .*vCorrectedLambda_z) + log(denst(vResidualT./ sqrt (vVar .* vCorrectedLambda_z), vLambda_z))
	    				-0.5*log(vSigmaSquared_u .*vCorrectedLambda_u) + log(denst(vResidualT./ sqrt (vSigmaSquared_u .* vCorrectedLambda_u), vLambda_u));
	   s_vLogL = s_vLogL + vLogLAdditionT;
	}
	
	decl vLikelihoodKernel, vPrior, vPosteriorKernelValue; 
	vLikelihoodKernel = exp(s_vLogL - s_PLL);
	vPrior = vIndicatorParamaterValues .* ((0.05*exp(-(vLambda_z-2)/20)).*(0.05*exp(-(vLambda_u-2)/20)));
	println(sumc(vIndicatorParamaterValues));
	vPosteriorKernelValue = vLikelihoodKernel .* vPrior;

	return vPosteriorKernelValue;
}

/*
**  Function:	calculate standard errors
**
**  Input: 		vThetaStar
**
**  Output: 	vStdErrors
*/

fSigmaStdError(const vThetaStar){

 		decl mHessian, mJacobian;
		Num2Derivative(fLogLike_LogReal_t_GARCH, vThetaStar, &mHessian);
		NumJacobian(fTransformBack, vThetaStar, &mJacobian);

		decl iNumberOfObservations, vStdErrors;
		s_mEstimatedCovMatrixMLE = invert(-mHessian);
		iNumberOfObservations 	= sizerc(s_vR);
		mHessian = mJacobian*invertsym(-iNumberOfObservations*mHessian)*mJacobian';
		vStdErrors 	= sqrt(diagonal(mHessian)');

		return 	vStdErrors;
}

main(){

	decl mData, vRV, vBV, vKV; 
	mData = loadmat("ReturnsOpenToClose.csv");
	s_vR = 100*mData[][1];
	vRV	= loadmat("RV.csv");		  				
	vBV	= loadmat("BV.csv");		 					
	vKV = loadmat("RK.csv");	  						
	s_vRV = vRV;												

	decl vRhulp, vRVhulp, iNumberOfObservations, vPStart, vP;
	s_dNumberOfParameters = 10;														
	vRhulp = s_vR;								 							
	vRVhulp = s_vRV;													
	iNumberOfObservations = rows(s_vR);												     			  	      
	vPStart = vP = <0.039938 ; 0.70111 ; 0.29256 ; -0.095942 ; 0.99598 ; -0.0066171 ; 0.099433 ; 0.21021 ; 7 ; 7>;	// Give starting values for the parameters. These are equal to optimal Gaussian estimates
	
	decl vPtrans, dLikelihood, ir;
	fTransform(&vPtrans, vP);															
	fLogLike_LogReal_t_GARCH(vPtrans, &dLikelihood, 0, 0);										
	ir = MaxBFGS(fLogLike_LogReal_t_GARCH, &vPtrans, &dLikelihood, 0, 1);		
	s_dMaxLogLikelihood = dLikelihood;													

	decl vStdErrors; 
	fTransformBack(&vP,vPtrans);															
	vStdErrors	= fSigmaStdError(vPtrans);										

	decl iBurnInDraws, iEffectiveNumberOfDraws, iTotalNumberOfDraws, iCandidateTDF, mCandidateDraws, vCandidatePdf, vPosteriorKernelValue, vImportanceWeights, dWeightLastAcceptedCandidateDraw, mMHDraws, iNumberOfAcceptedCandidateDraws; 
	iBurnInDraws = 1000;														
	iEffectiveNumberOfDraws = 25000; 														
	iTotalNumberOfDraws = iBurnInDraws + iEffectiveNumberOfDraws;  						
	iCandidateTDF = 10;													
	mCandidateDraws = zeros(iTotalNumberOfDraws, s_dNumberOfParameters); 						
	mCandidateDraws = fMultivariateTDraws(vP, iCandidateTDF, iTotalNumberOfDraws);	
	vCandidatePdf = fMultivariateT(vP, iCandidateTDF, mCandidateDraws);
	vPosteriorKernelValue = fGarchTPosteriorKernelMatrix(mCandidateDraws);
	vImportanceWeights = vPosteriorKernelValue ./ vCandidatePdf;
	dWeightLastAcceptedCandidateDraw = max(DBL_MIN, vImportanceWeights[0]);
	mMHDraws = mCandidateDraws;												
	iNumberOfAcceptedCandidateDraws = 0;													

	decl t, dU, dRatioAcceptedDraws, mPosteriorDraws, vPosteriorThetaMean, vPosteriorThetaVar;
	for(t=1; t<iTotalNumberOfDraws; t++){

		dU = ranu(1,1);
		if(dU < (vImportanceWeights[t] / dWeightLastAcceptedCandidateDraw)){			
			mMHDraws[t][] = mCandidateDraws[t][]; 	   								
			dWeightLastAcceptedCandidateDraw = vImportanceWeights[t];					
			iNumberOfAcceptedCandidateDraws++;												
		}else{
			mMHDraws[t][] = mMHDraws[t-1][];												
		}
	}
				
	dRatioAcceptedDraws = iNumberOfAcceptedCandidateDraws / (iTotalNumberOfDraws-1);
	mPosteriorDraws	= mMHDraws[iBurnInDraws:iTotalNumberOfDraws-1][];					
	vPosteriorThetaMean	= meanc(mPosteriorDraws);										
	vPosteriorThetaVar = varc(mPosteriorDraws);										

	println("Bayesian Estimation: ", vPosteriorThetaMean'~sqrt(vPosteriorThetaVar)');
	println("ML Estimation: ", vP~vStdErrors);

	SetDrawWindow("BAYES DENSITIES");
	DrawDensity(0, (mPosteriorDraws[][0])', {"(ii) Density omega"});
	DrawDensity(1, (mPosteriorDraws[][1])', {"(iii) Density beta"});
	DrawDensity(2, (mPosteriorDraws[][2])', {"(iv) Density gamma"});
	DrawDensity(3, (mPosteriorDraws[][3])', {"(v) Density xi"});
	DrawDensity(4, (mPosteriorDraws[][4])', {"(vi) Density phi"});
	DrawDensity(5, (mPosteriorDraws[][5])', {"(vii) Density tau_1"});
	DrawDensity(6, (mPosteriorDraws[][6])', {"(viii) Density tau_2"});
	DrawDensity(7, (mPosteriorDraws[][7])', {"(ix) Density sigma2_u"});
	DrawDensity(8, (mPosteriorDraws[][8])', {"(x) Density lambda_z"});
	DrawDensity(9, (mPosteriorDraws[][9])', {"(xi) Density lambda_u"});
	ShowDrawWindow();
}