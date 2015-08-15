#include <oxstd.h>
#include <oxprob.h>
#include <oxfloat.h>
#include <oxdraw.h>
#import <maximize>

fDieboldMariano(const vMeanDensityBayesianForecast, const vDensityFrequentistForecast);
fNeweyWestStandardError(const mX, const vResiduals);

main()
{
	decl v1, v2;

	v1 = loadmat("vAE_OC_RV_GARCH_t.xls");	
	v2 = loadmat("vAE_OC_RV_log_linear_t_RV_RealGARCH.xls");


   	println("RV Benchmark:");
	println("t-GARCH vs log real t-garch RV");
	fDieboldMariano(v1, v2);	
}

fNeweyWestStandardError(const mX, const vResiduals) 
{

	decl iTmp, iT, iK, iM, vIotaT, vOneToT, mAbs_i_minus_j,
	mWeights, mOmegaHat, mInvXX, mCovOLS_estimator, vNW_std_errors;

	iT 	= rows(mX);
	iK 	= columns(mX);

	iM = ceil(4*(iT/100)^(2/9) );
	vIotaT = ones(iT,1);
		
	vOneToT = zeros(iT,1);
	for (iTmp = 0;iTmp<iT; ++iTmp){
		vOneToT[iTmp][0] = iTmp+1;
	}

	mAbs_i_minus_j = fabs(vIotaT*vOneToT'-vOneToT*vIotaT');
	mWeights = (ones(iT,iT) - mAbs_i_minus_j*(1/iM)) .* (mAbs_i_minus_j .<= iM);
	mOmegaHat = (vResiduals * vResiduals') .* mWeights;

	mInvXX = invertsym(mX'*mX);
	mCovOLS_estimator = (iT/(iT-iK)) * mInvXX * (mX'*mOmegaHat*mX) * mInvXX;
	vNW_std_errors = sqrt(diagonal(mCovOLS_estimator)');
	
	return vNW_std_errors;
}

fDieboldMariano(const vForecast1, const vForecast2)
{
	decl vLossDifferential, dMean, dNeweyWestStandardError, dDieboldMarianoStatistic; 
	
	/*--- Define the loss differential (KLIC) ---*/
	vLossDifferential = log(vForecast1) - log(vForecast2);
	
	/*--- Recalculate the HAC variance of the loss differential ---*/
	dMean = meanc(vLossDifferential);
	dNeweyWestStandardError = fNeweyWestStandardError(ones(sizerc(vLossDifferential),1), vLossDifferential - dMean);	  

	/*--- Retrieve the Diebold Mariano statistic ---*/
	dDieboldMarianoStatistic = dMean/dNeweyWestStandardError;
	println("DM-statistic: ", dDieboldMarianoStatistic, "\n\n" );
}

