/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Simulate HF data and compute RV_sec, RV_5min and RK for simulated data
**
**  Date:
**    	18/01/2015
**
**  Author:
**	  	Tamer Dilaver, Koen de Man & Sina Zolnoor
**
**	Supervisor:
**		L.H. Hoogerheide & S.J. Koopman
**
*/

#include <oxstd.h>
#include <oxfloat.h>
#include <oxdraw.h>
#include <oxprob.h>

/*
**  Function:	Simulate stock prices with normal distribution for the noise
**
**  Input: 		iNumberOfSecondsInDay: integer with number of seconds in day
**				iNumberOfDays: integer with number of days
**				dSigmaSquared_U: double, with sigma^2 of noise U
**				dSigmaSquared_Z: double, with sigma^2 of random walk Z_t
**
**  Output: 	mPrices: Matrix with stock prices per second
*/

fSimulateHFPricesNormal(const iNumberOfSecondsInDay, const iNumberOfDays, const dSigmaSquared_U, const dSigmaSquared_Z){

	decl mZ, mU, dYZero;
	mZ = sqrt(dSigmaSquared_Z)*rann(iNumberOfSecondsInDay, iNumberOfDays);
	mU = sqrt(dSigmaSquared_U)*rann(iNumberOfSecondsInDay, iNumberOfDays);
	dYZero = 100; //`clean' price to start with
	
	decl mPrices, mY;
	mY = zeros(iNumberOfSecondsInDay, iNumberOfDays);
	mY[0][0] = dYZero;
	mPrices = zeros(iNumberOfSecondsInDay, iNumberOfDays);
	mPrices[0][0] = mY[0][0] + mU[0][0]; //start price
	
	decl i,j;
	for(i=0; i<iNumberOfDays; i++){
		if(i>0){ //check whether we are on a new day
			mY[0][i] = mY[iNumberOfSecondsInDay-1][i-1];  //no overnight jumps, hence take last (log) price of previous data
			mPrices[0][i] = mY[0][i] + mU[0][i];
		}
		for(j=1; j<iNumberOfSecondsInDay; j++){
			mY[j][i] = mY[j-1][i]*exp(-0.5*dSigmaSquared_Z+mZ[j][i]); //Geometric Brownian motion
			mPrices[j][i] = mY[j][i] + mU[j][i]; //`dirty' (log) price
		}
	}
	return mPrices;
}

/*
**  Function:	Simulate stock prices with Student's-t distribution for the noise
**
**  Input: 		iNumberOfSecondsInDay: integer with number of seconds in day
**				iNumberOfDays: integer with number of days
**				dSigmaSquared_U: double, with sigma^2 of noise U
**				dSigmaSquared_Z: double, with sigma^2 of random walk Z_t
**				dNu: double, with degrees of freedom
**
**  Output: 	mPrices: Matrix with stock prices per second
*/

fSimulateHFPricesStudentT(const iNumberOfSecondsInDay, const iNumberOfDays, const dSigmaSquared_U, const dSigmaSquared_Z, const dNu){

	decl mZ, mU, dYZero;
	mZ = sqrt(dSigmaSquared_Z)*rann(iNumberOfSecondsInDay, iNumberOfDays);
	mU = sqrt(dSigmaSquared_U)*sqrt((dNu-2)/dNu)*rant(iNumberOfSecondsInDay, iNumberOfDays, dNu);
	dYZero = 100; //price to start with
	
	decl mPrices, mY;
	mY = zeros(iNumberOfSecondsInDay, iNumberOfDays);
	mY[0][0] = dYZero;
	mPrices = zeros(iNumberOfSecondsInDay, iNumberOfDays);
	mPrices[0][0] = mY[0][0] + mU[0][0]; //start price
	
	decl i,j;
	for(i=0; i<iNumberOfDays; i++){
		if(i>0){   //check whether we are on a new day
			mY[0][i] = mY[iNumberOfSecondsInDay-1][i-1]; //no overnight jumps, hence take last (log) price of previous data
			mPrices[0][i] = mY[0][i] + mU[0][i];
		}
		for(j=1; j<iNumberOfSecondsInDay; j++){
			mY[j][i] = mY[j-1][i]*exp(-0.5*dSigmaSquared_Z+mZ[j][i]); //Geometric Brownian motion
			mPrices[j][i] = mY[j][i] + mU[j][i]; //`dirty' (log) price
		}
	}
	return mPrices;
}

/*
**  Function:	Compute return for each second
**
**  Input: 		mPrices: Matrix with prices			
**
**  Output: 	mReturns: Matrix with seconds returns
*/

fComputeSecondReturns(const mPrices){

	decl iPricesPerDay, iNumberOfDays;
	iPricesPerDay = rows(mPrices);
	iNumberOfDays = columns(mPrices);

	decl i, mReturns;
	mReturns = zeros(iPricesPerDay-1, iNumberOfDays);
	for(i=0; i<iNumberOfDays; i++){
		mReturns[][i] = log(mPrices[1:][i])-log(mPrices[:iPricesPerDay-2][i]); 	 //compute return
	}
	return mReturns;
}

/*
**  Function:	Compute five minute returns
**
**  Input: 		mPrices: Matrix with prices			
**
**  Output: 	mFiveMinuteReturns: Matrix with five minute returns
*/

fComputeFiveMinuteReturns(const mPrices){

	decl iNumberOfPrices, iNumberOf5MinutePrices, iNumberOfDays;
	iNumberOfPrices = rows(mPrices); //equal to total number of seconds in day
	iNumberOf5MinutePrices = (iNumberOfPrices/300)+1; //on a certain day, this is how many five minute prices there are
	iNumberOfDays = columns(mPrices); //number of day

	decl vHelp, vNewFiveMinutes, i;
	vHelp = 1|zeros(300-1,1);
	vNewFiveMinutes = vHelp;
	for(i=1; i<iNumberOf5MinutePrices-1; i++){
		vNewFiveMinutes = vNewFiveMinutes | vHelp;   
	}
	vNewFiveMinutes[iNumberOfPrices-1] = 1;
	
	decl mFiveMinutePrices;
	mFiveMinutePrices = zeros(iNumberOf5MinutePrices, iNumberOfDays);
	mFiveMinutePrices = selectifr(mPrices, vNewFiveMinutes);

	decl mFiveMinuteReturns, j;
	mFiveMinuteReturns = zeros(iNumberOf5MinutePrices-1, iNumberOfDays);
	for(j=0; j<iNumberOfDays; j++){
		mFiveMinuteReturns[][j] = log(mFiveMinutePrices[1:][j])-log(mFiveMinutePrices[:iNumberOf5MinutePrices-2][j]); //compute returns
	}
	
	return mFiveMinuteReturns;
}

/*
**  Function:	Compute realized variance
**
**  Input: 		vReturns: Matrix with returns
**
**  Output: 	vRV: Vector with realized variances for each column (=day)
*/

fComputeRV(const mReturns){

	decl vRV;
	vRV = sumsqrc(mReturns);
	
	return vRV;
}

/*
**  Function:	Compute Parzen kernel
**
**  Input: 		dX: Double, argument of Parzen function f(.)			
**
**  Output: 	dK: Double, dK = f(dX)
*/

fComputeParzenKernel(const dX){

	decl dK;

	if(fabs(dX) >= 0 &&fabs(dX) <= 0.5){
		dK = 1-6*dX^2 + 6 * fabs(dX)^3;
	}else if(fabs(dX) > 0.5 && fabs(dX) <= 1){
		dK = 2 * (1 - fabs(dX))^3;
	}else{
		dK = 0;
	}
	return dK;
}

/*
**  Function:	Compute gamma function
**
**  Input: 		vReturns: Vector of Returns
**				iH: Integer, Range of summation 
**
**  Output: 	dGamma_iH: Double, output of funtion gamma_H 
*/

fComputeGamma(const vReturns, const iH){

	decl iNumberOfTrades, dGamma_iH; 
	iNumberOfTrades = rows(vReturns);
	dGamma_iH = sumc(vReturns[:iNumberOfTrades-fabs(iH)-1] .* vReturns[fabs(iH):]);
	
	return dGamma_iH;                 	
}

/*
**  Function:	Compute omega^2	using RV_dense
**
**  Input: 		mPrices: Matrix with prices for each second and each day
**
**  Output: 	vOmegaSquared: vector of omega^2 for each day
*/

fComputeOmegaSquared(const mPrices){
	decl iNumberOfTradingDays, iNumberOfObservations, iQ, mLogPrices;
	mLogPrices = log(mPrices);
	iNumberOfTradingDays = columns(mLogPrices);
	iNumberOfObservations = rows(mLogPrices);
	iQ = 25; //number of time series

	decl t, i, j, mTimeSeries, vRV, vRVDense, vOmegaSquared;
	vRV = zeros(iQ,1);
	mTimeSeries = zeros(iQ,floor(iNumberOfObservations/iQ));
	vRVDense = zeros(iNumberOfTradingDays,1);
	vOmegaSquared = zeros(iNumberOfTradingDays,1);
	for(t = 0; t < iNumberOfTradingDays; t++){
		for(i = 0; i < iQ; i++){
			for(j = 0; j < iNumberOfObservations-iQ-i; j = j + iQ){
				mTimeSeries[i][j/iQ] = 100 * (mLogPrices[j+i+iQ][t]-mLogPrices[j+i][t]);
			}
		vRV[i] = fComputeRV((mTimeSeries[i][])');	
		}
	vRVDense[t] = meanc(vRV);
	vOmegaSquared[t] = vRVDense[t]/(2*iNumberOfObservations);
	}
	return vOmegaSquared;
}

/*
**  Function:	Compute RV_SPARSE
**
**  Input: 		mPrices: Matrix with prices for each second and each day
**
**  Output: 	vRVSparse: vector of estimated IV for each day
*/

fComputeRVSparse(const mPrices){

	decl iNumberOfTradingDays, iNumberOfObservations, iQ, mLogPrices;
	mLogPrices = log(mPrices);
	iNumberOfTradingDays = columns(mLogPrices);
	iNumberOfObservations = rows(mLogPrices);
	iQ = 1200; //number of time series

	decl t, i, j, mTimeSeries, vRV, vRVSparse;
	vRV = zeros(iQ,1);
	mTimeSeries = zeros(iQ,floor(iNumberOfObservations/iQ));
	vRVSparse = zeros(iNumberOfTradingDays,1);
	for(t = 0; t < iNumberOfTradingDays; t++){
		for(i = 0; i < iQ; i++){
			for(j = 0; j < iNumberOfObservations-iQ-i; j = j + iQ){
				mTimeSeries[i][j/iQ] = 100 * (mLogPrices[j+i+iQ][t]-mLogPrices[j+i][t]);
			}
		vRV[i] = fComputeRV((mTimeSeries[i][])');	
		}
	vRVSparse[t] = meanc(vRV);
	}
	return vRVSparse;
}

/*
**  Function:	Compute realized kernel (RK)
**
**  Input: 		mPrices: Matrix with prices for each second and each day
**
**  Output: 	vKV: Vector with kernel variance (=realized kernel)	for each day
*/

fComputeRK(const mPrices){

	decl iNumberOfTotalTradingDays;
	iNumberOfTotalTradingDays = columns(mPrices);
	
	decl iTotalNumberOfPricesPerDay, mReturns;
	iTotalNumberOfPricesPerDay = rows(mPrices);
	mReturns = 100 * (log(mPrices[1:][])-log(mPrices[:iTotalNumberOfPricesPerDay-2][]));

	decl vOmegaSquared, vRVSparse;
	vOmegaSquared = fComputeOmegaSquared(mPrices);
	vRVSparse = fComputeRVSparse(mPrices);
	
	decl vKsiSquared;
	vKsiSquared = vOmegaSquared ./ vRVSparse;
	
	decl vHStar, dCStar, vNumberOfTradesPerDay;
	vNumberOfTradesPerDay =  iTotalNumberOfPricesPerDay * ones(iNumberOfTotalTradingDays,1);
	dCStar = 3.5134; 
	vHStar = ceil(dCStar * (vKsiSquared .^(2/5) .* vNumberOfTradesPerDay .^(3/5)));
	
	decl t, h, vKV, vReturnsForDayT;
	vKV = zeros(iNumberOfTotalTradingDays,1);
	for(t = 0; t < iNumberOfTotalTradingDays; t++){
		vReturnsForDayT = mReturns[][t];	
		for(h = -vHStar[t]; h <= vHStar[t]; h++){
			if(vNumberOfTradesPerDay[t] > vHStar[t]+1){
				vKV[t] = vKV[t] + fComputeParzenKernel(h/(vHStar[t]+1)) * fComputeGamma(vReturnsForDayT, h);
			}
		}
	}
	return vKV;
}

main(){

	decl iNumberOfDays, dSigmaSquared_U, iNumberOfSecondsInDay, dSigmaSquared_Z, dNu; 
	iNumberOfDays = 2500; //Number of days: T
	dSigmaSquared_U = 0.0001^2; //vary with this number
	iNumberOfSecondsInDay = 6.5*60*60;
	dSigmaSquared_Z = 0.01^2/iNumberOfSecondsInDay;	//keep this constant
	dNu = 5; //degrees of freedom
	
	decl mPrices, mReturns, mFiveMinuteReturns;  
	mPrices = fSimulateHFPricesNormal(iNumberOfSecondsInDay, iNumberOfDays, dSigmaSquared_U, dSigmaSquared_Z); //draw from normal
//	mPrices	= fSimulateHFPricesStudentT(iNumberOfSecondsInDay, iNumberOfDays, dSigmaSquared_U, dSigmaSquared_Z, dNu); //draw from Student's-t
	mReturns = 100*fComputeSecondReturns(mPrices);
	mFiveMinuteReturns = 100*fComputeFiveMinuteReturns(mPrices);

	decl vRV;
	vRV = fComputeRV(mFiveMinuteReturns); //compute RV for 5 minute returns
	println(moments(vRV'));
	vRV = fComputeRV(mReturns);	//compute RV for second returns
	println(moments(vRV'));
	
	decl vKV;	
	vKV = fComputeRK(mPrices); 	//compute RK for second returns
	println(moments(vKV));	
}