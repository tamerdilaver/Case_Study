/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Compute RV_sparse as described in Barndorff-Nielsen et al. (2009) for estimation of IV
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

/*
**  Function:	Compute realized variance
**
**  Input: 		vReturns: Matrix with returns
**
**  Output: 	vRV: Vector with realized variances for each column (=day)
*/

fComputeRV(const vReturns){

	decl vRV;
	vRV = sumsqrc(vReturns);
	
	return vRV;
}

main(){
	decl mPrices, mLogPrices;
	mPrices = loadmat("SecondPrices_all.csv");	 //load data
	mLogPrices = log(mPrices);

	decl iNumberOfTradingDays, iNumberOfObservations, iQ;
	iNumberOfTradingDays = columns(mLogPrices);
	iNumberOfObservations = rows(mLogPrices);
	iQ = 1200; //number of time series per day

	decl t, i, j, mTimeSeries, vRV, vRVSparse;
	vRV = zeros(iQ,1);
	mTimeSeries = zeros(iQ,floor(iNumberOfObservations/iQ));
	vRVSparse = zeros(iNumberOfTradingDays,1);
	for(t = 0; t < iNumberOfTradingDays; t++){	//compute RVSparse for each day
		for(i = 0; i < iQ; i++){ //total of iQ time series
			for(j = 0; j < iNumberOfObservations-iQ-i; j = j + iQ){ //take steps of 1200 in picking the returns
				mTimeSeries[i][j/iQ] = 100*(mLogPrices[j+i+iQ][t]-mLogPrices[j+i][t]); //put returns in matrix
			}
		vRV[i] = fComputeRV((mTimeSeries[i][])');//compute RV for each of the 1200 time series	
		}
	vRVSparse[t] = meanc(vRV); //take mean of RVs for each day
	}
	savemat("RVSparse.csv",vRVSparse);
}