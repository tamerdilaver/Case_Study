/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Compute omega^2 as described in Barndorff-Nielsen et al. (2009)
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

fComputeRV(const mReturns){

	decl vRV;
	vRV = sumsqrc(mReturns);
	
	return vRV;
}

main(){
	decl mData, vPrices, vLogPrices, vDates;
	mData = loadmat("PricesEveryTrade_all.csv"); //load prices and dates of every trade
	vDates = mData[][0];
	vPrices = mData[][1];
	vLogPrices = log(vPrices); //take log prices

	decl vRVDense, iNumberOfTotalTradingDays;
	iNumberOfTotalTradingDays = columns(unique(vDates)); //load number of trading days (here:1699)
	vRVDense = zeros(iNumberOfTotalTradingDays,1); 

	decl iTotalNumberOfPricesAllDays, vNewTradingDay, vIndexNumbersNewTradingDay;
	iTotalNumberOfPricesAllDays = rows(vPrices);
	vNewTradingDay = 1|(vDates[1:] .!= vDates[:iTotalNumberOfPricesAllDays-2])|1; //vector with ones for each trading day
	vIndexNumbersNewTradingDay = vecrindex(vNewTradingDay); // compute indices of new trading days

	decl iQ, mTimeSeries, iNumberOfTradesOnDayT;
	iQ = 25; //take each iQth trade
	
	decl t, i, j, vRV, vOmegaSquared;
	vRV = zeros(iQ,1);
	vRVDense = zeros(iNumberOfTotalTradingDays,1);
	vOmegaSquared = zeros(iNumberOfTotalTradingDays,1);
	
	for(t = 0; t < iNumberOfTotalTradingDays; t++){
		iNumberOfTradesOnDayT = vIndexNumbersNewTradingDay[t+1]-vIndexNumbersNewTradingDay[t];
		mTimeSeries = zeros(iQ,floor(iNumberOfTradesOnDayT/iQ));
		for(i = 0; i < iQ; i++){
			for(j = 0; j < iNumberOfTradesOnDayT-iQ-i; j = j + iQ){	
			mTimeSeries[i][j/iQ] = 100 * (vLogPrices[j+i+iQ+vIndexNumbersNewTradingDay[t]]-vLogPrices[j+i+vIndexNumbersNewTradingDay[t]]);
			}
		vRV[i] = fComputeRV((mTimeSeries[i][])');	
		}
	vRVDense[t] = meanc(vRV);
	vOmegaSquared[t] = vRVDense[t]/(2*iNumberOfTradesOnDayT);
	}
	savemat("OmegaSquared.csv",vOmegaSquared);
}