/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**	Compute realized kernel 
**
**  Date:
**    	16/01/2015
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

main(){

	decl vRVSparse, vOmegaSquared; //load previously computed matrices
	vRVSparse = loadmat("RVSparse.csv");
    vOmegaSquared = loadmat("OmegaSquared.csv");

	decl mData, vPrices, vLogPrices, vDates;
	mData = loadmat("PricesEveryTrade_all.csv"); //load data
	vDates = mData[][0];
	vPrices = mData[][1];
	vLogPrices = log(vPrices); //take log prices

	decl iNumberOfTotalTradingDays;
	iNumberOfTotalTradingDays = columns(unique(vDates)); //extract all trading days 
	
	decl iTotalNumberOfPrices, vReturns;
	iTotalNumberOfPrices = rows(vLogPrices);
	
	decl vKsiSquared;
	vKsiSquared = vOmegaSquared ./ vRVSparse; //compute ksi

	decl vNewTradingDay, vIndexNumbersNewTradingDay, vNumberOfTradesPerDay, k;
	vNewTradingDay = 1|(vDates[1:] .!= vDates[:iTotalNumberOfPrices-2])|1;	//create vector with new trading days
	vIndexNumbersNewTradingDay = vecrindex(vNewTradingDay);
	vNumberOfTradesPerDay = zeros(iNumberOfTotalTradingDays,1);

	for(k = 0; k < iNumberOfTotalTradingDays; k++){
		vNumberOfTradesPerDay[k] =	vIndexNumbersNewTradingDay[k+1]-vIndexNumbersNewTradingDay[k];	//compute number of trades per day
	}
	
	decl vHStar, dCStar;
	dCStar = 3.5134; //optimal c 
	vHStar = ceil(dCStar * (vKsiSquared .^(2/5) .* vNumberOfTradesPerDay .^(3/5))); //compute optimal H
	
	decl t, h, vRK, vLogPricesForDayT, vReturnsForDayT;
	vRK = zeros(iNumberOfTotalTradingDays,1);
	for(t = 0; t < iNumberOfTotalTradingDays; t++){	//compute RK for each day
		vLogPricesForDayT = vLogPrices[vIndexNumbersNewTradingDay[t]:(vIndexNumbersNewTradingDay[t+1]-1)];	 //compute log prices
		vReturnsForDayT = 100*(vLogPricesForDayT[1:] - vLogPricesForDayT[:rows(vLogPricesForDayT)-2]);	//compute returns
		for(h = -vHStar[t]; h <= vHStar[t]; h++){
			if(vNumberOfTradesPerDay[t] > vHStar[t]+1){ //in practice this is always satisfied
				vRK[t] = vRK[t] + fComputeParzenKernel(h/(vHStar[t]+1)) * fComputeGamma(vReturnsForDayT, h);
			}
		}
	}
	savemat("RK.csv", vRK);	//save results
}