/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Compute 5-minute and second returns
**
**  Date:
**    	15/01/2015
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
**  Function:	Returns vector with ones where new day starts
**
**  Input: 		vDate: Vector with dates
**
**  Output: 	vNewDay: Vector with ones and zeros. One indicates new day and zero otherwise
*/

fCheckNewDay(const vDate){

	decl iN;
	iN = rows(vDate);

	decl vNewDay;
	vNewDay = zeros(iN,1); //initialize vector
	vNewDay[0] = 1; //Data always starts with a new day
	vNewDay[1:] = (vDate[1:] .!= vDate[:iN-2]); //compare dates

	return vNewDay;
}

/*
**  Function:	Returns vector with returns
**
**  Input: 		vDate: Vector with dates
**
**  Output: 	vNewDay: Vector with returns
*/

fComputeReturnsEveryTrade(const vPrices){

	decl iN, vReturnsEveryTrade;
	iN = rows(vPrices);
	vReturnsEveryTrade = log(vPrices[1:])-log(vPrices[:iN-2]); //compute returns

	return vReturnsEveryTrade;
}

/*
**  Function:	Create matrix where for every second the date and time are created
**
**  Input: 		vDate: Vector with dates
**
**  Output: 	mDatePlusTimeMatrix: Matrix for every second the date and price
*/

fMakeDatePlusTimeMatrix(const vDate){

	decl vUniqueDates, iUniqueDates, iTotalNumberOfSeconds, iNumberOfSecondsInDay;
	vUniqueDates = unique(vDate)';
	iUniqueDates = rows(vUniqueDates);
	iNumberOfSecondsInDay = 60*60*6.5; //these are the number of returns in a certain day
	iTotalNumberOfSeconds = rows(vUniqueDates)*iNumberOfSecondsInDay; // # of days * # of seconds per day
	
	decl vDateForEverySecond;
	vDateForEverySecond = zeros(iTotalNumberOfSeconds+iUniqueDates,1); //vector which indicates for every trade what date it is

	decl i, k;
	i = 0;
	for(k = 0; k < iTotalNumberOfSeconds+iUniqueDates; k = k + (iNumberOfSecondsInDay)+1){
		vDateForEverySecond[k:k+iNumberOfSecondsInDay] = vUniqueDates[i]; //give date for every trade	
		i = i +1;
	}

	decl mDatePlusTimeMatrix;
	mDatePlusTimeMatrix = zeros(iTotalNumberOfSeconds+1*iUniqueDates,7);
	mDatePlusTimeMatrix[][0] =  vDateForEverySecond; //date total
	mDatePlusTimeMatrix[][1] =  floor(vDateForEverySecond/10000);	//year
	mDatePlusTimeMatrix[][2] =	floor((vDateForEverySecond-floor(vDateForEverySecond/10000)*10000)/100); //month
	mDatePlusTimeMatrix[][3] = 	vDateForEverySecond-floor(vDateForEverySecond/100)*100;	//day

	decl vSecond, j;
	vSecond = range(0,59)';
	for(j=1; j < iNumberOfSecondsInDay/60; j = j + 1){
		vSecond = vSecond | range(0,59)'; 
	}
	vSecond = vSecond | 0; //give seconds for every trade

	decl vMinuteRangePerDay, vMinutesPerHour, vMinutesPerDay, vFirstHalfHour;
	vMinuteRangePerDay = range(0,59)';
	vMinutesPerHour = zeros(60*60, 1); 

	decl m, n;
	m = 0; 
	for(n = 0; n < 60*60; n = n + 60){
		vMinutesPerHour[n:n+60-1] = vMinuteRangePerDay[m]; 	
		m = m + 1;
	}

	vFirstHalfHour = vMinutesPerHour[0.5*rows(vMinutesPerHour):];
 	vMinutesPerDay = vMinutesPerHour; //give minutes for every trade 

	for(decl l = 0; l < 5; l++) {
		vMinutesPerDay = vMinutesPerDay | vMinutesPerHour; 
	}

	vMinutesPerDay = vFirstHalfHour | vMinutesPerDay | 0;


	decl vHoursPerDay, vHourRangePerDay;
	vHoursPerDay = zeros(iNumberOfSecondsInDay-30*60,1);
	vHourRangePerDay = range(10,15);
	
	decl r, s;
	r = 0; 
	for(s = 0; s < iNumberOfSecondsInDay-30*60; s = s + 60*60){
		vHoursPerDay[s:s+60*60-1] = vHourRangePerDay[r]; 	
		r = r + 1;
	}

	decl vHourFirstHalfHour;
	vHourFirstHalfHour = 9*ones(30*60,1);
	vHoursPerDay = vHourFirstHalfHour | vHoursPerDay|16;  //give hours for every trade 

	decl mTimePerDay;
	mTimePerDay = vHoursPerDay ~ vMinutesPerDay ~ vSecond; //make out of hours, minutes and seconds a matrix

	decl t;
	for(t = 0; t < rows(vUniqueDates)-1; t = t + 1){
		mTimePerDay = mTimePerDay | (vHoursPerDay ~ vMinutesPerDay ~ vSecond); //put everything below each other
	}

	mDatePlusTimeMatrix[][4:] = mTimePerDay;

	return mDatePlusTimeMatrix;
}

/*
**  Function:	Create matrix where for every second the price is determined
**
**  Input: 		mDateTime: Matrix computed in fMakeDatePlusTimeMatrix
**				mData: Loaded data with information for every trade
**
**  Output: 	mDateTimePrice: mDateTime concatenated with prices for every second
*/

fComputePriceMatrixEverySecond(const mDateTime, const mData){

	decl vDateTrades, vYearTrades, vMonthTrades, vDayTrades, vHoursTrades, vMinutesTrades, vSecondsTrades;
	vDateTrades = mData[][0];
	vYearTrades =  floor(vDateTrades/10000);	//year
	vMonthTrades = floor((vDateTrades-floor(vDateTrades/10000)*10000)/100); //month
	vDayTrades = vDateTrades-floor(vDateTrades/100)*100;	//day
	vHoursTrades = mData[][1]; //hours
	vMinutesTrades = mData[][2]; //minutes
	vSecondsTrades = mData[][3]; //seconds

	decl vYearEverySecond, vMonthEverySecond, vDayEverySecond, vHoursEverySecond, vMinutesEverySecond, vSecondsEverySecond;
	vYearEverySecond = mDateTime[][1];	//year
	vMonthEverySecond = mDateTime[][2];	//month
	vDayEverySecond = mDateTime[][3];	//day
	vHoursEverySecond = mDateTime[][4]; //hours
	vMinutesEverySecond = mDateTime[][5]; //minutes
	vSecondsEverySecond = mDateTime[][6]; //seconds

	decl vPricesTrades, vPricesEverySecond;
	vPricesTrades = mData[][4];
	vPricesEverySecond = zeros(rows(vYearEverySecond),1);
	vPricesEverySecond[0] = vPricesTrades[0]; 

	decl i, j;
	j = 1;
	for(i = 1; i < rows(vYearEverySecond); i = i + 1){
		if( j<rows(vPricesTrades)){
		//check here whether on a certain time if there is a trade and put price in a vector
			if(vYearEverySecond[i] == vYearTrades[j] && vMonthEverySecond[i] == vMonthTrades[j] && vDayEverySecond[i] == vDayTrades[j] && vHoursEverySecond[i] == vHoursTrades[j] && vMinutesEverySecond[i] == vMinutesTrades[j] && vSecondsEverySecond[i] == vSecondsTrades[j]){
				vPricesEverySecond[i] = vPricesTrades[j];
				j = j + 1;
			}else{	//if there is no trade on a given second, put previous price for that second
				vPricesEverySecond[i] = vPricesEverySecond[i-1];
			}
		}else{ //last trade, last day
		vPricesEverySecond[i] = vPricesEverySecond[i-1];
		}
	}
	
	decl mDateTimePrice;
	mDateTimePrice = mDateTime ~ vPricesEverySecond; //make new matrix which includs prices

	return mDateTimePrice;
}

/*
**  Function:	Extract prices from mDateTimePrice
**
**  Input: 		mDateTimePrice: Matrix computed in fComputePriceMatrixEverySecond
**
**  Output: 	mPriceEverySecond: Matrix with prices where on columns are days and each row is a new second
*/

fComputePricesEverySecond(const mDateTimePrice){

	decl iUniqueDays, iNumberOfSecondsPerDay, mPriceEverySecond;
	iUniqueDays = columns(unique(mDateTimePrice[][0]));
	iNumberOfSecondsPerDay = 6.5*60*60+1; //number of possible prices on a day
	mPriceEverySecond = zeros(iNumberOfSecondsPerDay, iUniqueDays); 
	mPriceEverySecond = shape(mDateTimePrice[][7],iNumberOfSecondsPerDay,iUniqueDays); //each day is a column, each second is a row

	return mPriceEverySecond;
}

/*
**  Function:	Compute return for each second on the day
**
**  Input: 		mDateTimePrice: Matrix computed in fComputePriceMatrixEverySecond
**
**  Output: 	mReturnEverySecond: Matrix with returns where on columns are days and each row is a new second
*/

fComputeReturnEverySecond(const mDateTimePrice){

	decl iUniqueDays, iNumberOfSecondsPerDay, mPriceEverySecond;
	iUniqueDays = columns(unique(mDateTimePrice[][0]));
	iNumberOfSecondsPerDay = 6.5*60*60+1; //possible number of prices on a day
	mPriceEverySecond = zeros(iNumberOfSecondsPerDay, iUniqueDays); 
	mPriceEverySecond = shape(mDateTimePrice[][7],iNumberOfSecondsPerDay,iUniqueDays);	

	decl mReturnEverySecond, i;
	mReturnEverySecond = zeros(iNumberOfSecondsPerDay-1, iUniqueDays); //1 observation lost
	for(i = 0; i < iUniqueDays; i = i + 1){
			mReturnEverySecond[][i] = log(mPriceEverySecond[1:][i])-log(mPriceEverySecond[:iNumberOfSecondsPerDay-2][i]); //compute returns
	}

	return mReturnEverySecond;
}

/*
**  Function:	Compute return for each ... minutes on the day (for example 5)
**
**  Input: 		mDateTimePrice: Matrix computed in fComputePriceMatrixEverySecond
** 				iHowManyMinutesReturns: Indicates how many minutes returns you want
**
**  Output: 	mReturns: Matrix with returns for ... minutes
*/

fComputeDotsMinutesReturnsPerDay(const mDateTimePrice, const iHowManyMinutesReturns){

	decl iNumberOfSecondsInMinutesReturns, iNumberOfDays, iNumberOfPricesPerHour, iPricesPerday;
	iNumberOfSecondsInMinutesReturns = iHowManyMinutesReturns * 60;
	iNumberOfDays = columns(unique(mDateTimePrice[][0]));
	iNumberOfPricesPerHour = 60/iHowManyMinutesReturns; 
	iPricesPerday = 6.5 * iNumberOfPricesPerHour+1;

	decl vPricesTotal, mPricesTotal;
	vPricesTotal = mDateTimePrice[][7];
	mPricesTotal =	shape(vPricesTotal, rows(vPricesTotal)/iNumberOfDays, iNumberOfDays);

	decl s, t, mPrices;
	mPrices = zeros(iPricesPerday,iNumberOfDays);
	for(s = 0; s < iNumberOfDays; s = s + 1){
		for(t = 0; t < iPricesPerday; t = t + 1){
			mPrices[t][s] = mPricesTotal[t*iNumberOfSecondsInMinutesReturns][s];
		}
	}
	decl mReturns, i;
	mReturns = zeros(iPricesPerday-1, iNumberOfDays); //lose 1 observation a day

	for(i = 0; i < iNumberOfDays; i = i + 1){
			mReturns[][i] = log(mPrices[1:][i])-log(mPrices[:iPricesPerday-2][i]); //compute returns
	}
	
    return mReturns;
}

main(){
	decl mData, iN;
	mData = loadmat("2008_cleaned.csv");
	iN = rows(mData);	 

	decl vPrices, vReturnsEveryTrade;
	vPrices = mData[][4];
	vReturnsEveryTrade = fComputeReturnsEveryTrade(vPrices);
	savemat("PricesEveryTrade2008.csv", mData[][0]~vPrices); //include overnight returns
   
	decl vDate, mDateTime;
	vDate = mData[][0];
	mDateTime = fMakeDatePlusTimeMatrix(vDate);

	decl mDateTimePrice, mPricesEverySecond;
	mDateTimePrice = fComputePriceMatrixEverySecond(mDateTime, mData);
	mPricesEverySecond = fComputePricesEverySecond(mDateTimePrice); 
	savemat("SecondPrices2008.csv", mPricesEverySecond);

	decl iHowManyMinutesReturns, mReturns;
	iHowManyMinutesReturns = 5;
	mReturns = fComputeDotsMinutesReturnsPerDay(mDateTimePrice, iHowManyMinutesReturns); //compute five minute returns
	savemat("FiveMinuteReturns2008.csv", mReturns);
}