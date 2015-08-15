/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**	Compute bipower variation and realized variance
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

/*
**  Function:	Compute bipower variation
**
**  Input: 		mReturns: Matrix with returns 
**
**  Output: 	vBV: Vector with bipower variations, computed for each column (=day)
*/

fComputeBipowerVariation(const mReturns){

	decl dMu, mRhelp, iNumberOfDays, iRows;
	dMu = sqrt(2/M_PI); 
	iNumberOfDays = columns(mReturns);
	iRows = rows(mReturns);
	mRhelp = zeros(iRows,iNumberOfDays);
	vBV= zeros(iNumberOfDays);

	decl i, j;
	for(i=0; i<columns(mReturns); i++){
		for (j=1; j<rows(mReturns); j++){
			mRhelp[j][i] = fabs(mReturns[j-1][i]) * fabs(mReturns[j][i]);	
		}
	}

	decl vBV;
	vBV = sumc(mRhelp[1:][])/dMu^2;
	
	return vBV;
}

main(){
	decl mFiveMinuteReturns;
	mFiveMinuteReturns = loadmat("FiveMinuteReturns_All.csv");
	mFiveMinuteReturns = 100 * mFiveMinuteReturns;

	decl vRV; //realized variance
	vRV = sumsqrc(mFiveMinuteReturns);	//compute RV

	decl mDataCloseToClose;
	mDataCloseToClose = loadmat("ReturnsCloseToClose.csv");

	decl vDate, vYear, vMonth, vDay;
	vDate = mDataCloseToClose[][0];
	vYear 	= floor(vDate/10000); //compute year from vDate							
	vMonth 	= floor((vDate-floor(vDate/10000)*10000)/100); //compute month from vDate		
	vDay 	= vDate-floor(vDate/100)*100; //compute day from vDate	
	vDate = dayofcalendar(vYear, vMonth, vDay); //date that Ox understands

	decl vBV; //bipower variation
	vBV = fComputeBipowerVariation(mFiveMinuteReturns); //compute BV
}