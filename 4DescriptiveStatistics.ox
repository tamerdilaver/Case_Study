/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Compute some useful statistics
**
**  Date:
**    	12/01/2015
**
**  Author:
**	  	Tamer Dilaver, Koen de Man & Sina Zolnoor
**
**	Supervisor:
**		L.N. Hoogerheide & S.J. Koopman
**
*/
#include <oxstd.h>
#include <oxfloat.h>

/*
**  Function:	Print first four moments, maximum and minimum
**
**  Input: 		vReturns: Vector with (daily) returns
*/

fDescriptiveStatistics(const vReturns){

	println("Moments: ", moments(vReturns));  //give number of observations and first four moments
	println("Maximum: ", maxc(vReturns));
	println("Minimum: ", minc(vReturns));
}


main(){

	decl mDataCloseToClose, mDataOpenToClose;
	mDataCloseToClose = loadmat("ReturnsCloseToClose.csv");	 //load return data
	mDataOpenToClose = loadmat("ReturnsOpenToClose.csv");

	decl vDate, vYear, vMonth, vDay;
	vDate = mDataCloseToClose[][0];
	vYear 	= floor(vDate/10000); //extract year from date number								
	vMonth 	= floor((vDate-floor(vDate/10000)*10000)/100); //extract month from date number		
	vDay 	= vDate-floor(vDate/100)*100; //extract day from date number	
	vDate = dayofcalendar(vYear, vMonth, vDay); //create vector with dates that Ox can handle

	decl vReturnsCloseToClose, vReturnsOpenToClose;
	vReturnsCloseToClose = mDataCloseToClose[][1];
	vReturnsOpenToClose = mDataOpenToClose[][1];
	
	println("Descriptive statistics for Close-to-Close returns: ");
	fDescriptiveStatistics(100*vReturnsCloseToClose);

	println("Descriptive statistics for Open-to-Close returns: ");
	fDescriptiveStatistics(100*vReturnsOpenToClose);

//	DrawTMatrix(1, 100*vReturnsCloseToClose', {"Close-to-close"}, vDate');
//	DrawTMatrix(0, 100*vReturnsOpenToClose', {"Open-to-close"}, vDate');
//	ShowDrawWindow();
}