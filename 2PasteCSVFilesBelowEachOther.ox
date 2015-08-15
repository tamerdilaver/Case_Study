/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**	Paste all cleaned date files below each other and save in one CSV file
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

main(){
	decl mData2008, mData2009, mData2010, mData2011, mData2012, mData2013, mData2014;

	//load data
	mData2008 = loadmat("2008_Cleaned.csv");
	mData2009 = loadmat("2009_Cleaned.csv");
	mData2010 = loadmat("2010_Cleaned.csv");
	mData2011 = loadmat("2011_Cleaned.csv");
	mData2012 = loadmat("2012_Cleaned.csv");
	mData2013 = loadmat("2013_Cleaned.csv");
	mData2014 = loadmat("2014_Cleaned.csv");

	//save data below each other
	savemat("Cleaned_all.csv", mData2008| mData2009| mData2010| mData2011| mData2012| mData2013| mData2014);
}