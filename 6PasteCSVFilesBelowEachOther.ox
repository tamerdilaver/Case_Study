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
#include <oxdraw.h>

main(){
	decl mData2008, mData2009, mData2010, mData2011, mData2012, mData2013, mData2014;

	mData2008 = loadmat("FiveMinuteReturns2008.csv");
	mData2009 = loadmat("FiveMinuteReturns2009.csv");
	mData2010 = loadmat("FiveMinuteReturns2010.csv");
	mData2011 = loadmat("FiveMinuteReturns2011.csv");
	mData2012 = loadmat("FiveMinuteReturns2012.csv");
	mData2013 = loadmat("FiveMinuteReturns2013.csv");
	mData2014 = loadmat("FiveMinuteReturns2014.csv");
	savemat("FiveMinuteReturns_all.csv", mData2008 ~ mData2009 ~ mData2010 ~ mData2011 ~ mData2012 ~ mData2013 ~ mData2014);
	

	mData2008 = loadmat("SecondPrices2008.csv");
	mData2009 = loadmat("SecondPrices2009.csv");
	mData2010 = loadmat("SecondPrices2010.csv");
	mData2011 = loadmat("SecondPrices2011.csv");
	mData2012 = loadmat("SecondPrices2012.csv");
	mData2013 = loadmat("SecondPrices2013.csv");
	mData2014 = loadmat("SecondPrices2014.csv");
	savemat("SecondPrices_all.csv", mData2008 ~ mData2009 ~ mData2010 ~ mData2011 ~ mData2012 ~ mData2013 ~ mData2014);
   
	mData2008 = loadmat("ReturnsEveryTrade2008.csv");
	mData2009 = loadmat("ReturnsEveryTrade2009.csv");
	mData2010 = loadmat("ReturnsEveryTrade2010.csv");
	mData2011 = loadmat("ReturnsEveryTrade2011.csv");
	mData2012 = loadmat("ReturnsEveryTrade2012.csv");
	mData2013 = loadmat("ReturnsEveryTrade2013.csv");
	mData2014 = loadmat("ReturnsEveryTrade2014.csv");
	savemat("ReturnsEveryTrade_all.csv", mData2008 | mData2009 | mData2010 | mData2011 | mData2012 | mData2013 | mData2014);

	mData2008 = loadmat("PricesEveryTrade2008.csv");
	mData2009 = loadmat("PricesEveryTrade2009.csv");
	mData2010 = loadmat("PricesEveryTrade2010.csv");
	mData2011 = loadmat("PricesEveryTrade2011.csv");
	mData2012 = loadmat("PricesEveryTrade2012.csv");
	mData2013 = loadmat("PricesEveryTrade2013.csv");
	mData2014 = loadmat("PricesEveryTrade2014.csv");
	savemat("PricesEveryTrade_all.csv", mData2008 | mData2009 | mData2010 | mData2011 | mData2012 | mData2013 | mData2014);
}