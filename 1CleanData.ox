/*
**	Case Study Financial Econometrics 4.3 
**
**  Purpose:
**  	Clean data using the cleaning and filtering steps as described in Barndorff-Nielsen et al. (2009)
**
**  Date:
**    	12/01/2015
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
**  Function:	Replace a character in a string with another character
**
**  Input: 		aString: String
**				sCharToBeReplaced: Indicates the 'old' character that you wish to replace by another character
**				sCharReplaced = Indicates the 'new' character
**
**  Output: 	aString: Same as input aString, except character that needs to be replaced is replaced.
*/

fReplaceCharInString(const asString, const sCharToBeReplaced, const sCharReplaced){

	decl sOutput, iLengthOfString; 
	sOutput = asString[0]; //Initialize output of string
	iLengthOfString = sizerc(sOutput); //Count total number of characters in string

	decl i;
	for (i = 0; i < iLengthOfString; i++){ 
		if(sOutput[i] == sCharToBeReplaced){ //check whether the character you wish to replace is present in the string and on which location
			sOutput[i]= sCharReplaced; //replace character
		}
	}
	
	asString[0] = sOutput; //return string with possible replacement of 'old' character by 'new' character
}

/*
**  Function:	Read Data file and perform cleaning step 2 (not traded on same exchange) and 3 (other sale condition)
** 				as described in report.
**
*/

fReadingStringData(){

	decl mData, iN;
	mData = loadmat("2014.csv"); //Load data per year to count number of original rows available in 'dirty' data file
	iN = rows(mData);

	decl sNameFileIn, sNameFileOut;	//determine the names of the file we want to read in (sNameFileIn) and the file to write the cleaned data in (sNameFileOut)
	sNameFileIn = "2014.csv";
	sNameFileOut = "2014_Cleaned.csv";

	decl sFileIn, sFileOut; 
	sFileIn = fopen(sNameFileIn, "r"); //Open file to read from
	sFileOut = fopen(sNameFileOut, "w"); //Open file to write in

	decl iTellerQ, iTellerSaleCondition;
	iTellerQ = 0; //counts number of observations on the 'right' exchange 
	iTellerSaleCondition = 0; //counts number of observations with 'right' sale condition

	decl sString, iLengthOfString;

	while(!(feof(sFileIn))){ //read as long as there are lines to read
		fscan(sFileIn, "%z", &sString); //put line in sString		  		
		iLengthOfString = sizerc(sString); //check how many characters are present in the string at hand
		
		if(sString[iLengthOfString-1] == 'Q'){ //check if stock is traded on the 'right' exchange
			iTellerQ = iTellerQ +1;	//counts number of observations on the 'right' exchange 
			if((sString[iLengthOfString-3] == 'E' ) || (sString[iLengthOfString-3] == 'F') || (sString[iLengthOfString-3] == ',')){ //check if there is a sale condition available, if there is it needs to be denoted either by an'E' or 'F'
				fReplaceCharInString(&sString, ':', ',' ); //replace ":" by "," for later purposes
				iTellerSaleCondition = iTellerSaleCondition + 1;	 //counts number of observations with 'right' sale condition
				fprint(sFileOut, sString, "\n"); //if conditions are met, cleaned observation is put in sFileOut
			}
		}
	}
	
	fclose(sFileIn); //close files, they are not needed anymore and can be saved
	fclose(sFileOut);

	//print some useful statistics of number of original, filtered and dropped observations
	println("                    Nobs       filtered         dropped");						
	println("FilterStep 2       ", iN, "      ", iTellerQ, "           ", iN-iTellerQ);
	println("                    Nobs       filtered         dropped");						
	println("FilterStep 3       ", iTellerQ, "      ", iTellerSaleCondition, "           ", iTellerQ-iTellerSaleCondition);						
}

/*
**  Function:	Create a vector with dates that Ox can handle
**
**  Input: 		vDateNumber: vector with date numbers that Ox cannot handle
**
**  Output: 	vDate: vector with dates that Ox can handle
*/

fCreateDateVector(const vDateNumber){

	decl iN, vYear, vMonth, vDay;
	iN = sizerc(vDateNumber); //count number of date numbers in vector
	vYear = floor(vDateNumber/10000); //extract year from date number							
	vMonth = floor((vDateNumber-floor(vDateNumber/10000)*10000)/100); //extract month from date number		
	vDay = vDateNumber-floor(vDateNumber/100)*100; //extract day from date number							
	
	decl vDate, i;
	vDate = zeros(iN,1); //initialize new date vector
	for(i = 0; i < iN; i++){
		vDate[i] = dayofcalendar(vYear[i], vMonth[i], vDay[i]); //fill new date vector		
	}
	
	return vDate;
}

/*
**  Function:	Check whether prices are all strictly positive
**
**  Input: 		vPrice: Vector with trade prices
**				vClean: Vector with ones with rows that are clean and zeros that need to be dropped
**
*/

fCleaningandFilteringStep4(const vPrice, const avClean){

  avClean[0]= avClean[0] + (vPrice .> 0);
}

/*
**  Function:	Check whether there are corrected prices
**
**  Input: 		vPrice: Vector with indicator whether trade is corrected
**				vClean: Vector with ones with rows that are clean and zeros that need to be dropped
**
*/

fCleaningandFilteringStep5(const vCorr, const avClean){

	avClean[0] = avClean[0] + (vCorr .== 0);

}

/*
**  Function:	Compress trades to a maximum of one per second
**
**  Input: 		mX: Matrix with data
**
**  Output: 	mXUni: Matrix with a maximum of one trade per second
**
*/

fCleaningandFilteringStep6(const mX){

	decl iN, iCount, iUniT; 
	iN 		= rows(mX);	// number of observations
	iCount	= 1; //count number of trades per second
	iUniT	= 0; // Indicator that counts the number of unique seconds

	decl vSameTime, mXUni;
	vSameTime = zeros(rows(mX),1);	
	vSameTime[1:] = mX[1:][10] .== mX[:iN-2][10]; // give value 1 in vector if a trade has the same time as the next trade
	mXUni = zeros(iN - sumc(vSameTime), columns(mX)); //initialize matrix with maximum of one trade per second 

	for(decl k = 1; k < iN; k++){
		if(vSameTime[k] == 1){
			iCount ++;
		}else{	// compress trading data that occured at the same second
		mXUni[iUniT][:3] = mX[k-1][:3];
		mXUni[iUniT][4] = quantilec(mX[k-iCount:k-1][4]);  // get median price
		mXUni[iUniT][5] = sumc(mX[k-iCount:k-1][5]); // sum volumes
		mXUni[iUniT][6:10] = mX[k-1][6:10];
		iUniT ++;
		iCount = 1;
		}
	}
	if(iCount == 1){ //Last trade is a unique second
		mXUni[iUniT][] = mX[iN-1][];
	}else{ //Last trade is a not a unique second
		mXUni[iUniT][:3] = mX[iN-1][:3];
		mXUni[iUniT][4] = quantilec(mX[iN-1-iCount:iN-1][4]); // get median price
		mXUni[iUniT][5] = sumc(mX[iN-1-iCount:iN-1][5]); // sum volumes
		mXUni[iUniT][6:10] = mX[iN-1][6:10]; 
	}

	//print some useful statistics of number of original, filtered and dropped observations
	print("                    Nobs       filtered         dropped");						
	print("\nFilterStep 6       ", iN, "      ", int(iN - sumc(vSameTime)), "           ", int(sumc(vSameTime)), "\n");

	return mXUni;		
}

/*
**  Function:	Filter outliers by using 10*MAD
**
**  Input: 		mX: Matrix with data
**
**  Output: 	mXNew: Same matrix, but rows are deleted with outliers
**
*/

fCleaningandFilteringStep7(const mX){

	decl iN, dNumberDeviations, iRowNumber, dMedian, dMAD, vDate, vPrices, mXNew, dLags;
	iN = rows(mX);
	dNumberDeviations = 10;
	dLags =	25; //number of trades to compute MAD
	iRowNumber= 0;
	dMedian	= 0;
	dMAD = 0;
	vDate = mX[][0]; //Date			
	vPrices = mX[][4];	//Prices

	decl vNewDay, iT;
	vNewDay = ones(iN,1); 
	vNewDay[1:] = vDate[1:] .!= vDate[:iN-2];  //check when a new day begins
	iT = rows(vNewDay);
	
	decl j;
	for(j = 0; j < iT; j++){
	    if(vNewDay[j] == 1){
			vNewDay[j] = j;	//determine row numbers where new day begins
		}
	}
	vNewDay = selectifr(vNewDay, vNewDay.!=zeros(rows(vNewDay),1));	//select prices on a new day

	decl iPastDays, dFirstTradeOfDay, dLastTradeOfDay;
	iPastDays = 0; //number of days previously 'visited'
	dFirstTradeOfDay = 0; //index number of first trade of the day
	dLastTradeOfDay = vNewDay[iPastDays]-1;	 //index number of last trade of the day

	decl vThrowAway, k;
	vThrowAway = zeros(iN,1);
	for(k = 0; k < iN; k++){
		if(k == dFirstTradeOfDay){
			dMedian = quantilec(vPrices[k+1:min(k+dLags,dLastTradeOfDay)]); //compute median
			dMAD = meanc(fabs(vPrices[k+1:min(k+dLags,dLastTradeOfDay)]-dMedian));	//determine MAD value 
		}else if(k != dLastTradeOfDay){ // k is not last trade of the day
			dMedian = quantilec(vPrices[max(k-dLags,dFirstTradeOfDay):k-1]|vPrices[k+1:min(k+dLags,dLastTradeOfDay)]); //calculate median of 25 before + 25 after prices
			dMAD = meanc(fabs(vPrices[max(k-dLags,dFirstTradeOfDay):k-1]-dMedian)|fabs(vPrices[k+1:min(k+dLags,dLastTradeOfDay)]-dMedian));	//determine MAD value 
		}
		else{ //k == dLastTradeOfDay
			dMedian = quantilec(vPrices[max(k-dLags,dFirstTradeOfDay):k-1]); 
			dMAD = meanc(fabs(vPrices[max(k-dLags,dFirstTradeOfDay):k-1]-dMedian));	//determine MAD value
			dFirstTradeOfDay = dLastTradeOfDay + 1;	//set first trade of the day equal to index number of last trade of the previous day
			iPastDays = iPastDays + 1;
			if (iPastDays < rows(vNewDay) - 1) {
				dLastTradeOfDay = vNewDay[iPastDays]-1; //index number of last trade of a day
			}else{
				dLastTradeOfDay = iN-1; //very last day
			}
		}
		vThrowAway[k] = (vPrices[k] > dNumberDeviations*dMAD+dMedian	|| vPrices[k] < -dNumberDeviations*dMAD+dMedian);	//if outlier
	}
	
	mXNew = selectifr(mX, !vThrowAway); //select trades that are not outliers

	return mXNew;
}

main(){

	fReadingStringData(); //start off by reading the 'dirty' files by reading string data en performance of cleaning step 2 and 3. 

	decl mData;
	mData = loadmat("2014_Cleaned.csv"); //Load data where cleaning steps 1, 2 and 3 are already processed on 
	
	decl vDate, vTimeValues, mX;
	vDate = fCreateDateVector(mData[][0]); //create vector with dates oX can work with;
	vTimeValues = mData[][1] / 24 + mData[][2] / 24 / 60 + mData[][3] / 24 / 60 / 60; //extract the time of each trade	
	mX = mData~vDate~vTimeValues; //put it all together								  	  	
	
	decl vClean1;
	fCleaningandFilteringStep4(mX[][4], &vClean1); //Only select trades strictly larger than 0	           
	mX = selectifr(mX[][], vClean1); //select trades after cleaning step 4					

	decl vClean2;
	fCleaningandFilteringStep5(mX[][7], &vClean2); //Only select trades that are not corrected	  	
	mX= selectifr(mX[][], vClean2);	//select trades after cleaning step 5				

	mX =  mX ~ (mX[][0] + mX[][9]); //create unique value for each second and concatenate with mX   	
	mX = fCleaningandFilteringStep6(mX); //Create median prices when multiple trades within a second		  	  
	mX = fCleaningandFilteringStep7(mX); //Delete outliers

	savemat("2014_Cleaned.csv", mX); //save final cleaned data file
}