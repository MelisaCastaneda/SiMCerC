 // Cervix Screening Program CSP
// Marcel Greuter
// 22 july 2021 - initial version 01
// 24 july 2021 - version 02 sensitivities adapted according to values in Model01.pptx
// 1 oct 2021 - version 03 MG: changed and corrected transition probabilities
// 14 nov 2021 - trial version 04 MG: first cost-effectiveness version
// 19 nov 2021 - bug fixed in readCumDeathProb, version ready for testing
// 29 nov 2021 - final version 04 for testing... :)
// 31 jan 2022 - after revision of the flow chart, version 05 was build, now ready for testing...
// 4 feb 2022 - Version with new flow chart implemented
// 5 feb 2022 - Debug function added to check flow chart
// 6 feb 2022 -	survival as a function of FIGO state coded
// to do: false positives have to be added, compliance function has to be adapted, all parameters to an input file
// 11 feb 2022 - version 06. adapted det file
// 23 mrt 2022 - version 07. New flow chart
// 31 mrt 2022 - update FIGO distribution and Survival (IKNL data)
// 10 apr 2022 - Update survival as a function of stage and age: see routine TumorSurvival
// 28 apr 2022 - Update Confirmatory diagnosis
// 9 may 2022 - Check on opening output files, corrected CIN2 and CIN3 to HPV+ transition
// 1 july 2022 - Add a screening test at age 65 if HPV+
// 8 july 2022 - Version 11. Added number of detected stages in the output
// 13 july 2022 - Number of detected stages per age in the output, compliance function also includes age 65
// 21 july 2022 - Ndetected only counted for direct referrals. There was no Ndetected for age 65 because there are 8 screen ages, not 7: this has been corrected.
// 25 july 2022 - Version 12. False positives included
// 2 jan 2023 - Version 13b. Base case scenario
// 12 jan 2023 - Version 13_1. Scenario 1. Age of death added to output, independent of participation in screening
// 10 jul 2023 - Number of iterations include by Nit

#include <iostream> 
#include <cstdlib> 
#include <time.h>
#include <cmath>
#include <conio.h>

#define ScreeningOn 1// turn screening on or off
#define Nwomen 100000
#define debug 0
#define Nit 10

using namespace std; 

int CurrentState = 0 ;
int age = 0 ;
float cumDeathProb[101] ;
int NFalsePositive = 0 ;
float costs ;
FILE *fp2 ;
int FIGOstate ;

int DetScreenAge = 0 ;
int DetCompliance = 0 ;
int DetReferral = 0 ;
int DetState = -1 ;
int DetFIGO = -1 ;
float DetDeathAge = -1 ;
float DetNaturalDeathAge = -1 ;
int DetTreatment = 0 ;												// 13-01-23 added

int Ndetected[8][6] ; 												// detected lesions for 8 ages: 30, 35..65 and 6 stages: HPV-, HPV+, ... cancer 

double Random()
{
	return (double)rand() / (double)((unsigned)RAND_MAX + 1 ) ;
}

void HPVmin( void )													// CurrentState = 0
{
	double rnd = Random() ;
	if ( ( age >= 15 ) and ( age <= 24 ) and ( rnd < 0.074 ) )
		CurrentState = 1 ;
	if ( ( age >= 25 ) and ( age <= 34 ) and ( rnd < 0.098 ) )
		CurrentState = 1 ;
	if ( ( age >= 35 ) and ( age <= 44 ) and ( rnd < 0.050 ) )
		CurrentState = 1 ;
	if ( ( age >= 45 ) and ( age <= 54 ) and ( rnd < 0.036 ) )
		CurrentState = 1 ;
	if ( ( age >= 55 ) and ( rnd < 0.027 ) )
		CurrentState = 1 ;

}

void HPVplus( void )												// CurrentState = 1
{
	double rnd = Random() ;
	if ( rnd < 0.5034 ) 
		CurrentState = 0 ;
	else
	{
		if ( rnd < 0.5034+0.0610 ) 
			CurrentState = 2 ;
		else
		{
			if ( rnd < 0.5034+0.0610+0.0034 ) 
				CurrentState = 3 ;
			else
			{
				if ( rnd < 0.5034+0.0610+0.0034+0.0072  ) 
						CurrentState = 4 ;
			}
		}
	}
}

void CIN1( void )													// CurrentState = 2
{	
	double rnd = Random() ;
	if ( rnd < 0.4000 ) 
		CurrentState = 0 ;
	else
	{
		if ( rnd < 0.4000+0.1600 ) 
			CurrentState = 1 ;
		else
		{
			if( rnd < 0.4000+0.1600+0.0242 ) 
				CurrentState = 3 ;
			else
			{
				if( rnd < 0.4000+0.1600+0.0242+0.0047 ) 
					CurrentState = 4 ;
			}
		}
	}
}		

void CIN2( void )													// CurrentState = 3
{	
	double rnd = Random() ;
	if ( rnd < 0.0073 ) 
		CurrentState = 1 ; 
	else
	{
		if ( rnd < 0.0073+0.0244 ) 
			CurrentState = 2 ;
		else
		{
			if ( rnd < 0.0073+0.0244+0.0474 ) 
				CurrentState = 4 ;
		}
	}
}

void CIN3( void )													// CurrentState = 4
{	
	double rnd = Random() ;
	if ( rnd < 0.0025 ) 
		CurrentState = 1 ;
	else
	{
		if ( rnd < 0.0025+0.0074 ) 
			CurrentState = 2 ;
		else
		{
			if ( rnd < 0.0025+0.0074+0.0486 ) 
				CurrentState = 3 ;
			else
			{
				if ( rnd < 0.0025+0.0074+0.0486+0.04) 
					CurrentState = 5 ;						
			}
		}
	} 
	
	
// testing:	CurrentState = 5 ;
}

double TumorSurvival( int age, int cc ) // Survivalprobability is given by y=a*exp(bt)+c
{
	double y = Random() ;
	double t ;
	int j ;
	const double a[5][7] = { {  2.010455946,  1.587227382,   1.891873867,  4.999864589,  9.895077586,  50.0,         76.9734858 },
	                         { 13.85142061,	 12.08714041,   19.16786749,  37.30197723,  35.00436499,  107.2434008,   95.65658272 },	
							 { 26.45198395,	 33.74632466,   31.50795519,  51.74147368,  62.66441825,   62.33323225,  97.42502489 },	
							 { 53.52169224,	 47.22221442,   61.76505583,  57.07367603,  71.96671389,   80.67383325,  96.39160654 },
							 { 75.4008106,	 57.7777899,    80.00053313,  86.22845415,  84.5168083,    96.91882121,  98.86993521 } } ;
	const double b[5][7] = { { -0.114209428, -0.206465909,  -0.168599918, -0.800033383, -0.167779519,  -0.3,         -0.303121298 },
							 { -0.34222783,	 -0.204104291,  -0.131019886, -0.09323743,  -0.162942508,  -0.083446662, -0.340631829 },
							 { -0.450706196, -0.195386638,  -0.221972525, -0.165972035, -0.109460283,  -0.239899205, -0.29158698 },
							 { -0.619815651, -0.36628054,   -0.345993317, -0.724967985, -0.388531728,  -0.350882002, -0.725008406 },
							 { -1.107243477, -1.044161075,  -0.667587916, -0.792290275, -1.073242708,  -1.079547586, -1.776686193 } } ;
	const double c[5][7] = { { 97.98954405,  98.41277262,  98.10812614,   94.71832646,  90.10492241,   50.0,         23.02651419 },	
							 { 86.14857938,  87.91285959,  80.83213248,   62.6980227,   64.99563502,   -7.243400753,  4.343417338 },
							 { 73.54801605,  66.25367535,  68.49204482,   48.25852638,  37.33558174,   37.66676776,   2.574975107 },
							 { 46.47820389,  52.77768527,  38.23494425,   42.92632398,  28.03328611,   19.32616675,   3.608393459 },
							 { 24.5991894,   42.22210992,  19.99946687,   13.77154585,  15.4832913,     3.081277527,  1.130064789 } } ;
	if ( ( age >= 20 ) && ( age < 30 ) ) { j = 0 ; } ;
	if ( ( age >= 30 ) && ( age < 40 ) ) { j = 1 ; } ;
	if ( ( age >= 40 ) && ( age < 50 ) ) { j = 2 ; } ;
	if ( ( age >= 50 ) && ( age < 60 ) ) { j = 3 ; } ;
	if ( ( age >= 60 ) && ( age < 70 ) ) { j = 4 ; } ;
	if ( ( age >= 70 ) && ( age < 80 ) ) { j = 5 ; } ;
	if   ( age >= 80 )                   { j = 6 ; } ;
	y = y * 100 ;
	if ( y > c[cc][j] )
		t = log((y-c[cc][j])/a[cc][j])/b[cc][j] ;
	else
		t = 100 ;	
	if (debug) { printf( "Survival = %d %d %d %f %f", age, cc, j, y, t ); getch() ; }
	return t ;
}

void readCumDeathProb() 
{
	FILE *fp;
	if ((fp = fopen( "CumDeathProb.txt", "r"))==NULL) {
		std::cerr << std::endl << "Error: Cannot open input file CumDeathProb.txt." << std::endl;
		exit(1);
	}
	for (int i = 0; i < 101; i++) {									// Range should be to 101, was 100, bug :)
		fscanf(fp, "%f\n", &cumDeathProb[i]);
	}
	fclose( fp ) ;													// File was not closed in original version
	return;
}

float findNaturalDeathAge(float rndNaturalDeath) 					//Find the natural death age for this women
{
	float result;
	if ( rndNaturalDeath <= cumDeathProb[0] ) {
		return 0;
	}
	if ( rndNaturalDeath >= cumDeathProb[100] ) {
		return 100;
	}
	else {
		int i = 0;
		while ( !( ( cumDeathProb[i] <= rndNaturalDeath ) && ( cumDeathProb[i+1] >= rndNaturalDeath ) ) ) 
			i++;
		result = ( rndNaturalDeath - cumDeathProb[i] ) / ( cumDeathProb[i+1] - cumDeathProb[i] ) + i; //Linear interpolation
		if (result > 100) {
			printf("\nError: natural death age larger than 100: %f for probability %f and i %d", result, rndNaturalDeath, i);
			exit(1);
		}
		return result;
	}
}

int HPVTestGP( void )												// Perform HPV test bij GP
{
	double rnd1 = Random() ;
	double rnd2 = Random() ;
	int result = 0 ;
	 
	result = (  ( ( CurrentState >= 3 ) && ( rnd1 < 0.964 ) ) ); // Sensitivity  HPV GP
	if ( !result )
	{
		if ( rnd2 > 0.942 )                                        // false positive detected
		{
			result = 1 ;
			NFalsePositive++ ;
		}
	} 
	return result ;					
}

int HPVTestSS( void )												// Perform HPV test bij Self-sampling
{
	double rnd1 = Random() ;
	double rnd2 = Random() ;
	int result = 0 ; 
	
	result = (  ( ( CurrentState >= 3 ) && ( rnd1 < 0.929 ) ) ); // Sensitivity HPV SS	
 	if ( !result )
	{
		if ( rnd2 > 0.939 )                                        // false positive detected
		{
			result = 1 ;
			NFalsePositive++ ;
		}
	} 			
	return result ;	
}

int MethylationTest( int GP )										// Perform Methylation test
{
	double rnd1 = Random() ;
	double rnd2 = Random() ; 
	int result = 0 ;
	if ( GP )
		costs += 0 ;                                       	// costs of Methylation test of GP: TO BE CHANGED
	else
		costs += 0 ;											// costs of Methylation test of SS: TO BE CHANGED
	result = (  ( ( CurrentState >= 3 ) && ( rnd1 < 0.80 ) ) ) ;  	// sensitivity for Methylation test
	if ( !result )
	{
		if ( rnd2 > 0.85 )											// false positive detected
		{
			result = 1;
			NFalsePositive++ ;
		}
	} 
	return result ;	
}


int ControlCytologyTest( void )										// Perform Control Cytology
{
	double rnd1 = Random() ; 
	double rnd2 = Random() ;
	int result = 0 ;
	costs += 55.14 ;                                        		// costs of Control Cytology test
	result = (  ( ( CurrentState >= 3 ) && ( rnd1 < 0.864 ) ) ) ;  	// sensitivities for Control Cytology test
	if ( !result )
	{
		if ( rnd2 > 0.874 )											// false positive detected
		{
			result = 1;
			NFalsePositive++ ;
		}
	} 
	return result ;	
}

int ScreenAge( void )
{
	return ((age==30)||(age==35)||(age==40)||(age==45)||(age==50)||(age==55)||(age==60)||(age==65)) ;
}

int Compliance( void )												// Is the woman compliant in the screening?
{
	double rnd = Random() ; 
	return (( ( age == 30 ) && ( rnd < 0.441 ) ) ||					// NB: this has to be adapted.Melisa: I think we must used the monitor participation.
			( ( age == 35 ) && ( rnd < 0.493 ) ) ||
			( ( age == 40 ) && ( rnd < 0.560 ) ) ||
			( ( age == 45 ) && ( rnd < 0.587 ) ) ||
			( ( age == 50 ) && ( rnd < 0.604 ) ) ||			
			( ( age == 55 ) && ( rnd < 0.609 ) ) ||
			( ( age == 60 ) && ( rnd < 0.613 ) ) ||
			( ( age == 65 ) && ( rnd < 0.613 ) ) ) ;
}

int GPTest( void )													// GP performs HPV and Cytology test
{
	int test = 0 ;
	costs += 60.34 ;												// Cost of GP test
	test = HPVTestGP() ;
	return test ;
}

int SelfTest( void )												// Selftest for HPV
{
	int test = 0 ;
	costs += 105.07 ;												// Cost of self-test
	test = HPVTestSS() ;
	return test ;
}

int BiopsyTest( void )
{
	double rnd = Random() ;
	int result = 0 ;
	costs += 63.77 ;												// costs of biopsy
	result = ( CurrentState >= 3 ) ;
	return result ;									
}

int DetermineFIGOstate( void )
{
	int cc ;
	double rnd = Random() ;
	if ( age < 25 )
  	{
		if ( rnd < 0 )
  			cc = 0 ;
  		else
  		{
  			if ( rnd < 0+ 0.40 )
  				cc = 1 ;
  			else
  			{
  				if ( rnd < 0 + 0.40 + 0.60 )
  					cc = 2 ;
  				else
				{
					if( rnd < 0 + 0.40 + 0.60+0 )
						cc = 3;
					else
						cc = 4 ;
				}
  			}
		}
	}
	if ( ( age >= 25 ) && ( age < 35 ) )
  	{
		if ( rnd < 0.32 )
  			cc = 0 ;
  		else
  		{
  			if ( rnd < 0.32 + 0.48 )
  				cc = 1 ;
  			else
  			{
  				if ( rnd < 0.32 + 0.48 + 0.15 )
  					cc = 2 ;
				else
				{
					if (rnd < 0.32 + 0.48 + 0.15 +0.02)
						cc = 3;
					else
				  		cc = 4 ;
				}
  			}
		}
	}
	if ( ( age >= 35 ) && ( age < 50 ) )
  	{
		if ( rnd < 0.27)
  			cc = 0 ;
  		else
  		{
  			if ( rnd < 0.27 + 0.43 )
  				cc = 1 ;
  			else
  			{
  				if ( rnd < 0.27 + 0.43 + 0.20 )
  					cc = 2 ;
				else
				{
					if (rnd < 0.27 + 0.43 + 0.20 +0.04)
						cc = 3;
					else
						cc = 4 ;	
				}
  			}
		}
	}
	if ( ( age >= 50 ) && ( age < 65 ) )
  	{
		if ( rnd < 0.13)
  			cc = 0 ;
  		else
  		{
  			if ( rnd < 0.13 + 0.33 )
  				cc = 1 ;
  			else
  			{
  				if ( rnd < 0.13 + 0.33 + 0.30 )
  					cc = 2 ;
  				else
				{
					if (rnd < 0.13 + 0.33 + 0.30 +0.08)
						cc = 3;
					else
						cc = 4 ;	
				}
  			}
		}
	}
	if ( ( age >= 65 ) && ( age < 75 ) )
  	{
		if ( rnd < 0.02 )
  			cc = 0 ;
  		else
  		{
  			if ( rnd < 0.02 + 0.25 )
  				cc = 1 ;
  			else
  			{
  				if ( rnd < 0.02 + 0.25 + 0.39 )
  					cc = 2 ;
  				else
				{
					if (rnd < 0.02 + 0.25 + 0.39 +0.11)
						cc = 3;
					else
						cc = 4 ;	
				}	
  			}
		}
	}
	if ( age >= 75 ) 
  	{
		if ( rnd < 0.03 )
  			cc = 0 ;
  		else
  		{
  			if ( rnd < 0.03 + 0.17 )
  				cc = 1 ;
  			else
  			{
  				if ( rnd < 0.03 + 0.17 + 0.42 )
  					cc = 2 ;
  				else
				{
					if (rnd < 0.03 + 0.17 + 0.42 +0.18)
						cc = 3;
					else
						cc = 4 ;	
				}		
  			}
		}
	}
	return cc ;
}

int ConfirmatoryDiagnosis( void )									// Do Confirmatory Diagnosis 
{
	int result = 0 ;
	if (debug) printf("CurrentSate=%d", CurrentState);
	if ( DetReferral == 1 )											// 21-07-22: only count the number of detected lesions at direct referral
		Ndetected[(age-30)/5][CurrentState]++;
	switch( CurrentState ) 
	{
		case 0: costs +=  344.85 ; break ;
		case 1: costs +=  344.85 ; break ;
		case 2: costs += 1073.60 ; break ;
  		case 3: costs += 1588.92 ; break ;
  		case 4: costs += 1860.90 ; break ;
  		case 5: FIGOstate = DetermineFIGOstate() ;
  				switch ( FIGOstate )
  				{
  					case 0: costs +=  6094.13 ; break ;
  					case 1: costs += 14451.27 ; break ;
  					case 2: costs += 14244.36 ; break ;
  					case 3: costs += 14244.36 ; break ;
					case 4: costs += 14244.36 ; break ;
			    }
			    if (debug) printf("FIGOstate=%d",FIGOstate);
  				break ;
	} 
	result = ( CurrentState >= 3 ) ;                            	// Confirmatory Diagnosis is positive when >= CIN2
	return result ;
}

void InitDetParam( void )
{
	DetScreenAge = 0 ;
	DetCompliance = 0 ;
	DetReferral = 0 ;
	DetState = -1 ;
	DetFIGO = -1 ;
	DetDeathAge = -1 ;
	DetNaturalDeathAge = -1 ;
	DetTreatment = 0 ;											// 13-01-23 added
}

int DoScreening( void )											// Perform screening
{
	int HPVtest = 0, CytologyTest = 0, stop = 0 ;
	int TumorFound = 0 ;	
	if ( ScreenAge() )
	{
		DetScreenAge = 1 ;
		if ( Compliance() ) 											// check whether woman complies to screening
		{
			DetCompliance = 1 ;
	 		if ( debug ) printf("1") ;
			if ( Random() < 0.91 )								// 91% goes to GP, 9% does self-test
			{
				if ( debug ) printf("2");
				HPVtest = GPTest() ;
				if ( HPVtest )
				{
					if (debug)printf("4");
					CytologyTest = MethylationTest( 1 ) ;
				}
				else 
				{
					if (debug)printf("5");
					stop = 1 ;
				}
			}
			else
			{
				if (debug) printf("3");
				HPVtest = SelfTest() ;
				if ( HPVtest ) 
				{
					if (debug)printf("6");
					CytologyTest = MethylationTest( 0 ) ;
				}
				else
				{
					if (debug)printf("7");
					stop = 1 ;
				}
			}
			if ( !stop )
			{
				stop = 1 ;
				if ( CytologyTest )  
				{                                 						// direct referral
					if (debug) printf("a");
					stop = Random() > 0.70 ;							// participate in confirmatory diagnosis?
					if ((debug) && (stop)) printf("c");
					if ((debug) && (!stop))printf("d");
					if (!stop)
						DetReferral = 1 ;
				}
				else													// indirect referral
				{
					if (debug) printf("b");
				}
			}
			if ( !stop )
			{
				if ( ConfirmatoryDiagnosis() )
				{
					if (debug) printf("l");
					if ( CurrentState == 5 )						// FIGO present and found
					{
						if (debug) printf("m");
						TumorFound = 1 ;
					}
					else
					{
						if (debug) printf("n");
						DetTreatment = 1 ;							// 13-01-23 added
//						CurrentState = 0 ;  						// Treatment of CIN2 and CIN3 // 13-01-23 commented out
					}
				}
				else
				{
					if (debug) printf("k");
				}
			}
		}
		else
		{
			if ( debug )
				printf( "0" ) ;
		}
	}
	return TumorFound ;
}


void Loop( int it )
{
	int woman ;
	int TumorFound ;
	float NaturalDeathAge, TumorDeathAge, DeathAge ;
	FILE *fp1 ;
	int death = 0 ;
	float TotalCost = 0 ;
	int TotalFalsePositive = 0 ;
	int i, j ;
	char buf[2], fn[80] ;
	
	if (!debug) printf( "Running " ) ;
	
	itoa( it, buf, 10 ) ;
  	strcpy( fn, "csp14_1" ) ;
  	strcat( fn, buf ) ;
  	strcat( fn, ".out" ) ;
	if ((fp1 = fopen( fn, "w"))==NULL) {
		printf( "Error: Cannot open output file ", fn ) ;
		exit(1);
	}
	itoa( it, buf, 10 ) ;
  	strcpy( fn, "csp14_1" ) ;
  	strcat( fn, buf ) ;
  	strcat( fn, ".det" ) ;
	if ((fp2 = fopen( fn, "w"))==NULL) {
		printf( "Error: Cannot open output file ", fn ) ;
		exit(1);
	}
	for ( i = 0 ; i < 8 ; i++ )
	for ( j = 0 ; j < 6 ; j++ )
		Ndetected[i][j] = 0 ;
	for ( woman = 0 ; woman < Nwomen ; woman++ )
	{
		if ((woman % 1000 == 0 )&& (!debug))
			printf( "." ) ;
		if (debug) printf( "%d", woman) ;
		fprintf( fp1, "%d\t", woman ) ;
		fprintf( fp2, "%d\t", woman ) ;
		age = 0 ;
		TumorFound = 0 ;
		CurrentState = 0 ;
		NaturalDeathAge = findNaturalDeathAge( Random() ) ;
		DeathAge = NaturalDeathAge ;
		fprintf( fp2, "%.2f\t", DeathAge ) ; 						// 12 jan 2023: Natural death age added for all women
		death = 0 ;
		costs = 0 ;
		NFalsePositive = 0 ;
//		while ( ( age < 100 ) && ( !TumorFound ) && ( !death ) )
		while ( ( age < 100 ) && ( !death ) )
		{
			switch ( CurrentState )
			{
				case 0: HPVmin() ; break ;
				case 1: HPVplus() ; break ;
				case 2: CIN1() ; break ;
				case 3: CIN2() ; break ;
				case 4: CIN3() ; break ;
			}
			fprintf( fp1, "%d\t", CurrentState ) ;
			if ( ScreeningOn )
			{
				InitDetParam() ;
				TumorFound = DoScreening() ;
				if ((TumorFound ) && ( CurrentState == 5 ) )				// Calculate survival FIGO 
				{
					TumorDeathAge = float(age) + TumorSurvival(age, FIGOstate) ;
					DeathAge = std::min( TumorDeathAge, NaturalDeathAge ) ; 
					death = 1 ;
//					fprintf( fp1, "%d ", CurrentState ) ;
					DetFIGO = FIGOstate ;
					DetDeathAge = TumorDeathAge ;
					DetNaturalDeathAge = NaturalDeathAge ;
//					fprintf( fp2, "%d\t%d\t%d\t%d\t%d\t%f\t%f", age, DetCompliance, DetReferral, DetState, DetFIGO, DetDeathAge, DetNaturalDeathAge ) ;
//					fprintf( fp2, "%d %d %d %d %.2f %.2f %d\n", woman, age, CurrentState, FIGOState, DeathAge, costs, NFalsePositive ) ;
					if (debug){printf ( "FOUND"); getch() ;}
				}	
				if ( DetScreenAge )
				{													// 13-01-23 added
					fprintf( fp2, "%d\t%d\t%d\t%d\t%d\t%.2f\t%.2f\t\%.2f\t", age, DetCompliance, DetReferral, CurrentState, DetFIGO, DetDeathAge, NaturalDeathAge, costs ) ;
					if ( DetTreatment )
					{												// 13-01-23 added
						CurrentState = 0 ;  						// Treatment of CIN2 and CIN3) // 13-01-23 added
						DetTreatment = 0 ;							// 13-01-23 added	
					}												// 13-01-23 added
				}													// 13-01-23 added
			}
			if (!death)
			{
				if ( ( CurrentState < 5 ) && ( age+1 > NaturalDeathAge ) )
				{
					death = 1 ;
					CurrentState = 6 ; 								// woman dies natural death
					fprintf(fp1,"%d\t", CurrentState ) ;
//					fprintf( fp2, "%d %d %d * %.2f %.2f %d\n", woman, age, CurrentState, NaturalDeathAge, costs, NFalsePositive ) ;
				}	
				else 
				{
					if ( CurrentState == 5 )					// Calculate survival FIGO 
					{
						TumorDeathAge = float(age) + TumorSurvival( age, FIGOstate) ;
						DeathAge = std::min( TumorDeathAge, NaturalDeathAge ) ; 
						death = 1 ;
						fprintf( fp1, "%d\t", CurrentState ) ;
//						fprintf( fp2, "%d %d %d %d %.2f %.2f %d\n", woman, age, CurrentState, FIGOState, DeathAge, costs, NFalsePositive ) ;
						if (debug){printf ( "FOUND"); getch() ;}
					}
					else
						age++ ;
				}
			}
//			fprintf( fp1, "%d ", CurrentState ) ;
		}
		fprintf( fp1, "\n" ) ;
		fprintf( fp2, "\n" ) ;
		if ( age >= 100 )
		{
			death = 1 ;
			CurrentState = 6 ; 										// woman dies natural death
//			fprintf( fp2, "%d %d %d * %.2f %.2f %d\n", woman, age, CurrentState, NaturalDeathAge, costs, NFalsePositive ) ;
		}	
		TotalCost += costs ;
		TotalFalsePositive += NFalsePositive ;
		if (debug) printf("\n") ;
	}
	fclose( fp1 ) ;
	fclose( fp2 ) ;
	printf( "\n" ) ;
	printf( "Total costs: %.2f\n", TotalCost ) ;
	printf( "Total FP: %d\n", TotalFalsePositive ) ;
	for ( i = 0 ; i < 8 ; i++ )
	{
		printf( "detected age %d:\t", i*5+30 ) ;	
		for ( j = 0 ; j < 6 ; j++)
			printf( "%d\t", Ndetected[i][j] ) ;
		printf( "\n" ) ;
	}
	
	printf( "Total FalsePositives: %d\n", TotalFalsePositive ) ;
}

int main()
{
	int it ;
	
	printf( "Cervix Screening Program CSP - scenario_I\n" ) ;
	printf( "Test version 14_1 - 10 jul 2023\n" ) ;
	printf( "Number of women: %d\n", Nwomen ) ;
	printf( "For every iteration i:\n" ) ;
	printf( "csp14_1i.out: for every woman: state per age\n" ) ;
	printf( "csp14_1i.det: for every woman: screening age, participate, referral, state, FIGO, tumor death, death age, costs\n" ) ;

	srand( time(NULL) ) ;												// Initialize random number generator
	readCumDeathProb() ;

	for ( it = 0 ; it < Nit ; it++ )
		Loop( it ) ;
}
