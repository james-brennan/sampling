
/*


	c++ implementation of sampling mcmc

	This computes

*/
#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()
#include <iostream>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>    // std::min_element, std::max_element
#include <gsl/gsl_statistics.h> // for percentiles calculation

using namespace std;

struct cloudiness_month 
{
    int month;
    int nXSize;
    int nYSize;
	vector<float> cloudField;
};

struct maskRaster 
{
	int nXSize;
	int nYSize;
	vector<float> mask;
};


struct sampleRaster 
{
	int nXSize;
	int nYSize;
	vector<int> raster;
};

/*
IO stuff
*/
vector<cloudiness_month> prepareData() {
	/*
	This method loads the data
	*/
	GDALDataset  *poDataset;
    GDALAllRegister();
	std::string filename;
    // GDALDataset*[] = datasets;

	vector<string> months = {"01", "02","03","04","05","06","07","08","09","10","11","12"};

	// want to store all date
	vector<cloudiness_month> cloudiness_field;

    for (int m=0; m<12; m++) 
    {
		long long month = m;
		filename = "0.25_MODCF_monthlymean_" + months[m] + ".tif";
		cout << "Loading " << filename << endl;
		const char * fc = filename.c_str();
    	poDataset = (GDALDataset *) GDALOpen( fc, GA_ReadOnly );
    	//cout << poDataset << endl;
    	/*
		
		Start data reading stuff

    	*/
    	GDALRasterBand  *poBand;
		int             nBlockXSize, nBlockYSize;
		int             bGotMin, bGotMax;
		double          adfMinMax[2];
		poBand = poDataset->GetRasterBand( 1 );
		GUInt16 *data;
		int   nXSize = poBand->GetXSize();
		int   nYSize = poBand->GetYSize();

		// data = (float *) CPLMalloc(sizeof(float)*nXSize*nYSize);

		data = (GUInt16 *) CPLMalloc(sizeof(GUInt16)*nXSize*nYSize);
		poBand->RasterIO( GF_Read, 0, 0, nXSize, nYSize, 
		                  data, nXSize, nYSize, GDT_UInt16, 
		                  0, 0 );
		// 
		cloudiness_month current_month;
		current_month = cloudiness_month();
		// 
		current_month.month = month;
		current_month.nXSize  = nXSize; 
		current_month.nYSize  = nYSize; 		
		for (int i=0; i<nXSize*nYSize; i++)
		{
			current_month.cloudField.push_back(data[i]*0.001*0.1);

		}
		cloudiness_field.push_back(current_month);
		// free memory
		CPLFree(data);
		GDALClose(poDataset);
    }
    return cloudiness_field;
}

maskRaster loadMask() 
{
	GDALDataset  *poDataset;
    GDALAllRegister();
	std::string filename;
	filename = "0.25_LW_mask.tif";
	const char * fc = filename.c_str();
	poDataset = (GDALDataset *) GDALOpen( fc, GA_ReadOnly );
	//cout << poDataset << endl;
	/*
	
	Start data reading stuff

	*/
	GDALRasterBand  *poBand;
	int             nBlockXSize, nBlockYSize;
	int             bGotMin, bGotMax;
	double          adfMinMax[2];
	poBand = poDataset->GetRasterBand( 1 );
	GUInt16 *data;
	int   nXSize = poBand->GetXSize();
	int   nYSize = poBand->GetYSize();

	// data = (float *) CPLMalloc(sizeof(float)*nXSize*nYSize);

	data = (GUInt16 *) CPLMalloc(sizeof(GUInt16)*nXSize*nYSize);
	//cout << data << endl;
	poBand->RasterIO( GF_Read, 0, 0, nXSize, nYSize, 
	                  data, nXSize, nYSize, GDT_UInt16, 
	                  0, 0 );
	// 
	maskRaster LW_MASK;
	LW_MASK = maskRaster();
	LW_MASK.nXSize = nXSize;
	LW_MASK.nYSize = nYSize;
	// 
	int value;
	const int BAD = 65535;
	for (int i=0; i<nXSize*nYSize; i++)
	{
		value = data[i];
		if (value == BAD)
		 {
			value = 0.0;
		}
		else {
			value = 1.0;
		}
		LW_MASK.mask.push_back(value);

	}
	return LW_MASK;
}

sampleRaster loadNumSamples(std::string filename) 
{
	/*
	Load the dataset containing the number 
	of samples expected for a sensor
	*/
	GDALDataset  *poDataset;
    GDALAllRegister();
	const char * fc = filename.c_str();
	poDataset = (GDALDataset *) GDALOpen( fc, GA_ReadOnly );
	//cout << poDataset << endl;
	/*
	
	Start data reading stuff

	*/
	GDALRasterBand  *poBand;
	int             nBlockXSize, nBlockYSize;
	int             bGotMin, bGotMax;
	double          adfMinMax[2];
	poBand = poDataset->GetRasterBand( 1 );
	GInt32 *data;
	int   nXSize = poBand->GetXSize();
	int   nYSize = poBand->GetYSize();

	// data = (float *) CPLMalloc(sizeof(float)*nXSize*nYSize);

	data = (GInt32 *) CPLMalloc(sizeof(GInt32)*nXSize*nYSize);
	//cout << data << endl;
	poBand->RasterIO( GF_Read, 0, 0, nXSize, nYSize, 
	                  data, nXSize, nYSize, GDT_Int32, 
	                  0, 0 );
	// 
	sampleRaster nSamples;
	nSamples = sampleRaster();
	nSamples.nXSize = nXSize;
	nSamples.nYSize = nYSize;
	// 
	int value;
	// const int BAD = -999;
	for (int i=0; i<nXSize*nYSize; i++)
	{
		value = data[i];
		nSamples.raster.push_back(value);
		//cout << value << " +" << endl;
	}
	return nSamples;
}


vector<double> doSamplingModel(vector<float> cldsXY,
					 int MaskValue, double daysBetweenSamples)
{
	/*
	The key part
	Does sampling on location x y

	Depending on the value of daysBetweenSamples
	we have to extend from sampling once per day up to N times per day

	*/
	const int N = 2; // number of iterations

	
	/*
	Now need to be clever and figure out how we do the sampling
	based upon the samplingRate. If the sampling rate is greater than 
	once per day we have to sample at more than 366 times...
	Approach for now is to use a set of samplingRates which correspond to
	1,2,3,4,5 times per day
	*/
	// dSampleRates[5] = {1.0, 0.5, 0.333, 0.25, 0.2}; // 1 .. 2 .. 3 .. 4 .. 5 samples per day
	/*
	Alternative approach is to do this manually
	*/
	int wegetSamples = 1;
	if (daysBetweenSamples < 1.0)
	{	
		double weneedSamples = 365.0 / daysBetweenSamples;
		wegetSamples = div((int)weneedSamples, 365).quot;

	}
	else {
		wegetSamples = 1;
	}
	vector<double> daysSample;
	const int mdays[12] = {31*wegetSamples,
						   28*wegetSamples,
						   31*wegetSamples,
						   30*wegetSamples,
						   31*wegetSamples,
						   30*wegetSamples,
						   31*wegetSamples,
						   31*wegetSamples,
						   30*wegetSamples,
						   31*wegetSamples,
						   30*wegetSamples,
						   31*wegetSamples};
	vector<int> clearDays;
	// so start as basis as 365 add if more per day and duplicate for 3 years
	clearDays.resize(365*wegetSamples*3);
	vector<int> observations;
	observations.resize(365*wegetSamples*3);


	/*
	now decide when we are observing the surface
	-- This is a regular thing, eg every day or every 5 days etc
	- to add necessary randomness vary the first sampletime
	*/

	int firstSample = 0; // for now let's not -- not sure how to make this work 
	for (int s=0; s<365*wegetSamples*3; s++) {
		/*
		Essentially we can get the nth of every opportunity 
		by getting the times when the quotient is zero of the division
		wegetSamples is the rate at which we get an observation
		*/
		int r = div(s, wegetSamples).rem;
		//cout << "remainder is:" << r << endl;
		if (r == 0.0) {
			// no remainder therefore a sample
			observations[s] = 1;
		}
		//cout << observations[s] << " -=" << endl;
	}

	/*
	trying to optmise random number generator
	--- that worked like crazy!
	*/
	random_device rd{}; // use to seed the rng 
	mt19937 rng{rd()}; // rng

	/*
	Make the bernoulli distributions
	*/
	vector<bernoulli_distribution> bernies;
	for (int m=0; m<12; m++) {
		bernoulli_distribution aDist(1.0-cldsXY[m]);
		bernies.push_back(aDist); 
	}

	for (int n=0; n<N; n++) {
	    int moffset = 0;
	    // now want to fill sample from the bernies
	    for (int m=0; m<12; m++)
	    {
	    	//cout << m <<"=======" << endl;
	    	// make a bernie
	    	// bernoulli_distribution mDist;
	    	bernoulli_distribution mDist = bernies[m];
	        //bernoulli_distribution mDist(1.0-cldsXY[m]); 
	    	// take n draws according to length of month
	    	for (int day=0; day<mdays[m];day++) 
	    	{
	    		moffset++;
	    		//cout << day << endl;
	    		if (mDist(rng))
	    		{
	    			// rememebr we want to duplicate before the year and after
	    			clearDays[moffset]=1;
	    			clearDays[(365*wegetSamples)-1+moffset]=1;
	    		    clearDays[(365*wegetSamples*2)-1+moffset]=1;
	    		}
	    		else
	    		{
	    			clearDays[moffset]=0;
	    			clearDays[(365*wegetSamples)-1+moffset]=0;
	    		    clearDays[(365*wegetSamples*2)-1+moffset]=0; 
	    		}
	    	}
	    }
	    //cout << "made sample" << endl;
	    // print out the sample
	    // vary DOB
	    int dob = 55*wegetSamples;
	    /*
		We now want to find the first observation before/after the 
		DOB that is not cloudy

	    */
	    // count days before and after to an observation
	    int before=0;
	    int after=1;
	    int bidx;
	    int aidx;
	    while (before < ((365*wegetSamples+dob)))
	    {
	    	//cout << sample[364+dob-before] << endl;
	    	// count days back before until a clear day
	    	bidx = 365*wegetSamples+dob-before-1;
	    	if ((clearDays[bidx]==0) ||
	    			(observations[bidx]==0)) {
	    		before++;
	    		//cout << before << endl;
	    	}
	    	else {
	    		// they are both 1
	    		break;
	    	}
	    }
	    while (after < (dob+(365*wegetSamples))) 
	    {
	    	// count days back before until a clear day
	    	aidx = 365*wegetSamples+dob+after-1;
	    	if ((clearDays[aidx]==0) ||
	    		(observations[aidx]==0)) {
	    		after++;
	    	}
	    	else {
	    		// they are both 1
	    		break;
	    	}
	    }
	    // return days between before and after
	    double days;
	    //cout << "samplesperday: " << wegetSamples << endl;
	    //cout << "before and after: " << before << " " << after << endl;
	    days = (double) after/wegetSamples + (double) before/wegetSamples;
	    // remember to convert between samplingRate and actual days
	    //cout << days << endl;
	    daysSample.push_back(days);
	}
    return daysSample;
}


vector <float> getCloudiness(vector<cloudiness_month> clouds, int loc)
 {
 	// extract for location
 	vector <float> cloudiness; 
 	for (int m=0; m<12; m++) 
 	{
 		// figure out location in flat array
 		int xoffset = clouds[m].nXSize; 
 		//cout << "xoffset is: "<< xoffset << endl;
 		//cout << clouds[m].cloudField.size() << endl;
 		cloudiness.push_back(clouds[m].cloudField[loc]);
 	} 
 	return cloudiness;
}

int getMaskValue(maskRaster msk, int loc)
 {
 	// extract for location
 	int value;
 	int xoffset = msk.nXSize; 
 	value = msk.mask[loc];
 	return value;
}

int getSampleCounts(sampleRaster nSamples, int loc)
 {
 	// extract for location
 	int value;
 	int xoffset = nSamples.nXSize; 
 	value = nSamples.raster[loc];
 	return value;
}

vector<double> runModel( vector<cloudiness_month> clds, maskRaster msk, 
					  sampleRaster numSamples) 
{
	/*
	This sets up the output and runs the model on the spatial grid

	Returns 10% percentiles of the sample
	*/
	vector<double> quantiles;
	int prog=0;
	double perc;
	for (int loc=0; loc<clds[0].nXSize * clds[0].nYSize; loc++)
	{

			int maskV=0;
			int nSamples=0;
			vector<double> out;
			vector<float> localCloudiness;
			double daysBetweenSamples;
			// only do if landmask is 1
			maskV = getMaskValue(msk, loc);
			nSamples = getSampleCounts(numSamples, loc);
			//cout << nSamples << endl;
			//cout << "mask is: " << maskV << endl;
			if ((maskV ==1) & (nSamples > 0.0))
			{
				cout << "processing: " << loc << " (Completed " << 100.0*prog/(clds[0].nXSize * clds[0].nYSize)<< "%%)" << endl;
				// retrieve cloudiness
				localCloudiness =  getCloudiness(clds, loc);
				// retrieve number of samples across the year
				/* convert the number of samples to a repeat rate
					eg in terms of days between sucessive samples	
					
					Now obviosuly this will not be an int 
					so we have a problem about how
					to implement this when only using sampling once per day

					-- my suggestion for now is that we can increase the size
					of the counting from one per day up to N per day quite easily
						-- this will give us the best accuracy possible really
				*/
				daysBetweenSamples = (double) 365.0/nSamples; 
				out = doSamplingModel(localCloudiness,
								 maskV, daysBetweenSamples);
				/*
				Calculate some statistics
				mean, std and 10 percentiles with gsl
				*/
				sort(out.begin(), out.end()); 
				double sum = accumulate(out.begin(), out.end(), 0.0);
				double mean = sum / out.size();	
				/*
				Calculate percentiles
				*/
				for (double p=0; p<1.1; p+=0.1) {
					perc = gsl_stats_quantile_from_sorted_data(out, 1,out.size(), p);
					cout << "percentile value: " << p << " " << perc << endl;
					quantiles.push_back(perc);
				}
				prog++;
			}
			else
			{
				/* 
				still need to add a value to ensure 
				ouput raster lines up nicely
				*/
				output.push_back(-999.0);
				prog++;
			}
	}

	return output;
} 

void saveDS(vector<double> outData, string filename) {
	/*
	Saves the ouput as a gdal array

	uses the mask as basis for output
	*/
	GDALRasterBand *inBand;
	GDALRasterBand *outBand;
	std::string srcfilename;
	srcfilename = "0.25_LW_mask.tif";
	const char * srcfc = srcfilename.c_str();
	// get size info
	int nXSize;
	int nYSize;
	GDALDataset *poSrcDS = (GDALDataset *) GDALOpen( srcfc, GA_ReadOnly );
	inBand = poSrcDS->GetRasterBand( 1 );
	nXSize = inBand->GetXSize();
	nYSize = inBand->GetYSize();
	double inGeoTransform[6];
	poSrcDS->GetGeoTransform( inGeoTransform );

	GDALDataset *poDstDS;
	const char *pszFormat = "GTiff";
    GDALDriver *poDriver;
    char **papszMetadata;
    poDriver = GetGDALDriverManager()->GetDriverByName(pszFormat);
    if( poDriver == NULL ) {
        exit( 1 );
    }
    //std::string dstfilename;
	//dstfilename = "mean.tif";
	const char * sdstfc = filename.c_str();

	char **papszOptions = NULL;
	poDstDS = poDriver->Create( sdstfc, nXSize, nYSize, 1, GDT_Float32, 
	                            papszOptions );
		
	//OGRSpatialReference oSRS;
	//char *pszSRS_WKT = NULL;
	//GDALRasterBand *poBand;


	poDstDS->SetGeoTransform( inGeoTransform );
	//oSRS.SetWellKnownGeogCS( "NAD27" );
	//oSRS.exportToWkt( &pszSRS_WKT );
	poDstDS->SetProjection( poSrcDS->GetProjectionRef() );


	// get size info
	// do some writing
	float ouputRaster[nXSize*nYSize];
	// fill the ouput array with input
	for (int i=0; i<nXSize*nYSize; i++) {
		ouputRaster[i]= outData[i];
	}
	outBand = poDstDS->GetRasterBand(1);
	outBand->RasterIO( GF_Write, 0, 0, nXSize, nYSize, 
	                  ouputRaster, nXSize, nYSize, GDT_Float32, 0, 0 );   
	/* Once we're done, close properly the dataset */
	//GDALClose( (GDALDatasetH) poDstDS );
	GDALClose( (GDALDatasetH) poDstDS );
	return;
}

int main(int argc, char* argv[])
 {
	
    if (argc != 3) {
        std::cerr << "* Observation Model *\n\n" 
                  << "Program which computes theoretical dob uncertainty from sampling and cloudiness.\n"
                  << "Outputs location and counts.\n\n"
                  << "Usage: " << argv[0] << " [ samples.tif ]" 
                  << "  [ out_file.tif ] "<< std::endl;
        return(1);
    }
    // associate command line args with files
    string samples_filename = argv[1];
    string outfilename = argv[2];
	/*
	Run the loader
	*/
	// prepare the input cloud datasets
	vector<cloudiness_month> clds;
	clds = prepareData();
	//cout << "Loaded cloud data" << endl;
	// Load land-water mask
	maskRaster msk;
	msk = loadMask();
	//cout << "m size:" << msk.nXSize << msk.nYSize << endl;
	//cout << "Loaded mask" << endl;
	// Load the sensor nSamples
	sampleRaster nSamples;
	nSamples = loadNumSamples(samples_filename);
	// try the sampler
	vector<double> results; // the bit we like
	// make a dummy results
	//for (int i=0; i<msk.nXSize*msk.nYSize; i++){
	//	results.push_back(i);
	//}
	results = runModel( clds, msk, nSamples);
	saveDS(results, outfilename);
	return 0;
}



