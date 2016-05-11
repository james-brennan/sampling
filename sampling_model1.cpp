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
		cout << months[m]<< endl; 
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
		//cout << data << endl;
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
		cout << "here" << endl;
		for (int i=0; i<nXSize*nYSize; i++)
		{
			current_month.cloudField.push_back(data[i]*0.001*0.1);

		}
		cloudiness_field.push_back(current_month);
		// free memory
		CPLFree(data);
		GDALClose(poDataset);
		cout << "done" << m << endl;
    }
    cout << "ready to return" << endl;
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


vector<int> doSamplingModel(vector<float> cldsXY,
					 int MaskValue, int sampleRepeatRate)
{
	/*
	The key part
	Does sampling on location x y
	*/
	const int N = 100; // number of iterations
	vector<int> daysSample;
	const int mdays[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
	int sample[365*3] = { }; ; // store the cloudiness

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
	    			sample[moffset]=1;
	    			sample[364+moffset]=1;
	    		    sample[364*2+moffset]=1;
	    		}
	    		else
	    		{
	    			sample[moffset]=0;
	    			sample[364+moffset]=0;
	    		    sample[364*2+moffset]=0; 
	    		}
	    	}
	    }
	    //cout << "made sample" << endl;
	    // print out the sample
	    // vary DOB
	    int dob = 55;
	    // count days before and after to an observation
	    int before=0;
	    int after=1;
	    while (before < ((365+dob)))
	    {
	    	//cout << sample[364+dob-before] << endl;
	    	// count days back before until a clear day
	    	if (sample[364+dob-before]==0) {
	    		before++;
	    		//cout << before << endl;
	    	}
	    	else {
	    		break;
	    	}
	    }
	    while (after < (dob+(365))) 
	    {
	    	// count days back before until a clear day
	    	if (sample[364+dob+after]==0) {
	    		after++;
	    	}
	    	else {
	    		break;
	    	}
	    }
	    // return days between before and after
	    int days;
	    days = after + before;
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

vector<double> runModel( vector<cloudiness_month> clds, maskRaster msk, 
					  int sampleRepeatRate) 
{
	/*
	This sets up the output and runs the model on the spatial grid
	*/
	vector<double> output;
	for (int loc=0; loc<clds[0].nXSize * clds[0].nYSize; loc++)
	{

			int maskV;
			int sampleRepeatRate = 1;
			vector<int> out;
			vector<float> localCloudiness;
			
			// only do if landmask is 1
			maskV = getMaskValue(msk, loc);
			//cout << "mask is: " << maskV << endl;
			if (maskV ==1)
			{
				cout << "processing: " << loc << endl;
				localCloudiness =  getCloudiness(clds, loc);
				maskV = getMaskValue(msk, loc);
				// cout << "land mask is: " << maskV << endl; 
				out = doSamplingModel(localCloudiness,
								 maskV, sampleRepeatRate);
				// get mean of the output
				// from
				// http://stackoverflow.com/questions/7616511/calculate-mean-and-standard-deviation-from-a-vector-of-samples-in-c-using-boos
				double sum = accumulate(out.begin(), out.end(), 0.0);
				double mean = sum / out.size();
				//for (int i=0; i<out.size(); i++)
				//{
				//	sum += out[i];
				//	count++;
				//}
				output.push_back(mean);
			}
			else
			{
			/* 
			still need to add a value to ensure 
			ouput raster lines up nicely
			*/
			output.push_back(-999.0);
			}
	}

	return output;
} 

void saveDS(vector<double> outData) {
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
    std::string dstfilename;
	dstfilename = "mean.tif";
	const char * sdstfc = dstfilename.c_str();

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

int main() {
	/*
	Run the loader
	*/
	// prepare the input cloud datasets
	vector<cloudiness_month> clds;
	clds = prepareData();
	cout << "Loaded cloud data" << endl;
	maskRaster msk;
	msk = loadMask();
	cout << "m size:" << msk.nXSize << msk.nYSize << endl;
	cout << "Loaded mask" << endl;
	// try the sampler
	vector<double> results; // the bit we like
	// make a dummy results
	//for (int i=0; i<msk.nXSize*msk.nYSize; i++){
	//	results.push_back(i);
	//}
	results = runModel( clds, msk, 1);
	saveDS(results);
	return 0;
}



