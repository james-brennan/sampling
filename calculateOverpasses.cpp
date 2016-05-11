
/*

    Script used to calculate overpasses

    Takes an of locations csv and swath csv
*/
#define _USE_MATH_DEFINES

#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>
#include <fstream>
#include <cmath>

#include "gdal_priv.h"
#include "cpl_conv.h" // for CPLMalloc()

using namespace std;

// Define some consts
const double R_earth = 6378137.0; // earth radius [m]

/*
    Store the swaths in a more convenient structure..
*/
struct swath 
{
    int orbit_id;
    vector<long> utc;
    vector<double> grLongitude;
    vector<double> grLatitude;
    vector<double> swWidth;
};

/*
    Store information about which pixels
    to process in a maskRaster
*/
struct maskRaster 
{
    int nXSize;
    int nYSize;
    double geoTransform[6];
    vector<float> mask;
};




vector<swath>  convert_to_swaths(vector< vector<double> >  &swaths2D) {
    /*
    Convert into a easier structure
    
    Approach:
    want to collect all rows and store into a struct
    until the id changes...
    */

    /*
    1. Count number of seperate swaths and store lenghts of each...
    */
    std::vector<swath> allSwaths;

    int id = 0;
    std::vector<int> id_counts;
    int entry_id=0;
    int id_number_elements = 0;
    // loop over each swath
    for (int i=0; i < swaths2D.size(); i++) {
        id_number_elements++;
        entry_id = swaths2D[i][3];
        if (id != entry_id) 
        {
            id++;
            id_counts.push_back(id_number_elements);
            // reset to zero
            id_number_elements = 0;
        }
    }
    int offset = 0;
    // make swaths
    for (int i=0; i<id; i++) {
        /*
        Create a swath...
        */
        swath aSwath = swath();
        aSwath.orbit_id = i;
        /*
        Collect the data...
        */ 
        int utc;
        double Longitude;
        double Latitude;
        double swathWidth;
        int number_of_entries = id_counts[i];
        for (int rec=0; rec < number_of_entries; rec++) {
            // collect entries
            utc = swaths2D[offset][0];
            Latitude = swaths2D[offset][1];
            Longitude = swaths2D[offset][2];
            swathWidth = swaths2D[offset][4];
            // add to the struct vectors...
            aSwath.utc.push_back(utc);
            aSwath.grLatitude.push_back(Latitude);
            aSwath.grLongitude.push_back(Longitude);
            aSwath.swWidth.push_back(swathWidth);
            // add one to offset
            offset++;
        }
        // add this swath to the swaths vec
        allSwaths.push_back(aSwath);
    }
    return allSwaths;
}



/*
    File loading stuff
    ==================
*/
vector< vector<double> > load_swath_file(string filename) {
    // Set up storage vector
    vector< vector<double> > swStore;
    vector<double> swRow;
    // load file
    ifstream stream(filename);
    // file elements
    string line;
    string item;
    // loop over lines
    for (string line; getline(stream, line); )
    {
        istringstream ln(line);
        while (getline(ln, item, ','))
        {
            // get item in row
            swRow.push_back(atof(item.c_str()));
        }
        swStore.push_back(swRow);
        swRow.clear();
    }
    return swStore;
}

vector< vector<double> > load_pixels_file(string filename) 
{
    // Set up storage vector
    vector< vector<double> > swStore;
    vector<double> swRow;
    // load file
    ifstream stream(filename);
    // file elements
    string line;
    string item;
    // loop over lines
    for (string line; getline(stream, line); ) 
    {
        istringstream ln(line);
        while (getline(ln, item, ','))
         { 
            // get field in line
            swRow.push_back(atof(item.c_str())); 
        }
        swStore.push_back(swRow);
        swRow.clear();
    }
    return swStore;
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
    double tmp[6] = {};
    poDataset->GetGeoTransform(LW_MASK.geoTransform);
    //LW_MASK.geoTransform = tmp;
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

/*

    Save output
*/
void saveDS(vector<int> outData, string outfilename) {
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
    dstfilename = outfilename;
    //dstfilename = "mean.tif";
    const char * sdstfc = dstfilename.c_str();

    char **papszOptions = NULL;
    poDstDS = poDriver->Create( sdstfc, nXSize, nYSize, 1, GDT_Int32, 
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
    int ouputRaster[nXSize*nYSize];
    // fill the ouput array with input
    for (int i=0; i<nXSize*nYSize; i++) {
        ouputRaster[i]= outData[i];
    }
    outBand = poDstDS->GetRasterBand(1);
    outBand->RasterIO( GF_Write, 0, 0, nXSize, nYSize, 
                      ouputRaster, nXSize, nYSize, GDT_Int32, 0, 0 );   
    /* Once we're done, close properly the dataset */
    //GDALClose( (GDALDatasetH) poDstDS );
    GDALClose( (GDALDatasetH) poDstDS );
    return;
}



/*
    Geodetic
    ========
*/

double haversine(double lat1, double lon1, double lat2, double lon2) {
    /*
    Haversine implementation to get distance between
    two locations in metres
    */
    // convert to radians
    lon1 = lon1 * (M_PI / 180.0); 
    lat1 = lat1 * (M_PI / 180.0); 
    lon2 = lon2 * (M_PI / 180.0); 
    lat2 = lat2 * (M_PI / 180.0); 
    // haversine formula 
    double dlon = lon2 - lon1;
    double dlat = lat2 - lat1;
    double a = pow(sin(dlat/2.0),2.0) + cos(lat1) * cos(lat2) * pow(sin(dlon/2.0),2.0);
    double c = 2 * asin(sqrt(a));
    double dist = c * R_earth;
    return dist;
}

double latitudeDistance(double lat0, double lat1) {
    // compute distance between two latitudes
    double dlat = abs(lat1 - lat0);
    double distance  = pow((lat1 - lat0),2.0);
    return distance;

}


int count_overpasses(double lat, double lon, std::vector<swath> &swaths ) {
    /*
        Calculate the number of overpasses for a location
    */
    int counts = 0;

    /* iterate through the swaths and find:
        1. nearest latitude location
        2. Whether the coord at [1.] is within 0.5 * swathWidth or this earthPoint.
    */

    std::vector<swath>::iterator sw;
    for (sw = swaths.begin(); sw != swaths.end(); ++sw) {
        /*
        Iterate through all the latitudes of the swath to find the nearest
        */
        int min_latitude_idx = 0;
        double dist = 0;
        double mindist = 1e12;
        int idx = 0;
        // loop over all latitudes in swath
        for (int j=0; j <  sw->grLatitude.size(); j++)
        {
            /*
                Find index of nearest latitude 
            */
            double orbitLat = sw->grLatitude[j];
            dist = latitudeDistance( lat, orbitLat );
            //cout << dist << endl;
            if (dist < mindist) 
            {
                // set min distance and set right latitude index...
                mindist = dist;
                min_latitude_idx = idx;
            }
            // increment index counter
            idx++;
        }
        /*
        Now have the minimum distance latitude location for this orbit
        Want to check whether this location to our point is within 0.5* swathWidth
        */
        //cout << "min index is: " << min_latitude_idx << endl;
        double orbitLatitude =  sw->grLatitude[min_latitude_idx];
        double orbitLongitude = sw->grLongitude[min_latitude_idx];
        /*
        Perform distance calculation
        */
        double distance;
        distance = haversine(lat, lon, orbitLatitude, orbitLongitude);
        //cout <<  "Metres distance between locations (" << lat << ", " << lon << ") and (" << orbitLatitude << "," << orbitLongitude << ") is: " << distance << endl;
        //cout << sw->orbit_id << " " << min_latitude_idx << " " << mindist << endl; 
        /*
        Check whether this distance is less than 0.5 swath
        */
        if (distance < 0.5 * sw->swWidth[min_latitude_idx]) {
            // add one to count
            counts++;
        }
    }
    return counts;
}

/*
    Main
    ====
*/

int main(int argc, char* argv[])
{
    
    if (argc != 3) {
        std::cerr << "* Count Overpasses *\n\n" 
                  << "Program which counts satellite overpasses for location.\n"
                  << "Outputs location and counts.\n\n"
                  << "Usage: " << argv[0] << " [ orbits_file ]" 
                  << "  [ out_file ] "<< std::endl;
        return(1);
    }
    // associate command line args with files
    string swathFile = argv[1];
    string outfilename = argv[2];

    // Load swath data file
    vector< vector<double> > swaths_f = load_swath_file(swathFile);

    // put into structures
    vector<swath> swaths = convert_to_swaths(swaths_f);  

    // Load the pixel locations
    // vector< vector<double> > locations = load_pixels_file(coordinates);
    // Load the LW mask
    maskRaster msk;
    msk = loadMask();

    vector<int> counts; // store the counts

    /*
    Do counts
    */
    double lat;
    double lon;
    int land;
    int count;
    int fill = -999;
    int i = 0;
    cout << "raster size is (x,y): " << msk.nXSize << " " << msk.nYSize << endl;
    for (int y=0; y < msk.nYSize; y++) 
    {
        for (int x=0; x < msk.nXSize; x++)
        {
            /*
            Calculate lat long location using the geoTransform
            */
            double yi = (double) y;
            double xi = (double) x;
            lon = msk.geoTransform[0] + msk.geoTransform[1]*xi + msk.geoTransform[2]*yi;
            lat = msk.geoTransform[3] + msk.geoTransform[4]*xi + msk.geoTransform[5]*yi;
            // want to measure from middle of the location too
            lon += 0.5 * msk.geoTransform[1];
            lat += 0.5 * msk.geoTransform[5];
            land = msk.mask[x + y * msk.nXSize];
            if (land == 1.0) {
                count =  count_overpasses(lat, lon, swaths);
                cout << lat << " " << lon <<  " "<<count << endl;
                counts.push_back(count);
            }
            else 
            {
                // is water so put a fill value
                cout << lat << " " << lon <<  " "<< fill << endl;
                counts.push_back(fill);
            }
            // increment i
            i++;
        }
    }
    /*
    Now let's save it out to file
    */
    cout << "did counts" << endl;
    saveDS(counts, outfilename);
    return 0;
}

