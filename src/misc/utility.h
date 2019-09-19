#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <string>
#include <vector>
#include <utility>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <queue>
#include <cmath>
#include <fstream>
#include <cassert>
#include <unistd.h>

#include "spectrum.h"

using namespace std;

class Utility
{
    public:

        //CONSTANTS for spectrum in .sptxt file
        const string NUMPEAKS_SPTXT = "NumPeaks";
        const string NAME_SPTXT = "Name";
        const string PRECURSORMZ_SPTXT = "PrecursorMZ";
        const string BINARYFILEOFFSET_SPTXT = "BinaryFileOffset";

        //CONSTANTS for spectrum in .msp file
        const string NUMPEAKS_MSP = "Num peaks";
        const string NAME_MSP = "Name";
        const string PRECURSORMZ_MSP = "Parent";
        //const string BINARYFILEOFFSET_MSP = "BinaryFileOffset";
        
        //CONSTANTS for SpectraST search result in .pep.xml file
        const string SPECTRUM_XML = "spectrum=\"";
        const string LIB_FILE_OFFSET_XML = "lib_file_offset";
        const string VALUE_XML = "value";

        //CONSTANTS for spectrum in .mgf file
        const string TITLE_MGF = "TITLE"; 
        const string CHARGE_MGF = "CHARGE";
        const string PEPMASS_MGF = "PEPMASS";
        const string BEGIN_IONS_MGF = "BEGIN IONS";
        const string END_IONS_MGF = "END IONS";


        //key: binary_file_offset
        //value: idx of a spectrum in 'spectra_unknown', extracted from .sptxt file
        unordered_map<long long, int> hash_offset;

        //http://stackoverflow.com/questions/19929681/c-global-variable-declaration
        unordered_map<int, Spectrum> spectra_all;

        unordered_map<int, Spectrum> spectra_unknown;

        //'spectrum query title' to 'lib_file_offset' in '.pep.xml' file
        unordered_map<string, long long> title_to_offset;
        //map<string, long long> title_to_offset;

        //'TITLE' to 'spectrum idx' in '.mgf' file
        unordered_map<string, int> title_to_spectrum_idx;

        /**
         * Read spectra from .sptxt file and store them in var 'spectra' 
         *
         *  unordered_map records with key 'spectrum id', value 'spectrum class instance'
        */
        unordered_map<int, Spectrum> readSpectraFromSptxt(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &inten_threshold, const double &scale, const double &bin_size, const bool &setBinarize, const string &opts4SpecEmbedCollision, const bool verbose=false);

        void convertSptxtToMGF(const string &inputFile, const string &mgfFile, const bool toFilter=false, const double inten_threshold=0);

        /**
         * Read spectra from .msp file and store them in var 'spectra_unknown' 
         *
         *  unordered_map records with key 'spectrum id', value 'spectrum class instance'
        */
        unordered_map<int, Spectrum> readSpectraFromMGF(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &inten_threshold, const double &scale, const double &bin_size, const bool &setBinarize,const string &opts4SpecEmbedCollision,  const bool verbose=false); 

        /*
         * read spectra from .mgf file while each spectrum keep topK peaks
         */
        unordered_map<int, Spectrum> readSpectraFromMGFKeepTopK(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &scale, const double &bin_size, const bool &setBinarize, const string &opts4SpecEmbedCollision, const int topK, const bool verbose=false); 
        /*
         * choose TOP K highest peaks among RANGE window
         */
        unordered_map<int, Spectrum> readSpectraFromMGFByTingChen(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &scale, const double &bin_size, const bool &setBinarize, const string &opts4SpecEmbedCollision, const int topK, const double window_size, const bool verbose=false); 

        /**
         * Read spectra from .msp file and store them in var 'spectra_unknown' 
         *
         *  unordered_map records with key 'spectrum id', value 'spectrum class instance'
        */
        unordered_map<int, Spectrum> readSpectraDBFromMGF(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &inten_threshold, const double &scale, const double &bin_size, const bool &setBinarize,const string &opts4SpecEmbedCollision,  const bool verbose=false); 

        void normalizeMGF(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &inten_threshold, const double &scale, const double &bin_size, const bool &setBinarize, const string &mgfFile); 

        /**
         * Read spectra from .msp file and store them in var 'spectra_unknown' 
         *
         *  unordered_map records with key 'spectrum id', value 'spectrum class instance'
        */
        unordered_map<int, Spectrum> readSpectraFromMSP(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &inten_threshold, const double &scale, const double &bin_size, const bool &setBinarize, const string &opts4SpecEmbedCollision); 

        /**
         * Read search result of SpectraST, and
         *      stores into unordered_map<spectrum title, lib_file_offset>
         *
         */
        unordered_map<string, long long> readSpectraSTXML(const string &inputFile, const bool verbose=false);
        //map<string, long long> readSpectraSTXML(const string &inputFile);

        /*
         * Normalize the intensity of a spectrum
         * Convert m/z into a interval starting from 0
         */
        vector<pair<int, double> > normalizeAndConvert(const vector<pair<double, double> > &data, const double &min_mz, const double &precision, const double &max_intensity, const double &scale, const bool &setBinarize, const string &opt4SpecEmbedCollision);
 

        vector<string> readStrFromFile(const string &inputFile);

        static int get_bin_by_mass_idx(const double &precursor_mz, const int &bin_size);

        static Spectrum merge2spectra(const Spectrum & s1, const Spectrum & s2);

        static Spectrum filterSpectrum(const Spectrum & s, const int topK, const int window_size);
};

#endif //UTILITY_H
