#ifndef NEAREST_NEIGHBOR_H
#define NEAREST_NEIGHBOR_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <queue>
#include <random>
#include <string>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include "spectrum.h"
#include "utility.h"

using namespace std;

class Nearest_Neighbor{
    public:
        /*
         * N:   # of random numbers
         * min: the minimum potentially generated value
         * max: the maximum potentially generated value
         */
        static vector<double> randomVec(int N, double min, double max, string fileName, const bool unitVec= false, const bool verbose=false);

        /*
         * L hash tables
         * m hash functions, per hash table
         */
        static vector<vector<vector<double> > > getLHashFunctions(const int L, const int m, const int randomSize, const double min_random, const double max_random, const bool unitVec = false, const bool verbose=false);


        /**
         * Input:
         *  lhs:    mass spectrum peak data, i.e. <m/z-to-int, normalized-intensity>
         *  rhs:    one hash function
         * Output:
         *  1, if dot product greater than 0
         *  0, else
         *
         */

        //cannot name this function hash, cauz http://stackoverflow.com/questions/15037705/reference-to-function-is-ambiguous
        static int hash_custom(const vector<pair<int, double> > &lhs, const vector<double> &rhs);

        static string hash_custom(const vector<pair<int, double> > &lhs, const vector<double> &rhs, const double b, const double r);


        /**
         * Input: 
         *  lhs: one hash table, i.e. m hash functions
         *  rhs: one spectrum
         * Output:
         *  an integer representing the hash value of the input spectrum
         *
         */
        static int hash_custom(const vector<vector<double> > &lhs, const vector<pair<int, double> > &rhs);

        static string hash_custom(const vector<vector<double> > &lhs, const vector<pair<int, double> > &rhs, const double b, const double r);

        /**
         *  Input:
         *      hash_table: map index to mass spectrum idx
         *      hash_functions: L hash funcions to get the key of hash_table 
         *      spectrum:   mass spectrum, with first value being mz-to-int, second normalized intensity 
         *
         *  Output:
         *      a list of indices of matched similar mass spectra 
         *
         */
        static vector<int> query(const vector<unordered_map<int, vector<unordered_map<int, vector<int> > > > > &hash_tables, const vector<vector<vector<double > > > &hash_functions, const Spectrum &spectrum, const double& mz_tol, const double &bin_size);

        static vector<int> query(const vector<unordered_map<int, vector<int> > > &hash_tables, const vector<vector<vector<double> > > &hash_functions, const Spectrum &spectrum, const double & mz_tol, const double &bin_size );
        

        static vector<int> query(const vector<unordered_map<string, vector<int> > > &hash_tables, const vector<vector<vector<double> > > &hash_functions, const Spectrum &spectrum, const double b, const double r, const double & mz_tol, const double &bin_size );

        //serves as a brute force method, return top K records and its similarity
        static vector<pair<double, int> > query(const unordered_map<int, Spectrum> &spectra_all, const Spectrum &spectrum, const double &mz_tol, const string &dist_measure, const int K, const bool filterByMZ = true, const bool filterByCharge = true);
        
        //all kinds of distances
        static double distance(const vector<pair<int, double> > &lhs, const vector<pair<int, double> > &rhs, const string &name);

        static double distance(const vector<vector<double> > &lhs, const vector<vector<double> > &rhs, const string &name);

        static vector<pair<int, double> >  normalizeToUnit(const vector<pair<int, double> > &peaks);

        static vector<pair<int, double> >  binAndNormalizeToUnit(const vector<pair<int, double> > &peaks, const int binSize, const double maxMZ, const double precision);     

};
#endif
