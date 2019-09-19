#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <map>
#include <vector>


#include "spectrum.h"
#include "utility.h"
#include "nearest_neightbor.h"

using namespace std;

int main()
{
    
    //CONSTANTS
    //const int randomSize = 100000;
    const int randomSize = 100001;
    const double min_random = -1, max_random = 1;
    const int L = 30; //L hash tables
    const int m = 9; //m hash functions per table
    const bool setBinarize = true; //inten = 1, for peak, 0, otherwise
    const bool hashFuncUnitVec = true;// hash random spectrum unit length or not
    const bool verbose = false; //save intermediate files

    int counter = 0;

    //preprocessing for getting L hash tables
    cout << "generating hash functions ..." << endl;
    const vector<vector<vector<double> > > &hash_functions = Nearest_Neighbor::getLHashFunctions(L, m, randomSize, min_random, max_random, hashFuncUnitVec, verbose);

     
    //read in a list of class instance of SPECTRUM, as a dataset
    
    const double min_mz = 0, max_mz = 2000, precision = 0.02, inten_threshold = 0.01, scale = 10000, bin_size = 2, mass_tolerance = 0.5;

    Utility utility;
    
    //normalize peak's intensity to a scale
    //utility.normalizeMGF("./data/Adult_Liver_Gel_Velos_11_f01.mgf", min_mz, max_mz, precision, inten_threshold, scale, bin_size, setBinarize, "./data/Adult_Liver_Gel_Velos_11_f01_norm.mgf");
    //cout << "press CTRL+Z now!!" << endl;

    cout << "reading Human.sptxt ..." << endl;
    //cout << "reading Human10k.sptxt ..." << endl;
    
    utility.spectra_all = utility.readSpectraFromSptxt("./data/Human.sptxt", min_mz, max_mz, precision, inten_threshold, scale, bin_size, setBinarize);
    const unordered_map<int, Spectrum> &spectra_all = utility.spectra_all;
    const unordered_map<long long, int> &hash_offset = utility.hash_offset;

    //cout << "converting .sptxt to .mgf..." << endl;
    //utility.convertSptxtToMGF("./data/Human.sptxt", "./data/HumanFromSPTXT_v2.mgf", true, 0.01);
    //cout << "conversion done!!" << endl;

    cout << endl << endl;
    cout << "Human.sptxt size:\t" << spectra_all.size() << endl;
    cout << "hash_offset size:\t" << hash_offset.size() << endl;


    //vector<unordered_map<int, vector<unordered_map<int, vector<int> > > > > hash_tables(L); // L hash tables, where in each hash table key is hash_custom(dot product), value is vector<vector<spectrum ids with same hashed key> >, where the outerloop contains charge, and spectrum ids
    vector<unordered_map<int, vector<int> > > hash_tables(L); 
    cout << endl << endl;
    cout << "mapping known spectra..." << endl;
    counter = 0;
    for (const auto &spectrum : spectra_all)
    {
        //loop through each hash function, and map accordingly
        for (int i = 0; i < int(hash_functions.size()); ++i)
        {
            int idx = Nearest_Neighbor::hash_custom(hash_functions[i], spectrum.second.peaks);     
            hash_tables[i][idx].push_back(spectrum.first);
        }

        if (0 == (++counter % 50000))
            cout << counter << "...\t" << flush;
    }
    cout << counter << flush;
    
    cout << endl << endl;

    //  and put spectra into a map <hash_custom, mass spectrum id>
    
    /*
    cout << endl << endl;
    cout << "mapping known spectra ..." << endl;
    counter = 0;
    for (const auto &spectrum : spectra)
    {
        //loop through each hash function, and map accordingly
        for (int i = 0; i < int(hash_functions.size()); ++i)
        {
            int idx = Nearest_Neighbor::hash_custom(hash_functions[i], spectrum.second.peaks);     

            int charge = spectrum.second.charge;

            assert(charge > 0 && charge < 10);
            
            if (charge > (int)hash_tables[i][idx].size())
                hash_tables[i][idx].resize(charge);

            //hash_tables[i][idx][charge - 1].push_back(spectrum.first);
            hash_tables[i][idx][charge - 1][spectrum.second.bin_by_mass_idx].push_back(spectrum.first);
        }
        if (0 == (++counter % 50000))
            cout << counter << "...\t" << flush;
    }
    cout << counter << flush;
    
    cout << endl << endl;
    */

    cout << "reading unknown spectra ..." << endl;
    //read in a list of unknown spectra
    utility.spectra_unknown = utility.readSpectraFromMGF("./data/Adult_Liver_Gel_Velos_11_f01_norm.mgf", min_mz, max_mz, precision, inten_threshold, scale, bin_size, setBinarize);

    const unordered_map<int, Spectrum> &spectra_unknown = utility.spectra_unknown;
    const unordered_map<string, int> title_to_spectrum_idx = utility.title_to_spectrum_idx;

    cout << endl << endl;
    cout << "unknown spectra in .mgf size:\t" << spectra_unknown.size() << endl;

    //read in the search result of SpectraST
    utility.title_to_offset = utility.readSpectraSTXML("./data/Adult_Liver_Gel_Velos_11_f01.pep.xml");
    const unordered_map<string, long long> &title_to_offset = utility.title_to_offset;
    //const map<string, long long> &title_to_offset = utility.title_to_offset;

    cout << endl << endl;
    cout << "spectrast result file #search:\t" << title_to_offset.size() << endl;
    
    cout << endl << endl;

    cout << "querying ..." << endl;
    //query unknown spectra 
    int hit = 0;
    int brute_hit = 0;
    counter = 0;
    int continued = 0;
    vector<string> missed;

    
    ofstream fout("liver_brute_output_binary_cos70p.txt");
    for (const auto &search_query : title_to_offset)
    {
        const string search_title = search_query.first;
        const long long search_offset = search_query.second;

        //cout << search_title << endl;

        if (title_to_spectrum_idx.find(search_title) == title_to_spectrum_idx.end())
        {
            cout << "search title not in title_to_spectrum_idx!!" << endl;
            exit(-1);
        }

        const int search_spec_idx = title_to_spectrum_idx.at(search_title);
        //cout << "search_spec_idx: " << search_spec_idx << endl;
        
        const auto &unknown_spectrum = spectra_unknown.at(search_spec_idx);

        int charge = unknown_spectrum.charge;
        assert(charge > 0 && charge < 10);

        /*
        const vector<int> &similar_spectra = Nearest_Neighbor::query(hash_tables, hash_functions, unknown_spectrum, mass_tolerance, bin_size);

        cand_num.push_back(similar_spectra.size());
        candidate_num += similar_spectra.size();
        
        //http://stackoverflow.com/questions/5437643/strange-behaviour-of-stdcout-in-linux
        if (0 == ++counter % 500)
            cout << counter << "..." << flush;


        if (find(similar_spectra.begin(), similar_spectra.end(), hash_offset.at(search_offset)) != similar_spectra.end())
            ++hit;
        else
            missed.push_back(search_title);
            */

        vector<pair<double, int> >  similar_cos = Nearest_Neighbor::query(spectra_all, unknown_spectrum, mass_tolerance,  "cosine", 10); 

        if (0 == ++counter % 500)
            cout << counter << "..." << flush;
        
        if(!similar_cos.empty() && similar_cos.back().first > 0.7)
        {
            ++continued;
            continue;
        }

        fout << search_title << endl;

        /*
        fout << "brute force most similar idx:" << endl;
        fout << "COSINE, == NORMALIZED DOT PRODUCT" << endl;
        for (auto &pr : similar_cos)
            fout << "idx:\t" << pr.second << " cosine distance:\t" << pr.first << endl;*/
        
        /*fout << "EUCLIDEAN" << endl;
        vector<pair<double, int> >  similar_euc = Nearest_Neighbor::query(spectra_all, unknown_spectrum, mass_tolerance,  "euclidean", 10); 
        for (auto &pr : similar_euc)
            fout << "idx:\t" << pr.second << " euclidean distance:\t" << pr.first << endl;
        fout << "spectraST idx:\t" << hash_offset.at(search_offset)  << ", with euclidean distance" << Nearest_Neighbor::distance(spectra_all.at(hash_offset.at(search_offset)).peaks, unknown_spectrum.peaks, "euclidean") << endl;*/

        fout << "spectraST idx:\t" << hash_offset.at(search_offset)  << ", with cosine distance" << Nearest_Neighbor::distance(spectra_all.at(hash_offset.at(search_offset)).peaks, unknown_spectrum.peaks, "cosine") << endl;


        const vector<int> &similar_spectra = Nearest_Neighbor::query(hash_tables, hash_functions, unknown_spectrum, mass_tolerance, bin_size);

        fout << "lsh algorithm has cand #:\t" << similar_spectra.size() << endl;
        for (const auto & ele : similar_spectra)
            fout << ele <<", ";
        fout << endl << endl;

        if(find(similar_spectra.begin(), similar_spectra.end(), hash_offset.at(search_offset)) != similar_spectra.end())
            ++hit;


        for (const auto & pr : similar_cos)
            if (hash_offset.at(search_offset) == pr.second)
            ++brute_hit;
    }

    cout << endl << endl;
    cout << "hit #:\t" << hit << endl;
    cout << "brute_hit #:\t" << brute_hit << endl;
    cout << "total #:\t" << (spectra_unknown.size() - continued)<< endl;


    /*
    cout << endl << endl;
    cout << L << " hash tables, with " << m << " hash functions each" << endl; 
    cout << "total # unknown spectra:\t" << spectra_unknown.size() << endl;
    cout << "hits:\t" << hit << endl; 
    cout << "hit ratio:\t" << (hit * 1.0/spectra_unknown.size()) << endl; 
    cout << "avg # candidates per query: " << candidate_num * 1. / spectra_unknown.size() << endl; 


    ofstream fout2("missed.txt");
    if (!missed.empty())
    {
        fout2 << "missed spectra title listed as follows: " << endl;
        for (unsigned i = 0; i < missed.size(); ++i)
            fout2 << missed[i] << endl;
    }

    fout2.close();
    fout.close();
    cout << endl;

    */

    return 0;
}


