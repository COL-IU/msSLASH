#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <random>
#include <algorithm>
#include <map>
#include <unordered_set>
#include <vector>
#include <unistd.h>
#include <chrono>
#include <bitset>

#include "spectrum.h"
#include "utility.h"
#include "nearest_neightbor.h"


//#define DEBUG
//#define WRITE_IO
const int SEED = 4057218;

using namespace std;


int main(int argc, char *argv[])
{
    double cos_threshold = stof(argv[1]);
    double precision = stof(argv[2]);

    const double min_random = -1, max_random = 1;
    const bool setBinarize = false; //true for inten = 1, for peak, 0, otherwise
    const bool verbose = false; //save intermediate ./out/files
    const double min_mz = 0, max_mz = 1800, /*precision = 1,*/ inten_threshold = 0.01, scale = 10000, bin_size = 2, mass_tolerance = 0.05;
    const string opts4SpecEmbedCollision = "sum";
    const int randomSize = int((max_mz - min_mz) / precision) + 1; 

    const int topK_segment = 5; //have try 5, now try 10
    const double window_size = 100;

    const int L = 1, m = 12; //hash table, hash functions
    constexpr int iterations = 80;


    //double cos_threshold = 0.5;
    //double cos_threshold = 0.9;

    cout << "------PARAMS-----" << endl;
    cout << "minmz: " << min_mz << endl;
    cout << "maxmz: " << max_mz << endl;
    cout << "precision: " << precision << endl;
    cout << "cos_threshold: " << cos_threshold << endl; 
    cout << "mass_tolerance: " << mass_tolerance << endl;
    cout << "opts4specEmbedCollision: " << opts4SpecEmbedCollision << endl;
    cout << "window_size: " << window_size << endl;
    cout << "topK per window: " << topK_segment<< endl;
    cout << "look back steps to stop: " << iterations << endl;
    cout << "------PARAMS-----" << endl;


    auto start_time = chrono::high_resolution_clock::now();
    Utility utility;
    cout << "reading unknown spectra" << endl;
    utility.spectra_unknown = utility.readSpectraFromMGFByTingChen("./data/Adult_Liver.mgf", min_mz, max_mz, precision, scale, bin_size, setBinarize, opts4SpecEmbedCollision, topK_segment, window_size);


    //unordered_map<int, Spectrum> spectra_unknown = utility.spectra_unknown;

    unordered_map<int, Spectrum > spectra_unknown; //only contains +charge spectra
    unordered_map<int, int> relations;//map c+2 index to original index
    const int charge = 2;
    int counter = 0;
    for (size_t i = 0; i < utility.spectra_unknown.size(); ++i)
        if (utility.spectra_unknown[i].charge == charge)
        {
            spectra_unknown[counter] = utility.spectra_unknown[i];
            relations[counter] = i;
            ++counter;
        }

    const unordered_map<string, int> title2idx_un = utility.title_to_spectrum_idx;
    cout << "now we have unknown spectra:\t" << utility.spectra_unknown.size()  << endl;
    cout << "now we have c+2 unknown spectra:\t" << spectra_unknown.size()  << endl;
    auto end_time = chrono::high_resolution_clock::now();
    auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    cout << " used time:\t" << elapsed_read << endl;

    //ABOVE CORRECT, GUARANTEED, May 12, 2017

    unordered_set<int> inCluster; // keep track of spectrum already in a cluster.
    unordered_map<int, unordered_set<int> > cluster; //

    ofstream out;

    start_time = chrono::high_resolution_clock::now();

    bool clustered = false;
    double max_cos_sim;
    int s1_id, s2_id;
    int max_i, max_j;

    do
    {
        clustered = false;
        max_cos_sim = -1;
        s1_id = s2_id = -1;
        max_i = max_j = -1;

        size_t len = spectra_unknown.size();
        for (size_t i = 0; i < len; ++i)
        {
            if (inCluster.count(i))
                continue;

            const auto & curSpectrum = spectra_unknown.at(i); 
            for (size_t j = i + 1; j < len; ++j)
            {
                if(inCluster.count(j))
                    continue;

                cout << "i: " << i << " j: " << j << endl << flush;

                const auto & candSpectrum = spectra_unknown.at(j);
                if (abs(curSpectrum.precursor_mz - candSpectrum.precursor_mz) > mass_tolerance)//used to be 2
                    continue;

                double cos_sim = 1 - Nearest_Neighbor::distance(curSpectrum.peaks, candSpectrum.peaks, "cosine");

                cout << "pass precursor range" << endl;

                if (cos_sim >= cos_threshold && cos_sim > max_cos_sim)
                {
                    max_cos_sim = cos_sim;
                    s1_id = i, s2_id = j;
                    max_i = i, max_j = j;
                }
            }
        }//for

        if (max_cos_sim >= 0)
        {
            const Spectrum & s1 = spectra_unknown.at(s1_id); 
            const Spectrum & s2 = spectra_unknown.at(s2_id); 
            Spectrum s_new = Utility::merge2spectra(s1, s2);
            Spectrum s_new_filtered = Utility::filterSpectrum(s_new, topK_segment, window_size);

            int index = counter;
            spectra_unknown[index] = s_new_filtered;

            inCluster.insert(s1_id);
            inCluster.insert(s2_id);

            //merge into a cluster
            if (cluster.count(s1_id))
            {
                cluster[index].insert(cluster[s1_id].begin(), cluster[s1_id].end());
                cluster.erase(s1_id);
            }
            else
                cluster[index].insert(s1_id);
            if (cluster.count(s2_id))
            {
                cluster[index].insert(cluster[s2_id].begin(), cluster[s2_id].end());
                cluster.erase(s2_id);
            }
            else
                cluster[index].insert(s2_id);

            ++counter;

            clustered = true;
            
            cout << "merging " << s1_id << " & " << s2_id << " ->" << counter << endl << flush;
        }
    }while(clustered);

    end_time = chrono::high_resolution_clock::now();
    elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    cout << "exhaustive cluster time:\t" << elapsed_read << endl;

    cout << "write cluster to files" << endl;
    string res_file = "./out/exhaustive_cluster.txt.cos" + to_string(cos_threshold)+"pre"+to_string(precision);
    out.clear();
    out.open(res_file);

    out << "ID" << "\t" << "Titles" << endl;

    int cluster_cnt = 0;
    start_time = chrono::high_resolution_clock::now();

    for (size_t i = 0; i < spectra_unknown.size(); ++i)
    {
        if (inCluster.count(i))
            continue;

        if (cluster.count(i))
            for (const auto & ele : cluster[i])
                out << spectra_unknown[ele].title << ".dta;";
        else
            out << spectra_unknown[i].title << ".dta" ;

        out << endl;
        ++cluster_cnt;
    }

    cout << "total cluster num:\t" << cluster_cnt << endl; 
    out.close();
    end_time = chrono::high_resolution_clock::now();
    elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time).count();
    cout << "write IO time:\t" << elapsed_read << endl;
    return 0;
}

