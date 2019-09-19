#include "nearest_neightbor.h"

vector<double> Nearest_Neighbor::randomVec(int N, double min, double max, string fileName, const bool unitVec, const bool verbose)
{
    assert(N > 0 && min < max);


    vector<double> result(N, 0);    
    random_device rd;
    mt19937 gen(rd());
    //uniform_real_distribution<> dis(min, max);
    normal_distribution<> dis(0, 1); //go with Gaussian distribution first

    for (int i = 0; i < N; ++i)
    {
        double rand = dis(gen);
        result[i] = rand;
    }

    // normzlie this hash vector to unit length
    if (unitVec)
    {
        double total = 0;
        for(const auto & ele : result)
            total += pow(ele, 2);
        total = sqrt(total);

        for (int i = 0; i < N; ++i)
            result[i] /= total;
    }
    
    if (verbose)
    {
        ofstream fout(fileName);

        for (int i = 0; i < N; ++i)
            fout << setprecision(4) << result[i] << endl;
        fout.close();//close the handler
    }

    return result;
}

/*
 * L hash tables
 * m hash functions, per hash table
 */
vector<vector<vector<double> > > Nearest_Neighbor::getLHashFunctions(const int L, const int m, const int randomSize, const double min_random, const double max_random, const bool unitVec, const bool verbose)
{
    vector<vector<vector<double> > > hash_functions;
    for (int i = 0; i < L; ++i)//L hash tables
    {

        //construct m hash functions first
        vector<vector<double>> m_hash_functions;
        for (int j = 0; j < m; ++j)//m hash functions
        {
            string fileName = "hash_table_" + to_string(i) + "_hash_function_" + to_string(j);
            m_hash_functions.push_back(randomVec(randomSize, min_random, max_random, fileName, unitVec, verbose));
        }
        hash_functions.push_back(m_hash_functions);
    }
    
    return hash_functions;
}

//lhs peak vector, rhs hash function
string Nearest_Neighbor::hash_custom(const vector<pair<int, double> > &lhs, const vector<double> &rhs, const double b, const double r)
{
    double sum = 0;

    for (const auto & peak : lhs)
        sum += peak.second * rhs[peak.first];

    int hash_value = int((sum + b) / r);

    return to_string(hash_value);
}
string Nearest_Neighbor::hash_custom(const vector<vector<double> > &lhs, const vector<pair<int, double> > &rhs, const double b, const double r)
{
    string key = "";
    for (const auto & hash_function : lhs)
    {
        string single_value = hash_custom(rhs, hash_function, b, r);        
        key += single_value;
    }
    return key;
}

vector<int> Nearest_Neighbor::query(const vector<unordered_map<string, vector<int> > > &hash_tables, const vector<vector<vector<double> > > &hash_functions, const Spectrum &spectrum, const double b, const double r, const double & mz_tol, const double &bin_size )
{
    vector<int> similar_spectra;
    for (int i = 0; i < int(hash_tables.size()); ++i)
    {
        const auto &hash_table = hash_tables[i];
        const auto &m_hash_functions = hash_functions[i];
        string index = hash_custom(m_hash_functions, spectrum.peaks, b, r);    

        if (0 != hash_table.count(index))
        {
            const vector<int> &temp = hash_table.at(index);
            similar_spectra.insert(similar_spectra.end(), temp.begin(), temp.end());
        }
    }
    
    return similar_spectra;
}


int Nearest_Neighbor::hash_custom(const vector<pair<int, double> > &lhs, const vector<double> &rhs)
{
    double sum = 0;// result of dot product, then divided by |lhs| and |rhs| to normalize

    //normalize each vector first
    
    //comment out on night, May 9, 2017
    //vector<pair<int, double> > lhs_norm = normalizeToUnit(lhs); 
    const vector<pair<int, double>> lhs_norm = lhs;

    //make sure that hash function is unit vector
    for (const auto & peak : lhs_norm)
        sum += peak.second * rhs[peak.first];

    //return sum > 0 ? 1 : 0;
    
    // May 9, 2017
    //cout << sum << endl;
    
    return sum >= 0 ? 1 : 0;
}

int Nearest_Neighbor::hash_custom(const vector<vector<double> > &lhs, const vector<pair<int, double> > &rhs)
{
    int hash_value = 0;
    string key = "";
    for (const auto & hash_function : lhs)
    {
        int single_value = hash_custom(rhs, hash_function);        
        hash_value = hash_value * 2 + single_value;
        key += to_string(single_value);
    }

    return hash_value;
}

vector<int> Nearest_Neighbor::query(const vector<unordered_map<int, vector<int> > > &hash_tables, const vector<vector<vector<double> > > &hash_functions, const Spectrum &spectrum, const double & mz_tol, const double &bin_size )
{
    vector<int> similar_spectra;
    for (int i = 0; i < int(hash_tables.size()); ++i)
    {
        const auto &hash_table = hash_tables[i];
        const auto &m_hash_functions = hash_functions[i];
        int index = hash_custom(m_hash_functions, spectrum.peaks);    

        if (0 == hash_table.count(index))
            continue;
        const vector<int> &temp = hash_table.at(index);
        similar_spectra.insert(similar_spectra.end(), temp.begin(), temp.end());
    }
    
    return similar_spectra;
}

vector<int> Nearest_Neighbor::query(const vector<unordered_map<int, vector<unordered_map<int, vector<int> > > > > &hash_tables, const vector<vector<vector<double > > > &hash_functions, const Spectrum &spectrum, const double &mass_tolerance, const double &bin_size)
{
    vector<int> similar_spectra;
   
    for (int i = 0; i < int(hash_tables.size()); ++i)
    {
        const auto &hash_table = hash_tables[i];
        const auto &hash_function = hash_functions[i];
        int index = hash_custom(hash_function, spectrum.peaks);    

        int idx = spectrum.bin_by_mass_idx;

        //check charge first
        if (spectrum.charge >= (int)hash_table.at(index).size())
            continue;


        //for conciseness
        int chargeIdx = spectrum.charge - 1;
        const auto &curr_hash_table = hash_table.at(index).at(chargeIdx);

        int prevIdx = Utility::get_bin_by_mass_idx(spectrum.precursor_mz - mass_tolerance, bin_size);
        int nextIdx = Utility::get_bin_by_mass_idx(spectrum.precursor_mz + mass_tolerance, bin_size);
        
        for (int j = prevIdx; j <= nextIdx; ++j)
        {
            if (curr_hash_table.find(j) == curr_hash_table.end())
                continue;
            const vector<int> &temp = curr_hash_table.at(j);
            similar_spectra.insert(similar_spectra.end(), temp.begin(), temp.end());
        }
    }
    return similar_spectra; 
}
vector<pair<int, double>> Nearest_Neighbor::normalizeToUnit(const vector<pair<int, double> > &peaks)
{
    vector<pair<int, double> > res;
    double total = 0;
    for(const auto & peak : peaks)
        total += pow(peak.second, 2);
    total = sqrt(total);

    for(const auto & peak : peaks)
        res.push_back(make_pair(peak.first, peak.second / total) );

    return res;
}

//to verify tomorrow
vector<pair<int, double> >  Nearest_Neighbor::binAndNormalizeToUnit(const vector<pair<int, double> > &peaks, const int binSize, const double maxMZ, const double precision)
{
    int binNum = int(maxMZ / binSize) + 1; // ???
    double binRatio = precision / binSize;

    vector<pair<int, double> > newPeaks;
    vector<double> temp(binNum, 0); //temp peaks
    
    for (const auto &peak : peaks)
    {
        int bin_idx = peak.first * binRatio;
        temp[bin_idx] = 1;
    }
    
    for (int i = 0; i < binNum; ++i)
        if (1 == temp[i])
            newPeaks.push_back(make_pair(i, temp[i]));

    return normalizeToUnit(newPeaks);
}

double Nearest_Neighbor::distance(const vector<pair<int, double> > &lhs, const vector<pair<int, double> > &rhs, const string &name)
{
    //https://en.wikipedia.org/wiki/Cosine_similarity
    
    assert("cosine" == name || "dot_product" == name || "euclidean" == name);

    double distance = 0;
    double similarity = 0;

    if("euclidean" == name)
    {
        //normalize each vector first
        vector<pair<int, double> > lhs_norm = normalizeToUnit(lhs); 
        vector<pair<int, double> > rhs_norm = normalizeToUnit(rhs); 
        distance = 0;

        unordered_map<int, double> rhs_norm_map;//store rhs into this map for random access
        unordered_map<int, bool> rhs_visited;//idx is visited or not

        for (const auto & ele: rhs_norm)
        {
            rhs_norm_map[ele.first] = ele.second;
            rhs_visited[ele.first] = false;
        }
        
        //calculate
        for (const auto & ele : lhs_norm)
        {
            if(rhs_norm_map.find(ele.first) != rhs_norm_map.end()) // if found
            {
                distance += pow(ele.second - rhs_norm_map[ele.first], 2);
                rhs_visited[ele.first] = true;
            }
            else
                distance += pow(ele.second, 2);
        }
            
        for (const auto & ele : rhs_norm)
            if (!rhs_visited[ele.first])
                distance += pow(ele.second, 2);

        distance = sqrt(distance);

    }
    //else if ("dot_product" == name)
    else if ("dot_product" == name)
    {
        //normalize each vector first
        vector<pair<int, double> > lhs_norm = normalizeToUnit(lhs); 
        vector<pair<int, double> > rhs_norm = normalizeToUnit(rhs); 
        similarity = 0;

        unordered_map<int, double> rhs_norm_map;//store rhs into this map for random access

        for (const auto & ele: rhs_norm)
            rhs_norm_map[ele.first] = ele.second;

        for (const auto &ele : lhs_norm)
            similarity += ele.second * rhs_norm_map[ele.first];

        distance = 1 - similarity;

    }

    else if ("cosine" == name)
    {
        similarity = 0;
        double magnitude_lhs = 0, magnitude_rhs = 0;
        unordered_map<int, double> rhs_map;//store rhs into this map for random access

        for (const auto & ele: rhs)
        {
            magnitude_rhs += ele.second * ele.second;
            rhs_map[ele.first] = ele.second;
        }
        magnitude_rhs = sqrt(magnitude_rhs);

        for (const auto &ele : lhs)
        {
            magnitude_lhs += ele.second * ele.second;
            similarity += ele.second * rhs_map[ele.first];
        }
        magnitude_lhs = sqrt(magnitude_lhs);
        
        similarity = similarity / (magnitude_lhs * magnitude_rhs); 

        distance = 1 - similarity;
    }

    return distance;
}


//brute froce method of query, return top k
vector<pair<double, int>> Nearest_Neighbor::query(const unordered_map<int, Spectrum> &spectra_all, const Spectrum &spectrum, const double &mz_tol, const string &dist_measure, const int k, const bool filterByMZ, const bool filterByCharge)
{
    vector<pair<double, int> > res; //top k

    //max heap, higher distance get popped first
    auto cmp = [](const pair<double, int> &lhs, const pair<double, int> &rhs){return (lhs.first < rhs.first);};

    priority_queue<pair<double, int>, vector<pair<double, int> >, decltype(cmp)> que(cmp);

    for (const auto &ele : spectra_all)
    {
        const int idx = ele.first;
        const auto &candidate_spectrum = ele.second;

        //filter by charge, then mz_tol
        if (filterByCharge && candidate_spectrum.charge != spectrum.charge)
            continue;

        if (filterByMZ && (candidate_spectrum.precursor_mz < (spectrum.precursor_mz - mz_tol) ||
                candidate_spectrum.precursor_mz > (spectrum.precursor_mz + mz_tol)))
            continue;

        //then proceed to distance calc.
        double dist = distance(candidate_spectrum.peaks, spectrum.peaks, dist_measure);        

        que.push(make_pair(dist, idx));
        if(int(que.size()) > k)
            que.pop();
    }

    while(!que.empty())
    {
        res.push_back(que.top());
        que.pop();
    }

    reverse(res.begin(), res.end());

    return res;
}
