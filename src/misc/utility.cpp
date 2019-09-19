#include "utility.h"

//Plan to discard
//read search result from .xml file

unordered_map<string, long long> Utility::readSpectraSTXML(const string &inputFile, const bool verbose)
//map<string, long long> Utility::readSpectraSTXML(const string &inputFile)
{
    unordered_map<string, long long> title_to_offset;
    //map<string, long long> title_to_offset;

    ifstream infile(inputFile);
    string line, temp;
    string title;
    long long offset;
    int spectra_idx = 0;

    while(getline(infile, line))
    {
        //if found 'spectrum_query spectrum="' in .pep.xml
        if (string::npos != line.find(SPECTRUM_XML))
        {
            int idx = line.find(SPECTRUM_XML) + SPECTRUM_XML.length();
            temp = line.substr(idx, line.find_first_of('"', idx) - idx);
            title = temp.substr(0, temp.find_last_of('.'));
            continue;
        }
        //if not found lib_file_offset 
        if (string::npos == line.find(LIB_FILE_OFFSET_XML)) 
            continue;

        int idx = line.find(VALUE_XML) + VALUE_XML.length() + 1 + 1;// 1 for '=', 1 for '"'
        offset = stoll(line.substr(idx, line.find_first_of('"', idx) - idx));

        title_to_offset[title] = offset;
        ++spectra_idx;

        if (verbose)
            if (0 == spectra_idx % 5000)
                cout << spectra_idx << "...\t" << flush;
    }
    if (verbose)
        cout << spectra_idx << "...\t" << endl << flush;

    infile.close();
    return title_to_offset;
}

//read spectra from .mgf file
unordered_map<int, Spectrum> Utility::readSpectraFromMGF(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &inten_threshold, const double &scale, const double &bin_size, const bool &setBinarize, const string &opts4SpecEmbedCollision,  const bool verbose)
{
    //record title to spectrum_idx
    title_to_spectrum_idx.clear();

    unordered_map<int, Spectrum> spectra;

    vector<pair<double, double> > raw_peaks; //raw peaks
    vector<pair<double, double> > filtered_peaks; //filtered peaks

    ifstream infile(inputFile);
    string line, temp;
    string title;
    int charge;
    int spectrum_id = 0;
    double m_over_z, intensity, precursor_mz;
    string peptide;

    //long long binary_file_offset;
    istringstream iss;

    while(getline(infile, line))
    {
        if (0 == line.find(PEPMASS_MGF))
        {
            precursor_mz = stod(line.substr(PEPMASS_MGF.length() + 1)); // 1 for =''
            continue;
        }

        if (0 == line.find(CHARGE_MGF))
        {
            charge = stoi(line.substr(CHARGE_MGF.length() + 1)); //update stoi usage method
            continue;
        }


        //read until 'TITLE'
        if (string::npos == line.find(TITLE_MGF))
            continue;

        title = line.substr(line.find('=') + 1); 
        title = title.substr(0, title.find_last_of('.'));

        //read peaks of size 'num_peaks'
        double max_intensity = 0;
        double total_intensity = 0;
        raw_peaks.clear();
        while(getline(infile, line))
        {
            //if found end ions
            if (string::npos != line.find(END_IONS_MGF))
                break;

            iss.clear();
            iss.str(line);

            iss >> m_over_z >> intensity;
            
            //m/z has to be in the range [min_mz, max_mz];
            if (m_over_z < min_mz || m_over_z > max_mz)
                continue;

            total_intensity += intensity;
            max_intensity = max(max_intensity, intensity);

            raw_peaks.push_back(make_pair(m_over_z, intensity));
        }

        //filter each peak using a inten_threshold * total_intensity
        filtered_peaks.clear();
        for (const auto & peak : raw_peaks)
            if (peak.second >= (inten_threshold *total_intensity))
            {
                filtered_peaks.push_back(peak);
            }


        //normalize intensity && convert idx
        const vector<pair<int, double> > &normalizedPeaks= normalizeAndConvert(filtered_peaks, min_mz, precision, max_intensity, scale, setBinarize, opts4SpecEmbedCollision);

        Spectrum spectrum(normalizedPeaks, charge, -1, precursor_mz, get_bin_by_mass_idx(precursor_mz, bin_size));
        //Spectrum spectrum(normalizedPeaks, charge, -1, precursor_mz);
        //spectrum.bin_idx = get_bin_by_mass_idx(precursor_mz, bin_size);

        spectra[spectrum_id] = spectrum;
        title_to_spectrum_idx[title] = spectrum_id;
        ++spectrum_id;

        if (verbose)
            if (0 == spectrum_id % 5000)
                cout << spectrum_id << "...\t" << flush;
        //cout << endl << endl;
        //sleep(5);
    }
    if(verbose)
        cout << spectrum_id << "...\t" << endl << flush;

    infile.close();

    return spectra;
}

unordered_map<int, Spectrum> Utility::readSpectraFromMGFKeepTopK(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &scale, const double &bin_size, const bool &setBinarize, const string &opts4SpecEmbedCollision, const int topK, const bool verbose) 
{
    //record title to spectrum_idx
    title_to_spectrum_idx.clear();

    unordered_map<int, Spectrum> spectra;

    vector<pair<double, double> > raw_peaks; //raw peaks
    vector<pair<double, double> > filtered_peaks; //filtered peaks

    ifstream infile(inputFile);
    string line, temp;
    string title;
    int charge;
    int spectrum_id = 0;
    double m_over_z, intensity, precursor_mz;
    string peptide;

    //long long binary_file_offset;
    istringstream iss;

    while(getline(infile, line))
    {
        if (0 == line.find(PEPMASS_MGF))
        {
            precursor_mz = stod(line.substr(PEPMASS_MGF.length() + 1)); // 1 for =''
            continue;
        }

        if (0 == line.find(CHARGE_MGF))
        {
            charge = stoi(line.substr(CHARGE_MGF.length() + 1)); //update stoi usage method
            continue;
        }


        //read until 'TITLE'
        if (string::npos == line.find(TITLE_MGF))
            continue;

        title = line.substr(line.find('=') + 1); 
        title = title.substr(0, title.find_last_of('.'));

        //read peaks of size 'num_peaks'
        double max_intensity = 0;
        double total_intensity = 0;
        raw_peaks.clear();
        while(getline(infile, line))
        {
            //if found end ions
            if (string::npos != line.find(END_IONS_MGF))
                break;

            iss.clear();
            iss.str(line);

            iss >> m_over_z >> intensity;
            
            //m/z has to be in the range [min_mz, max_mz];
            if (m_over_z < min_mz || m_over_z > max_mz)
                continue;

            max_intensity = max(m_over_z, max_intensity);
            total_intensity += m_over_z;
            
            raw_peaks.push_back(make_pair(m_over_z, intensity));
        }

        //min heap, lower intensity get popped first
        auto cmp = [](const pair<double, double> &lhs, const pair<double, double> &rhs){return (lhs.second > rhs.second);};

        priority_queue<pair<double, double>, vector<pair<double, double> >, decltype(cmp)> que(cmp);

        for (const auto & peak : raw_peaks)
        {
            que.push(peak);
            if ((int)que.size() > topK)
                que.pop();
        }

        //extract topK in each bin.
        filtered_peaks.clear();
        while(!que.empty())
            filtered_peaks.push_back(que.top()), que.pop();
         
        sort(filtered_peaks.begin(), filtered_peaks.end()); //sort according to mz, i.e. 1st element, ascendingly


        //normalize intensity && convert idx
        const vector<pair<int, double> > &normalizedPeaks= normalizeAndConvert(filtered_peaks, min_mz, precision, max_intensity, scale, setBinarize, opts4SpecEmbedCollision);

        Spectrum spectrum(normalizedPeaks, charge, -1, precursor_mz, get_bin_by_mass_idx(precursor_mz, bin_size));
        //Spectrum spectrum(normalizedPeaks, charge, -1, precursor_mz);
        //spectrum.bin_idx = get_bin_by_mass_idx(precursor_mz, bin_size);

        spectra[spectrum_id] = spectrum;
        title_to_spectrum_idx[title] = spectrum_id;
        ++spectrum_id;

        /*
        cout << "spectrum_id:\t" << spectrum_id << endl;
        cout << title << endl;
        for (const auto & peak : filtered_peaks)
            cout << peak.first << "\t" << peak.second << endl;

        cout << endl;
        cout << "enter to continue ..." << endl;
        cin.ignore();
        */

        if (verbose)
            if (0 == spectrum_id % 5000)
                cout << spectrum_id << "...\t" << flush;
        //cout << endl << endl;
        //sleep(5);
    }
    if(verbose)
        cout << spectrum_id << "...\t" << endl << flush;

    infile.close();

    return spectra;
}

//read spectra from .mgf file, top K per RANGE window
unordered_map<int, Spectrum> Utility::readSpectraFromMGFByTingChen(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &scale, const double &bin_size, const bool &setBinarize, const string &opts4SpecEmbedCollision, const int topK, const double window_size, const bool verbose) 
{
    //record title to spectrum_idx
    title_to_spectrum_idx.clear();

    unordered_map<int, Spectrum> spectra;

    vector<pair<double, double> > raw_peaks; //raw peaks
    vector<pair<double, double> > filtered_peaks; //filtered peaks

    map<int, vector<pair<double, double> > > corpus; // used for extract topK among window_size

    ifstream infile(inputFile);
    string line, temp;
    string title;
    int charge;
    int spectrum_id = 0;
    double m_over_z, intensity, precursor_mz;
    string peptide;

    //long long binary_file_offset;
    istringstream iss;

    while(getline(infile, line))
    {
        if (0 == line.find(PEPMASS_MGF))
        {
            precursor_mz = stod(line.substr(PEPMASS_MGF.length() + 1)); // 1 for =''
            continue;
        }

        if (0 == line.find(CHARGE_MGF))
        {
            charge = stoi(line.substr(CHARGE_MGF.length() + 1)); //update stoi usage method
            continue;
        }


        //read until 'TITLE'
        if (string::npos == line.find(TITLE_MGF))
            continue;

        title = line.substr(line.find('=') + 1); 
        title = title.substr(0, title.find_last_of('.'));

        //read peaks of size 'num_peaks'
        double max_intensity = 0;
        double total_intensity = 0;
        raw_peaks.clear();
        corpus.clear();
        while(getline(infile, line))
        {
            //if found end ions
            if (string::npos != line.find(END_IONS_MGF))
                break;

            iss.clear();
            iss.str(line);

            iss >> m_over_z >> intensity;
            
            //m/z has to be in the range [min_mz, max_mz];
            if (m_over_z < min_mz || m_over_z > max_mz)
                continue;

            max_intensity = max(m_over_z, max_intensity);
            total_intensity += intensity;
            int idx = int(m_over_z / window_size);
            corpus[idx].push_back(make_pair(m_over_z, intensity));
        }

        //min heap, lower intensity get popped first
        auto cmp = [](const pair<double, double> &lhs, const pair<double, double> &rhs){return (lhs.second > rhs.second);};

        for (const auto & bin : corpus) 
        {
            priority_queue<pair<double, double>, vector<pair<double, double> >, decltype(cmp)> que(cmp);

            int idx = bin.first;
            const auto &peaks = bin.second;
            for (const auto & peak : peaks)
            {
                que.push(peak);
                if ((int)que.size() > topK)
                    que.pop();
            }

            //extract topK in each bin.
            vector<pair<double, double> > res; //<mz, intensity>
            while(!que.empty())
                res.push_back(que.top()), que.pop();
             
            sort(res.begin(), res.end()); //sort according to mz, i.e. 1st element, ascendingly
            raw_peaks.insert(raw_peaks.end(), res.begin(), res.end());
        }


        //normalize intensity && convert idx
        const vector<pair<int, double> > &normalizedPeaks= normalizeAndConvert(raw_peaks, min_mz, precision, max_intensity, scale, setBinarize, opts4SpecEmbedCollision);

        Spectrum spectrum(normalizedPeaks, charge, -1, precursor_mz, get_bin_by_mass_idx(precursor_mz, bin_size));
        
        //added on May 21th, 2017
        spectrum.title = title;

        //Spectrum spectrum(normalizedPeaks, charge, -1, precursor_mz);
        //spectrum.bin_idx = get_bin_by_mass_idx(precursor_mz, bin_size);

        spectra[spectrum_id] = spectrum;
        title_to_spectrum_idx[title] = spectrum_id;
        ++spectrum_id;

        /*
        cout << "spectrum_id:\t" << spectrum_id << endl;
        cout << title << endl;
        for (const auto & peak : raw_peaks)
            cout << peak.first << "\t" << peak.second << endl;

        cout << endl;
        cin.ignore();*/

        if (verbose)
            if (0 == spectrum_id % 5000)
                cout << spectrum_id << "...\t" << flush;
        //cout << endl << endl;
        //sleep(5);
    }
    if(verbose)
        cout << spectrum_id << "...\t" << endl << flush;

    infile.close();

    return spectra;
}

//read spectra from .mgf file
unordered_map<int, Spectrum> Utility::readSpectraDBFromMGF(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &inten_threshold, const double &scale, const double &bin_size, const bool &setBinarize, const string &opts4SpecEmbedCollision, const bool verbose)
{
    //record title to spectrum_idx
    title_to_spectrum_idx.clear();

    unordered_map<int, Spectrum> spectra;

    vector<pair<double, double> > raw_peaks; //raw peaks
    vector<pair<double, double> > filtered_peaks; //filtered peaks

    ifstream infile(inputFile);
    string line, temp;
    string title;
    int charge;
    int spectrum_id = 0;
    double m_over_z, intensity, precursor_mz;
    string peptide;

    istringstream iss;

    while(getline(infile, line))
    {
        if (0 == line.find(PEPMASS_MGF))
        {
            precursor_mz = stod(line.substr(PEPMASS_MGF.length() + 1)); // 1 for =''
            continue;
        }

        if (0 == line.find(CHARGE_MGF))
        {
            charge = stoi(line.substr(CHARGE_MGF.length() + 1)); //update stoi usage method
            continue;
        }


        //read until 'TITLE'
        if (string::npos == line.find(TITLE_MGF))
            continue;

        title = line.substr(line.find('=') + 1); 
        title = title.substr(0, title.find_last_of('.'));
        peptide = title.substr(0, title.find('/'));

        //read peaks of size 'num_peaks'
        double max_intensity = 0;
        double total_intensity = 0;
        raw_peaks.clear();
        while(getline(infile, line))
        {
            //if found end ions
            if (string::npos != line.find(END_IONS_MGF))
                break;

            iss.clear();
            iss.str(line);

            iss >> m_over_z >> intensity;

            if (m_over_z < min_mz || m_over_z > max_mz)
                continue;

            total_intensity += intensity;
            max_intensity = max(max_intensity, intensity);

            raw_peaks.push_back(make_pair(m_over_z, intensity));
        }

        //normalize intensity && convert idx
        const vector<pair<int, double> > &normalizedPeaks= normalizeAndConvert(raw_peaks, min_mz, precision, max_intensity, scale, setBinarize, opts4SpecEmbedCollision);

        Spectrum spectrum(normalizedPeaks, charge, -1, precursor_mz, get_bin_by_mass_idx(precursor_mz, bin_size), peptide);
        //Spectrum spectrum(normalizedPeaks, charge, -1, precursor_mz);
        //spectrum.bin_idx = get_bin_by_mass_idx(precursor_mz, bin_size);

        spectra[spectrum_id] = spectrum;
        title_to_spectrum_idx[title] = spectrum_id;
        ++spectrum_id;

        if (verbose)
            if (0 == spectrum_id % 5000)
                cout << spectrum_id << "...\t" << flush;
        //cout << endl << endl;
        //sleep(5);
    }
    if(verbose)
        cout << spectrum_id << "...\t" << endl << flush;

    infile.close();

    return spectra;
}

//read spectra from .mgf file and normalize the peaks to a new .mgf file
void Utility::normalizeMGF(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &inten_threshold, const double &scale, const double &bin_size, const bool &setBinarize, const string &mgfFile)
{
    ifstream infile(inputFile);
    ofstream outfile(mgfFile);

    vector<pair<double, double> > raw_peaks; //raw peaks
    vector<pair<double, double> > filtered_peaks; //raw peaks

    if(!outfile)
    {
        cout << "write file open error!" << endl;
        exit(-1);
    }

    string line, temp;
    string title;
    int spectrum_id = 0;
    double m_over_z, intensity, precursor_mz;
    string peptide;

    istringstream iss;

    while(getline(infile, line))
    {
        //read until 'TITLE'
        if (string::npos == line.find(TITLE_MGF))
        {
            outfile << line << endl;
            continue;
        }

        outfile << line << endl;

        //read peaks of size 'num_peaks'
        double max_intensity = 0;
        double total_intensity = 0;
        raw_peaks.clear();
        while(getline(infile, line))
        {
            //if found end ions
            if (string::npos != line.find(END_IONS_MGF))
                break;

            iss.clear();
            iss.str(line);

            iss >> m_over_z >> intensity;
            max_intensity = max(max_intensity, intensity);
            total_intensity += intensity;

            raw_peaks.push_back(make_pair(m_over_z, intensity));
        }

        //filter each peak using a inten_threshold * total_intensity
        //  and the m/z has to be in the range [min_mz, max_mz];
        filtered_peaks.clear();
        for (const auto & peak : raw_peaks)
            if (peak.second >= (inten_threshold * total_intensity) &&
                    peak.first >= min_mz && peak.first <= max_mz)
            {
                filtered_peaks.push_back(peak);
            }

        //normalize intensity && convert idx
        for (const auto &peak : filtered_peaks)
        {
            double norm_inten = peak.second * scale / max_intensity;
            outfile << peak.first << "\t" << norm_inten << endl;
        }
        outfile << "END IONS" << endl;
        ++spectrum_id;

        if (0 == spectrum_id % 5000)
            cout << spectrum_id << "...\t" << flush;
        //cout << endl << endl;
        //sleep(5);
    }
    cout << spectrum_id << endl;

    outfile.flush();
    outfile.close();
    infile.close();

    return;
}


//read spectra from .msp file
unordered_map<int, Spectrum> Utility::readSpectraFromMSP(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &inten_threshold, const double &scale, const double &bin_size, const bool &setBinarize, const string &opts4SpecEmbedCollision)
{
    unordered_map<int, Spectrum> spectra;

    vector<pair<double, double> > raw_peaks; //raw peaks
    vector<pair<double, double> > filtered_peaks; //filtered peaks

    ifstream infile(inputFile);
    string line, temp;
    int num_peaks, charge;
    string peptide;
    int spectrum_id = 0;
    double m_over_z, intensity, precursor_mz;
    //long long binary_file_offset;
    istringstream iss;

    while(getline(infile, line))
    {
        //read until spectra charge
        if (0 == line.find(NAME_MSP))
        {
            iss.clear();
            iss.str(line);
            iss >> temp >> temp;

            peptide = temp.substr(0, temp.find("/"));
            charge = stoi(temp.substr(temp.find("/") + 1));
            //cout << "charge: " << charge << endl;
            continue;
        }


        //read until spectra PrecursorMZ
        if (string::npos != line.find(PRECURSORMZ_MSP))
        {
            int idx = line.find(PRECURSORMZ_MSP) + PRECURSORMZ_MSP.length() + 1;// 1 for '='
            precursor_mz = stod(line.substr(idx, line.find_first_of(' ', idx) - idx));

            //cout << "precursor_mz: " << precursor_mz << endl;
            continue;
        }

        //read until 'NumPeaks'
        if (string::npos == line.find(NUMPEAKS_MSP))
            continue;

        iss.clear();
        iss.str(line);
        iss >> temp >> temp >> num_peaks; 

        //cout << "NUMPEAKS:\t" << num_peaks << endl;
        //read peaks of size 'num_peaks'
        double max_intensity = 0;
        double total_intensity = 0;
        raw_peaks.clear();
        while(num_peaks--)
        {
            getline(infile, line); 
            iss.clear();
            iss.str(line);

            iss >> m_over_z >> intensity;
            max_intensity = max(max_intensity, intensity);
            total_intensity += intensity;
            //cout << m_over_z << "\t" << intensity << endl;

            raw_peaks.push_back(make_pair(m_over_z, intensity));
        }

        //filter each peak using an inten_threshold * total_intensity
        //  and the m/z has to be in the range [min_mz, max_mz];
        filtered_peaks.clear();
        for (const auto & peak : raw_peaks)
            if (peak.second >= (inten_threshold * total_intensity) &&
                    peak.first >= min_mz && peak.first <= max_mz)
            {
                filtered_peaks.push_back(peak);
                //cout << peak.first << "\t" << peak.second << endl;
            }
        //cout << "max intensity:\t" << max_intensity << endl;


        //normalize intensity && convert idx
        const vector<pair<int, double> > &normalizedPeaks= normalizeAndConvert(filtered_peaks, min_mz, precision, max_intensity, scale, setBinarize, opts4SpecEmbedCollision);

        Spectrum spectrum(normalizedPeaks, charge, -1, precursor_mz, get_bin_by_mass_idx(precursor_mz, bin_size), peptide);
        //Spectrum spectrum(normalizedPeaks, charge, -1, precursor_mz);
        //spectrum.bin_idx = get_bin_idx(precursor_mz, bin_size);

        spectra[spectrum_id] = spectrum;

        ++spectrum_id;

        //cout << endl << endl;
    }

    infile.close();

    return spectra;
}

void Utility::convertSptxtToMGF(const string &inputFile, const string &mgfFile, const bool toFilter, const double inten_threshold)
{
    ifstream infile(inputFile);
    ofstream fout(mgfFile);

    string line, temp;
    string peptide;
    int num_peaks, charge, spectrum_id = 0;
    double m_over_z, intensity, precursor_mz;

    istringstream iss;

    while(getline(infile, line))
    {
        //read until spectra charge
        if (0 == line.find(NAME_SPTXT))
        {
            iss.clear();
            iss.str(line);
            iss >> temp >> temp;

            peptide = temp.substr(0, temp.find("/"));
            charge = stoi(temp.substr(temp.find("/") + 1));
            continue;
        }

        //read until spectra PrecursorMZ
        if (0 == line.find(PRECURSORMZ_SPTXT))
        {
            iss.clear();
            iss.str(line);
            iss >> temp >> precursor_mz;
            continue;
        }

        //read until 'NumPeaks'
        if (string::npos == line.find(NUMPEAKS_SPTXT))
            continue;

        iss.clear();
        iss.str(line);
        iss >> temp >> num_peaks; 

        fout << "BEGIN IONS" << endl;
        fout << "TITLE=" << spectrum_id++ <<",sequence="<<peptide << endl;
        fout << "PEPMASS="<<precursor_mz<<endl;
        fout << "CHARGE="<<charge<<"+"<<endl;

        if (spectrum_id% 5000 == 0)
            cout << spectrum_id << "...\t" << flush;

        //if no filter
        if (!toFilter)
            //read peaks of size 'num_peaks'
            while(num_peaks--)
            {
                getline(infile, line); 
                iss.clear();
                iss.str(line);

                iss >> m_over_z >> intensity;
                fout << m_over_z << "\t" << intensity << endl;
            }
        else //else filter
        {
            double total_intensity = 0;
            vector<pair<double, double> > raw_peaks;
            while(num_peaks--)
            {
                getline(infile, line); 
                iss.clear();
                iss.str(line);

                iss >> m_over_z >> intensity;
                total_intensity += intensity;
                raw_peaks.push_back(make_pair(m_over_z, intensity));
            }

            for (const auto &peak : raw_peaks)
                if (peak.second >= inten_threshold * total_intensity)
                    fout << peak.first << "\t" << peak.second << endl;
        }

        fout << "END IONS" << endl;
        fout << endl;
    }
    cout << spectrum_id << endl;

    fout.close();
    infile.close();
}


//read spectra from .sptxt file
unordered_map<int, Spectrum> Utility::readSpectraFromSptxt(const string &inputFile, const double &min_mz, const double &max_mz, const double &precision, const double &inten_threshold, const double &scale, const double &bin_size, const bool &setBinarize, const string &opts4SpecEmbedCollision, const bool verbose)
{
    hash_offset.clear();

    unordered_map<int, Spectrum> spectra;

    vector<pair<double, double> > raw_peaks; 
    vector<pair<double, double> > filtered_peaks; 

    ifstream infile(inputFile);

    string line, temp;
    string peptide;
    int num_peaks, charge, spectrum_id = 0;
    double m_over_z, intensity, precursor_mz;
    long long binary_file_offset;

    istringstream iss;

    while(getline(infile, line))
    {
        //read until spectra charge
        if (0 == line.find(NAME_SPTXT))
        {
            iss.clear();
            iss.str(line);
            iss >> temp >> temp;

            peptide = temp.substr(0, temp.find("/"));
            charge = stoi(temp.substr(temp.find("/") + 1));
            continue;
        }

        //read until spectra PrecursorMZ
        if (0 == line.find(PRECURSORMZ_SPTXT))
        {
            iss.clear();
            iss.str(line);
            iss >> temp >> precursor_mz;
            continue;
        }

        //read until spectra BinaryFileOffset
        if (string::npos != line.find(BINARYFILEOFFSET_SPTXT))
        {
            int idx = line.find(BINARYFILEOFFSET_SPTXT) + BINARYFILEOFFSET_SPTXT.length() + 1;// 1 for '='
            binary_file_offset = stoll(line.substr(idx, line.find_first_of(' ', idx) - idx));

            continue;
        }

        //read until 'NumPeaks'
        if (string::npos == line.find(NUMPEAKS_SPTXT))
            continue;

        iss.clear();
        iss.str(line);
        iss >> temp >> num_peaks; 

        //cout << "NUMPEAKS:\t" << num_peaks << endl;
        //read peaks of size 'num_peaks'
        double max_intensity = 0;
        double total_intensity = 0;
        raw_peaks.clear();
        while(num_peaks--)
        {
            getline(infile, line); 
            iss.clear();
            iss.str(line);

            iss >> m_over_z >> intensity;
            max_intensity = max(max_intensity, intensity);
            total_intensity += intensity;

            raw_peaks.push_back(make_pair(m_over_z, intensity));
        }

        //filter each peak using a inten_threshold * total_intensity
        //  and the m/z has to be in the range [min_mz, max_mz];
        filtered_peaks.clear();
        for (const auto & ele : raw_peaks)
            if (ele.second >= (inten_threshold * total_intensity) &&
                    ele.first >= min_mz && ele.first <= max_mz)
            {
                filtered_peaks.push_back(ele);
            }


        //normalize intensity && convert idx
        vector<pair<int, double> > normalizedPeaks= normalizeAndConvert(filtered_peaks, min_mz, precision, max_intensity, scale, setBinarize, opts4SpecEmbedCollision);

        Spectrum spectrum(normalizedPeaks, charge, binary_file_offset, precursor_mz, get_bin_by_mass_idx(precursor_mz, bin_size), peptide);

        //spectrum.bin_by_mass_idx = get_bin_by_mass_idx(precursor_mz, bin_size);

        spectra[spectrum_id] = spectrum;
        hash_offset[binary_file_offset] = spectrum_id;
        ++spectrum_id;

        if (verbose)
            if (0 == spectrum_id % 10000)
                cout << spectrum_id << "...\t" << flush;
        //cout << endl << endl;
        //sleep(5);
    }
    if (verbose)
        cout << spectrum_id << "...\t" << endl << flush;

    infile.close();

    return spectra;
}


/*
 * Normalize the intensity of a spectrum, to a max of 'scale'
 * Convert m/z into a interval starting from min_mz
 */
vector<pair<int, double> > Utility::normalizeAndConvert(const vector<pair<double, double> > &peaks, const double &min_mz, const double &precision, const double &max_intensity, const double &scale, const bool &setBinarize, const string & opts4SpecEmbedCollision)
{
    //when spectra embedded into same entry of high dimension, choose either 'max' inten or 'sum' all intens
    assert(opts4SpecEmbedCollision == "max" || opts4SpecEmbedCollision == "sum");

    vector<pair<int, double> > normalizedPeaks;
    double mz, inten;
    int idx;

    for (auto & peak : peaks) 
    {
        mz = peak.first, inten = peak.second;
        idx = round((mz - min_mz) / precision);

        if (!normalizedPeaks.empty() && normalizedPeaks.back().first == idx)
        {
            if ("sum" == opts4SpecEmbedCollision)
                normalizedPeaks.back().second += inten;
            else if ("max" == opts4SpecEmbedCollision)
                normalizedPeaks.back().second = max(normalizedPeaks.back().second, inten);
        }
        else
            normalizedPeaks.push_back(make_pair(idx, inten));

    }
    double max_inten = 0;
    for (auto & peak : normalizedPeaks)
        max_inten = max(max_inten , peak.second);

    for (auto &peak : normalizedPeaks)
    {
        if (setBinarize)
            peak.second = 1.;
        else
            peak.second = (peak.second / max_inten) * scale;
    }
    /*for (auto & peak : peaks) 
    {
        mz = peak.first, inten = peak.second;
        idx = round((mz - min_mz) / precision);

        if (setBinarize)
            inten = 1.;
        else
            inten = (inten / max_intensity) * scale;

        if (!normalizedPeaks.empty() && normalizedPeaks.back().first == idx)
        {
            if ("sum" == opts4SpecEmbedCollision)
            {
                normalizedPeaks.back().second += inten;
                if (setBinarize)
                    normalizedPeaks.back().second = 1.;
            }
            else if ("max" == opts4SpecEmbedCollision)
            {
                normalizedPeaks.back().second = max(normalizedPeaks.back().second, inten);
            }
        }
        else
            normalizedPeaks.push_back(make_pair(idx, inten));

    }*/
    return normalizedPeaks;
}

int Utility::get_bin_by_mass_idx(const double &precursor_mz, const int &bin_size)
{
    return round(precursor_mz / bin_size);
}

vector<string> Utility::readStrFromFile(const string &inputFile)
{
    ifstream infile(inputFile);
    vector<string> contents;
    string line;

    while(getline(infile, line))
        if (!line.empty() && "" != line)
            contents.push_back(line);

    infile.close();
    return contents;
}

Spectrum Utility::merge2spectra(const Spectrum &s1, const Spectrum &s2)
{
    assert(s1.charge == s2.charge);

    Spectrum s_new;
    s_new.charge = s1.charge; //s1.charge == s2.charge, should be
    s_new.bin_by_mass_idx = -1;
    s_new.binary_file_offset = -1;
    s_new.precursor_mz = (s1.precursor_mz + s2.precursor_mz) / 2.;
    s_new.peptide = "";

    size_t ii1 = 0, ii2 = 0;
    size_t s1size = s1.peaks.size(), s2size = s2.peaks.size();
    
    while(ii1 < s1size || ii2 < s2size)
    {
        if (ii2 == s2size)
        {
            s_new.peaks.push_back(s1.peaks[ii1]);
            ++ii1;
            continue;
        }
        if (ii1 == s1size)
        {
            s_new.peaks.push_back(s2.peaks[ii2]);
            ++ii2;
            continue;
        }
        if (s1.peaks[ii1].first < s2.peaks[ii2].first)
        {
            s_new.peaks.push_back(s1.peaks[ii1]);
            ++ii1;
            continue;
        }
        if (s1.peaks[ii1].first > s2.peaks[ii2].first)
        {
            s_new.peaks.push_back(s2.peaks[ii2]);
            ++ii2;
            continue;
        }
        
        s_new.peaks.push_back(make_pair(s1.peaks[ii1].first, s1.peaks[ii1].second + s2.peaks[ii2].second));
        ++ii1;
        ++ii2;
    }

    return s_new;
}

Spectrum Utility::filterSpectrum(const Spectrum & s, const int topK, const int window_size)
{
    Spectrum s_new(s); 
    map<int, vector<pair<int, double> > > corpus;
    for (const auto & peak : s.peaks)
    {
        int idx = int(peak.first / window_size);
        corpus[idx].push_back(peak);
    }

    //min heap, lower intensity get popped first
    auto cmp = [](const pair<int, double> &lhs, const pair<int, double> &rhs){return (lhs.second > rhs.second);};

    s_new.peaks.clear();
    for (const auto & bin : corpus) 
    {
        priority_queue<pair<int, double>, vector<pair<int, double> >, decltype(cmp)> que(cmp);

        int idx = bin.first;
        const auto &peaks = bin.second;
        for (const auto & peak : peaks)
        {
            que.push(peak);
            if ((int)que.size() > topK)
                que.pop();
        }

        //extract topK in each bin.
        vector<pair<int, double> > res; //<mz-to-int, intensity>
        while(!que.empty())
            res.push_back(que.top()), que.pop();
         
        sort(res.begin(), res.end()); //sort according to mz, i.e. 1st element, ascendingly
        s_new.peaks.insert(s_new.peaks.end(), res.begin(), res.end());
    }

    return s_new;
}
