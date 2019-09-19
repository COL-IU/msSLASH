
#include <chrono>
#include <iostream>
#include <map>
#include <mutex>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <string>
#include <thread>

#include "../class/params.h"
#include "../class/spectrum.h"
#include "../class/core.h"
#include "../class/distance.h"
#include "../class/hyperparams.h"
#include "../utility/io.h"
#include "../utility/cmdparser.h"
#include "../utility/commons.h"

using namespace Core;
using namespace Utility;
using namespace std;

void SetParameters(const cli::Parser& parser, Core::HyperParams* params, vector<string>* files) {
  (*params).peak_intensity_rescale_method = parser.get<int>("r");
  (*params).precision = parser.get<float>("p");
  (*params).min_mz = parser.get<float>("l");
  (*files) = parser.get<vector<string> >("f");
}

void configure_parser(cli::Parser& parser) {
  parser.set_optional<int>("r", "rescale", 2, "[Int] PEAK_INTENSITY_RESCALE_METHOD (0:None | 1:Log | 2:Sqrt)");
  parser.set_optional<float>("p", "precision", 0.8, "[Float] Fragment precision");
  parser.set_optional<float>("l", "left_mz", 200, "[Float] left-most peak(mz) to consider");
  parser.set_required<vector<string> >("f", "files", "[Strs] Files holding MS2 spectra in .mgf format");
}

void readTSV(string fileName, unordered_map<string, string>& map_title_to_peptide, unordered_map<string, vector<string>>& map_peptide_to_titles) {
  ifstream fp(fileName);
  if (!fp.is_open()) {
    cout << "Error! Can not open" << fileName << endl;
    exit(-1);
  }
  string line;
  // Get rid of 1st row, i.e headers.
  getline(fp, line);
  string title, peptide;
  int i = 0;
  while (getline(fp, line)) {
    int pos = line.find('\t');
    title = line.substr(0, pos);
    peptide = line.substr(pos+1);
    map_title_to_peptide[title] = peptide;
    map_peptide_to_titles[peptide].emplace_back(title);
  }
  fp.close();
  return;
}

void readCluster(string fileName, unordered_map<int, vector<string>>& map_cluster_to_titles, unordered_map<string, int>& map_title_to_cluster) {
  ifstream fp(fileName);
  if (!fp.is_open()) {
    cout << "Error! Can not open" << fileName << endl;
    exit(-1);
  }
  string line;
  // Get rid of 1st row[ID Titles], i.e headers.
  getline(fp, line);
  istringstream iss;
  string title;
  int id;
  while (getline(fp, line)) {
    int pos = line.find('\t');
    id = stoi(line.substr(0, pos));
    iss.clear();
    iss.str(line.substr(pos));
    while (getline(iss, title, ';')) {
      map_cluster_to_titles[id].emplace_back(title);
      map_title_to_cluster[title] = id;
    }
  }
}

int main (int argc, char *argv[]) {
  cli::Parser parser(argc, argv);
  configure_parser(parser);
  parser.run_and_exit_if_error();

  Core::HyperParams params;
  vector<string> files;
  SetParameters(parser, &params, &files);
  cout << params << endl;

  // Read MGF files.
  auto start_time = chrono::high_resolution_clock::now();

  vector<Spectrum*> spectra;
  unordered_map<int, vector<int>> map_spectra_by_charge;
  unordered_map<int, vector<int>> map_spectra_by_precursor_mass;
  unordered_map<string, int> map_ms_titles_to_index;
  unordered_map<string, vector<int>> map_peptide_to_indices;


  int curr_file_spectra_size = 0;
  cout << "Loading spectra..." << endl;
  for (int i = 0; i < files.size(); ++i) {
    cout << "Loading spectra from file: [" << i+1 <<  "/" << files.size() << "]" << endl;
    cout << "file name: " << files[i] << endl;
    IO::ReadSpectraFromMGF(
        &spectra, 
        &map_spectra_by_charge,
        &map_spectra_by_precursor_mass,
        &map_peptide_to_indices,
        &map_ms_titles_to_index,
        &curr_file_spectra_size, 
        files[i], 
        params.mz_scale,
        params.min_mz, 
        params.max_mz, 
        params.precision, 
        params.select_topk,
        params.window_mz, 
        static_cast<PEAK_INTENSITY_RESCALE_METHOD>(params.peak_intensity_rescale_method),
        params.remove_precursor,
        true,
        true,
        false);
    cout << "curr file has " << curr_file_spectra_size << " spectra." << endl;
  }
  cout << endl;

  // Doing experiments from here.
  //auto peptide = "MREDYDSVEQDGDEPGPQR";
  //auto title="OrbiElite04783, Cmpd 15885, +MSn(1111.9666), 44.62 min";
  //auto peptide_ids = map_peptide_to_indices[peptide];
  //cout << "peptide ids length: " <<  peptide_ids.size() << endl;
  //auto peptide_id = peptide_ids[0];
  //int title_id = map_ms_titles_to_index[title];
  //cout << title << ": " << title_id << endl;
  //cout << peptide << ": " << peptide_id << endl;
  //Spectrum first = *spectra[title_id], second = *spectra[peptide_id];
  //cout << "cosine sim: " << 1 - Distance::cosine(first, second, params.precision) << endl;
  //cout << first << endl;
  //cout << second << endl;
  
  Spectrum first, second;
  float dist, sim;
  for (int i = 0; i < spectra.size(); ++i) {
    first = *(spectra[i]);
    // cout << "[i]: " << i << " , " << first << endl;  
    cout << "[i]: " << i << " , pepmass: " << first._precursor_mz << endl;
    cout << "title: " << first._title << " , peptide: " << first._peptide_raw << endl;  
    //if (i == 0 || i == 1) cout << first << endl;
    
    for (int j = i + 1; j < spectra.size(); ++j) {
      second = *(spectra[j]);
      //dist = Distance::cosine(first._filtered_peaks, 
      //                        second._filtered_peaks, params.precision);
      dist = Distance::cosine(first, second, params.precision);
      //cout << "+-+-+-+-" << endl;
      sim = 1 - dist;
      cout << "i: " << i << ", j: " << j << ", sim: " << sim << endl; 
      //cout << first._title << " vs " << second._title << ", sim: " << sim << endl;
    }
    cout << "+-+-+-+-" << endl;
  }
  
  // Releasing memory of spectra.
  start_time = chrono::high_resolution_clock::now();
  cout << "Releasing memory of spectra." << endl;
  for (int i = 0; i < spectra.size(); ++i) {
    delete spectra[i];
  }

  auto end_time = chrono::high_resolution_clock::now();
  auto elapsed_read = chrono::duration_cast<std::chrono::duration<double>>(
      end_time - start_time).count();
  cout << "Releasing memory takes: " << elapsed_read << endl;
  cout << endl;

  return 0;
}

