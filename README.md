## Introduction
msSLASH (standing for **m**ass **s**pectral **S**earching using **L**oc**A**lity **S**ensitive **H**ashing) was developed by [Lei Wang](mailto:wang558@indiana.edu), [Sujun Li](https://scholar.google.com/citations?user=y4keCocAAAAJ&hl=en) and [Haixu Tang*](https://www.sice.indiana.edu/all-people/profile.html?profile_id=308), for performing rapid MSMS spectral library searching while preserving sensitivity. Multithreading is enabled in this package.

## Prerequisites
g++ with version 5.1.0+ is required.

## Installation
1. `cd` to the main directory of msSLASH.
2. Type `./install.sh`
3. Two executable files will be placed under *bin* directory: *./msSLASH* for searching with locality sensitive hashing technique and *./bruteforce* for searching using naive method, respectively.  The `bruteforce` method literally compares input spectrum with each candidate of similar precursor and same charge, using the same scoring function as `msSLASH`, thus serving as benchmark

## Usage
- This process will search input spectra dataset against designated target library and the respective decoy library.  A csv file will be generated at the end to hold the searching results.
  1. `cd bin`
  2. Usage: `./msSLASH -e experimental_spectra -d decoy_library -l target_library [-c fragment_precision] [-t threads] [-n hash_func_num] [-i iteration] [-s min_similarity] [-a min_mz] [-r rescale_method] [-o output_file] [-m precursor_mass_tolerance].`
  3. Typical example: `./msSLASH -n 100 -i 10 -t 20 -l $target -e $exp -d $decoy -o $out_file -r 1 -m 0.05 --precision 0.5  --min_mz 0 -u`. You will find the searching result file with name `$out_file` if specified, otherwise you will find it in the same dir as the `$exp`, with suffix `.msSLASH.tsv`.
  4. Description
  Type `./msSLASH -h` to see full list of command options.
        * -u,    --unfragment
            > [Bool] Filter unfragmented ms2 if specified
    This parameter is optional. The default value is '0'.
        
        * -i,    --iteration
            > [Int] iteration for searching with LSH
     This parameter is optional. The default value is '100'.

        * -n,    --hash
            > [Int] Hash functions per hash table.
     This parameter is optional. The default value is '8'.

        * -t,    --threads
            > [Int] Threads to use.
     This parameter is optional. The default value is '20'.

        * -a,    --min_mz
            > [Float] min mz for peaks to consider.
     This parameter is optional. The default value is '200'.

        * -s,    --similarity
            > [Float] Minimum cosine similalrity when reporting a hit.
     This parameter is optional. The default value is '0'.
    
        * -c    --precision
            > [Float] fragment precision in Da.
    This parameter is optional. The default value is '0.500000'. 

        * -m    --precursor_mass_tolerance
            > [Float] precursor mass tolerance in Da.
    This parameter is optional. The default value is '0.050000'. 

        * -r    --rescale
            > [Int] PEAK_INTENSITY_RESCALE_METHOD.
   This parameter is optional. The default value is '1'. 
   
        * -o    --out_file
            > [String] out file contanining msSLASH searching results.
   This parameter is optional. The default value is ''.

        * -d    --decoy  (required)
            > [String] decoy library mgf file.

        * -e    --experimental  (required)
            > [String] experimental mgf file.

        * -l    --library       (required)
            > [String] library mgf file.


- It is very important to select proper paramters for msSLASH to search spectral library with high sensitivity and specificity. Specifically:
  1. `hash_func_num` controls the collision probability of two spectra in a single hash table, the more hash functions, the smaller collision probabiilty, the more buckets, the faster msSLASH program. 7~10 is a good starting point for MassIVE HCD spectral library, NIST CID library, and predicted human proteome library.
  2. `iteration` controls the number of hash tables to use with the aim to increase the probabiilty of two similar spectra to collide in at least one hash table, the more iterations, the higher collision probability, the less clusters, the slower mcCRUSH program. 100 iterations is a good starting point to play with.
  3. `rescale` specifies transformation that will be applied to each peaks.  0 for no transformation, 1 for log transformation, 2 for sqrt transformation.
  4. `precision` specifies the fragment precision in Da.  Preferably set to value equivalently 1Th, i.e. 0.5Da for doubly charged spectra, 0.33 Da for tribly charged spectra


## Citation
TBD

## Questions
Please contact Lei Wang (wang558@indiana.edu) for assistance.
## Acknowledgement
This work was supported by the NIH grant 1R01AI108888 and the Indiana University Precision Health Initiative.

