#ifndef SPARSEDATASETS_H
#define SPARSEDATASETS_H

#include <vector> 
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

#include "sample.h"

//TODO: put somewhere else:

namespace sparse{
	
template <typename TInt, typename TFloat>
class SparseData{
	public:
		std::vector<TFloat> values;
		std::vector<TInt> indices;
		std::vector<TInt> starts;
		std::vector<std::string> labels;
		
		TInt ndata;
		TInt dimension;
		
		void set_ndata_dimension(){
			//determine ndata and dimension from indices and starts
			ndata = starts.size() - 1;
			dimension = *std::max_element(indices.begin(), indices.end()) + 1;
		}
		
		void print() const{
			std::cout << ndata << "\t " << dimension << std::endl;
			for (TInt i = 0; i < starts.size() - 1; ++i){
				std::cout << i << "\t " << labels[i] << " \t " << starts[i] << " ---> " << starts[i+1] << std::endl;
				for (TInt j = starts[i]; j < starts[i+1]; ++j){
				
					std::cout << j << " ( " << indices[j] << " , " << values[j] << " )   " << std::flush;
				}
				std::cout << "\n--------------" << std::endl;
			}
			std::cout << "\n+       ++       ++       ++       ++       ++       +" << std::endl;
		}
		
		
		void print_compact() const{
			
			std::cout << "[" << ndata << ", " << dimension << "]" << std::endl;
			for (TInt i = 0; i < starts.size() - 1; ++i){

				std::cout << i << "\t";
				for (TInt j = starts[i]; j < starts[i+1]; ++j){
					std::cout << indices[j] << ":" << values[j] << " " << std::flush;
				}
				std::cout << std::endl;
			}
		}
	
			
			
		void write(const std::string & filename, bool dimheader = true){
			std::ofstream file(filename);
			if (dimheader == true){
				file << ndata << "\t" << dimension << "\n";
			}
			for (TInt i =0; i < this->ndata; ++i){
				
				file << labels[i] << " ";
				for (TInt j = starts[i]; j < starts[i + 1]; ++j){
					file << indices[j] << ":" << values[j] << " ";
				}
				file << "\n";
			}
			
			file.close();
		}
		
		void write_dense(const std::string & filename, bool dimheader = true){
			std::ofstream file(filename);
			if (dimheader == true){
				file << ndata << "\t" << dimension << "\n";
			}
			std::vector<TFloat> dvalues (this->dimension);
			for (TInt i =0; i < this->ndata; ++i){
				std::fill_n(dvalues.data(), this->dimension, 0);
				for (TInt j = starts[i]; j < starts[i + 1]; ++j){
					dvalues[indices[j]] = values[j];
				}
				for (TInt j = 0; j < this->dimension; ++j){
					file << dvalues[j] << "\t";
				}
				file << "\n";
			}
			
			file.close();
			
		}
		
		SparseData & operator= ( const SparseData & ) = default;	
		
		SparseData() : values {}, indices {},  starts{0}, labels {}, ndata(0), dimension(0) {
			
		}
		
		SparseData(
		std::vector<TFloat> && values, 
		std::vector<TInt> && indices, 
		std::vector<TInt> && starts, 
		std::vector<std::string> && labels):
		values(values), indices(indices), starts(starts), labels(labels) {
			set_ndata_dimension();
		}   
		
		SparseData(const std::string & filename, bool dimheader = false) : SparseData() {
			from_file(filename, dimheader);
		}
		
		//untested, but should work (to use for dimheader = false)
		SparseData(const std::vector<std::string> & filenames, bool dimheader = false) : SparseData() {
			from_files(filenames, dimheader);
		}
		
		
		//untested, but should work. 
		void from_files(const std::vector<std::string> & filenames, bool dimheader = false){
			for (auto fn : filenames){
				this->from_file(fn, dimheader);
			}
		}
			
		
		void from_file(const std::string & filename, bool dimheader = false){
			
			std::ifstream dfile(filename, std::ios_base::in);		
			std::string aline;
			TInt push_backs=this->starts[this->starts.size() - 1];
			std::string label_i;
			std::string fragment;
			std::string fragsplitter (":");
			
			
			TInt filendata;
			TInt filedim;
			
			short split_position;
			

			if (dfile.is_open()){
				
				if (dimheader == true){ //first line is `ndata dimension'
					std::getline(dfile, aline);
					std::stringstream ss(aline);
					ss >> filendata;
					ss >> filedim;
				}
				
				
				std::vector<std::string> datalines;
				while (!dfile.eof()){
					std::getline(dfile, aline);
					if (aline.compare("") != 0){
						datalines.push_back(aline);
						//std::cout << "-----------> " << aline << std::endl;
					}

				}
				
					
				
				std::cout << "shuffling input data ..." << std::flush;
				
				auto permuted_range = randomutil::sample::get_permuted_range(datalines.size(), rand);
				//shuffle the lines random order. Make vector of strings, one per line? Use random permutation to fill data?
				
				std::cout << "done " << std::endl;
				
				
				for (auto & ri : permuted_range){
					
					//std::cout << ri << " " << std::flush;
					auto line = datalines[ri];
					
					if (line.size() > 1){
						std::stringstream ss(line);
						ss >> label_i;
						while (ss){
							ss >> fragment;
							split_position = fragment.find_first_of(fragsplitter);
							if (ss){ //I don't understand why this if statement is necessary, but it is
				
								this->indices.push_back(std::stoi(fragment.substr(0, split_position)));
								this->values.push_back(std::stof(fragment.substr(split_position+1, fragment.size())));
								++push_backs;
								
							}
						}
						this->labels.push_back(label_i);	
						this->starts.push_back(push_backs);
					}
				}		
			}
			
			else{
				throw std::runtime_error(std::string("file `") + filename + "' not open (in sparsedatasets) ");
			}
			
			this->set_ndata_dimension();
		}
		

};

}

#endif
