/*
Copyright (c) 2015-2018 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <james.newling@gmail.com>
All rights reserved.

eakmeans is a library for exact and approximate k-means written in C++ and
Python. This file is part of eakmeans. See file COPYING for more details.

This file is part of eakmeans.

eakmeans is free software: you can redistribute it and/or modify
it under the terms of the 3-Clause BSD Licence. See
https://opensource.org/licenses/BSD-3-Clause for more details.

eakmeans is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See file
COPYING for more details.
*/

#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "optionsutil.h"
#include "pllkmeansfuncs_void.h"



//TODO: I should auto detect when a sparse algorithm is applied to dense data or vv, would make debugging a bit easier. 

optionsutil::Options getoptions(std::vector<std::string> & algorithms, std::vector<std::string> & initoptions){
	optionsutil::Options uopts;
	
	std::string algstring = algorithms[0];
	for (unsigned i = 1; i < algorithms.size(); ++i){
		algstring = algstring + ", " + algorithms[i];
	}
	
	std::string initstring = initoptions[0];
	for (unsigned i = 1; i < initoptions.size(); ++i){
		initstring = initstring + ", " + initoptions[i];
	}
	
	uopts.add("help", "h", "Display this message ", "s", "");
	
	uopts.add("algorithm", "alg",  
	"The algorithm which may be used, listed above. The complete list of algorithms (under development) includes ........... " + algstring, "s", "expNS");
	
	 //In practice, for exact k-means we suggest one of {expNS, syinNS, selkNS} for low, medium and high dimensions respectively. The complete set of options is " + algstring, "s", "expNS");

	uopts.add("nruns", "nr", "The number of runs (freshly seeded restarts) ", "s", "1");

	uopts.add("ncentroids", "nc", "How many clusters", "i", "");
	
	uopts.add("nthreads", "nth", "How many threads to use. If multiple runs, runs are still run serially but each run is internally parallelised", "i", "1");

	uopts.add("fprecision", "fpr", "Floating point precision, options are 32 (float) and 64 (double). Note that results with 32 and 64 may differ due to rounding in the 32-bit case, for reliable comparisons 64 should be used ", "i", "64");	
	
	uopts.add("cout_verbosity", "cver", "verbosity to terminal, options are 0 (silent) 1 (moderate) 2 (verbose) ", "i", "1");	
	
	uopts.add("file_verbosity", "fver", "verbosity to file, options are 0 (no writing) 2 (write round-by-round summary to soutfn) 3 (write verbose round-by-round data to voutfn) ", "i", "0");	
	
	uopts.add("datainfn", "din", "Path of input data, may contain a first line which is `nrows ncols', subsequent lines should be data, see examples directory. For dense data, each line is whitespace separated values, for sparse data each line should be `l ind:val ... ind:val, see examples directory ", "s", "");	

	
	uopts.add("coutfn", "cou", "Path where final centroids should be written", "s", "");	
	
	uopts.add("loutfn", "lou", "Path where final labels should be written", "s", "");	

	uopts.add("ioutfn", "iou", "Path where initialisation indices should be written", "s", "");
	
	uopts.add("soutfn", "sou", "Absolute path where round-by-round summaries should be written, as well as other bits and bobs. To be set only if file_verbosity \\in {2,3} ", "s", "");	
	
	uopts.add("moutfn", "mou", "Absolute path where multi-run summary should be written", "s", "");	

	uopts.add("voutfn", "vou", "Absolute path where round-by-round verbose output should be written. Useful for studying k-means dynamics, but may increase runtime. Should be be set iff file_verbosity \\in {3}", "s", "");	

	uopts.add("cinfn", "cin", "Absolute path from where initial centroids should be read, each line should be a centroid, values separated by whitespace", "s", "");	
	
	uopts.add("seed", "see", "Initialisation seed", "i", "time(NULL)");	

	
	uopts.add("ind0fn", "iin", "Absolute path from where indices of data for initialising centroids should be read. Each line of file should be an integer in [0, ndata) ", "s", "");	
	
	uopts.add("init0", "ini", std::string("How to initialise the data, currently should be one of {") + initstring  +"}. For minibatch and turbobatch algorithms, the initialisation is performed on the starting batch ", "s", "uniform");	

	uopts.add("maxiter", "mi", "The maximum number of iterations before terminating ", "i", "1e9");	
		
	uopts.add("maxtime", "mt", "Maximum time (in seconds) before terminating ", "f", "1e9");
	
	uopts.add("moutdir", "mdir", "Absolute path of dir where sout-like output should be written for multi-runs, filenames will be run0.txt...run_n.txt ", "s", "");	
	
	uopts.add("valinfn", "vin", "Absolute path to file containing validation data, that is data on which an mse will be computed but which will not be used in learning centroids. Same format as for datainfn. Note that computation of validation mse uses no accelerating algorithm, and can be quite slow. See also valperiod.", "s", "");
	
	uopts.add("valperiod", "vp", "Period at which mse on validation set is to be computed: mse computed in rounds 0 + i*valperiod.", "i", "5");
	
	uopts.add("minibatchsize", "mbs", "The amount of data per batch (per centroid update) using the minibatch algorithm, or the amount of data in the initial batch for growbatch-rho and turbobatch-rho algorithms ", "i", "125");
	
	uopts.add("cmsewritefn", "cmsefn", "The directory to write (clean) mses to. Note that the computation of the mses written here may be time consuming, but this time is not included in the final time tally (this is the only computation which is not included in runtime). Also note that if validation is provided, then the validation is used to compute the cmse, otherwise the training data is", "s", "");
	
	uopts.add("cmserate", "cmr", "The rate at which clean mses will be computed, in milliseconds. Specifically, if more than cmserate milliseconds have passed since the last cmse computation and the opportunity to compute the cmse presents itself, cmse is computed", "i", "1000");
	
	uopts.add("tbrho", "rho", "For growbatch-rho and turbobatch-rho, the doubling threshold for the median ratio (see arXiv:1602.02934  for details)", "s", "1e6");
		
	//uopts.add("gbphi", "gbp", "For turbobatch-rho, the doubling threshold for the median ratio (see arXiv:1602.02934  for details)", "s", "31.42");
	
	uopts.tail = "\nCompulsory options are {datainfn, ncentroids}. Must have at no more than one of {cinfn, ind0fn, init0}, if none will use seed to initialise uniformly. If nruns > 1, then  none of {cinfn, ind0fn, soufn, voutfn, file_verbosity} should be set. If nruns == 1, then moutfn should not be set. Optional independant of all are {algorithm, nthreads, fprecision, verbosity, seed, maxiter, maxtime, coutfn, loutfn, cout_verbosity}. ioutfn can be set iff cinfn is not set. file_verbosity == 0 iff none of {soutfn, voutfn} are set, file_verbosity == 2 iff soutfn set and voutfn not set set, file_verbosity == 3 iff all of {soutfn, voutfn} set ";
	
	
	return uopts;
}


bool algissparse(const std::string & alg){ 
	if (alg.find("sparse") != std::string::npos) {
    return true;
	}
	else{
		return false;
	}
}

bool algisgb(const std::string & alg){
	if (alg.find("gb") != std::string::npos) {
    return true;
	}
	else{
		return false;
	}
}

int main(int argc, char* argv[]){



	std::vector<std::string> onfront {"turbobatch-rho", "exp-ns", "selk-ns", "syin-ns"};
	std::vector<std::string> useralgorithms {"turbobatch-rho", "exp-ns", "selk-ns", "syin-ns", "elk-sn", "syin-sn", "ham", "ann", "minibatch-f", "standard", "selk-sn","elk-ns", "exp-sn", "yin-sn", "minibatch", "growbatch-rho"};
	
	
	std::map<std::string, std::string> mapfromusertoraw;
	
		
	mapfromusertoraw["turbobatch-rho"] =  "gbmse3v1";
	mapfromusertoraw["turbobatch-rho-s"] = "sparsegbmse3v1";
	mapfromusertoraw["growbatch-rho"] =  "gbmse";
	mapfromusertoraw["growbatch-rho-s"] = "sparsegbmsesimple";
	mapfromusertoraw["exp-ns"] =  "expNS" ;
	mapfromusertoraw["selk-ns"] =  "selkNS"; 
	mapfromusertoraw["syin-ns"] =   "syinNS";
	mapfromusertoraw["elk-sn"] =   "elkSN";
	mapfromusertoraw["syin-sn"] =   "syinSN";
	mapfromusertoraw["ham"] =  "ham";
	mapfromusertoraw["ann"] =   "ann";
	mapfromusertoraw["minibatch-f"] =  "minibatch";
	mapfromusertoraw["minibatch-f-s"] = "sparseminibatch";
	mapfromusertoraw["standard"] =  "exactsimplebatch";
	mapfromusertoraw["standard-s"] = "sparsesimple";
	mapfromusertoraw["selk-sn"] =  "selkSN";
	mapfromusertoraw["selk-sn-s"] = "sparsep3v0";
	mapfromusertoraw["elk-ns"] =   "elkNS";
	mapfromusertoraw["exp-sn"] =   "expSN";
	mapfromusertoraw["yin-sn"] =  "p17v6";
	mapfromusertoraw["minibatch"] =   "standardminibatch";
	mapfromusertoraw["minibatch-s"] = "sparsestandardminibatch";


	std::vector<std::string> issparsev {"growbatch-rho", "turbobatch-rho", "selk-sn",  "minibatch-f", "standard", "minibatch"};
	
                                        
	std::map<std::string, std::string> algdefs;
	algdefs["syin-ns"] = "Simplified Yinyang algorithm with ns-bounding, the fastest exact algorithm in dimensions ~ 7->50, memory O(NK) but lighter in memory that selk-ns";
	algdefs["exp-ns"] = "Exponion with ns-bounding, the fastest non-exact algorithm in dimension below ~ 7, memory O(N)";
	algdefs["selk-ns"] = "Simplified version of Elkan with ns-bounding, the fastest exact algorithm in dimensions greater that ~50, memory O(NK)";
	algdefs["turbobatch-rho"] = "Turbobatch algorithm [arXiv:1602.02934], fastest non-exact algorithm, a good choice in general, memory O(NK)";
	algdefs["standard"] = "Lloyd's algorithm, the simple exact algorithm";
	algdefs["selk-sn"] = "Simplified version of Elkan with sn-bounding, in general outperformed by selk-ns";
	algdefs["elk-sn"] = "Elkan, in general outperformed by selk-ns";
	algdefs["syin-sn"] = "Simplified Yinyang algorithm with sn-bounding, in general outperformed by syin-ns";	
	algdefs["ham"] = "Hamerly, in general outperformed by exp-ns";
	algdefs["ann"] = "Annulus, in general outperformed by exp-ns";
	algdefs["minibatch-f"] = "Fixed version of minibatch, in general outperformed by turbobatch-rho";
	algdefs["minibatch"] = "Unfixed version of minibatch, in general outperformed by turbobatch-rho";
	algdefs["exp-sn"] = "Exponion with ns-bounding, in general slightly outperformed by exp-ns";
	algdefs["elk-ns"] = "Elkan with ns-bounding, in general outperformed by selk-ns";
	algdefs["yin-sn"] = "Yinyang algorithm, in general outperformed by syin-ns";
	algdefs["growbatch-rho"] = "Growbatch algorithm [arXiv:1602.02934],  in general outperformed by turbobatch-rho O(N)";

		
	std::string whatisthis("\n\nK-Means clustering\n\nAlgorithms available (\033[1m\033[37mS\033[0mpare / \033[1m\033[37mD\033[0mense)  (\033[31mfastest\033[0m and \033[32mothers\033[0m) : \n" );
		
	//std::string whatisthis("");
	for (auto & k : useralgorithms){
		if (std::find(onfront.begin(), onfront.end(), k) != onfront.end())
		{
			whatisthis = whatisthis + "\033[31m  " +  k + " \033[0m : " + algdefs[k];
		}
		
		else
		{
			whatisthis = whatisthis + "\033[32m  " +  k + " \033[0m : " + algdefs[k];
		}
		
		if (std::find(issparsev.begin(), issparsev.end(), k) != issparsev.end()){
			whatisthis = whatisthis + " [\033[1m\033[37mSD\033[0m] \n";
		}
		
		else{
			whatisthis = whatisthis + " [\033[1m\033[37mD\033[0m] \n";
		}
	}
	
	whatisthis = whatisthis + "\nExample calls: \n-------------\n./bin/withblaskmeans -din examples/sparse_randdim10_train_noheader.txt  -nc 5 --seed 1012 -alg turbobatch-rho-s -cver 2\n./bin/withblaskmeans -nc 5  -din examples/dense_200_5_header.txt -alg exp-ns -cver 1\n-------------" ;
	
	
	
	std::vector<std::string> algorithms {"sta", "expSN", "expNS", "selkSN", "selkNS", "elkSN", "elkNS", "syinNS", "syinSN", "ham", "ann",
		"simple", "simplest", "exactsimplebatch", "p11v0", "p12v6", "p12v7", "p13v0", "p17v3", "p17v6", "p21v3", "p21v4", "p21v5", "p3v0", "p4v2", "p5v1", "p6v0", "minibatch", "standardminibatch", "sparsesimple", "sparseselkSN", "sparsep3v0", "gbsimple", "gbmse", "gbmse3v1", "mb3v0", "elkminibatch", "sparseminibatch", "sparsestandardminibatch", "sparsegbmsesimple", "sparsegbmse3v1"};
			
	std::vector<std::string> initoptions = {"uniform", "kmeans++"};
	
	auto uopts = getoptions(algorithms, initoptions);
	
	std::string algorithm = uopts.options["algorithm"].defval;
	size_t nruns = std::stoi(uopts.options["nruns"].defval);
	size_t nthreads = std::stoi(uopts.options["nthreads"].defval);	
	size_t fprecision = std::stoi(uopts.options["fprecision"].defval);
	int cout_verbosity = std::stoi(uopts.options["cout_verbosity"].defval);
	int file_verbosity = std::stoi(uopts.options["file_verbosity"].defval);
	std::string datainfn = uopts.options["datainfn"].defval;
	std::string coutfn = uopts.options["coutfn"].defval;
	std::string loutfn = uopts.options["loutfn"].defval;	
	std::string ioutfn = uopts.options["ioutfn"].defval;
	std::string soutfn = uopts.options["soutfn"].defval;	
	std::string voutfn = uopts.options["voutfn"].defval;	
	std::string moutfn = uopts.options["moutfn"].defval;	
	std::string moutdir = uopts.options["moutdir"].defval;	
	std::string cinfn = uopts.options["cinfn"].defval;
	std::string cmsewritefn = uopts.options["cmsewritefn"].defval;
	size_t seed = time(NULL);
	std::string ind0fn = uopts.options["ind0fn"].defval;
	std::string init0 = uopts.options["init0"].defval;
	size_t ncentroids = 0;
	size_t maxiter = 1e9;//std::numeric_limits<size_t>::max();
	double maxtime = 1e9;
	std::string valinfn = uopts.options["valinfn"].defval;
	size_t valperiod = std::stoi(uopts.options["valperiod"].defval);
	size_t minibatchsize = std::stoi(uopts.options["minibatchsize"].defval);
	size_t cmserate = std::stoi(uopts.options["cmserate"].defval);
	double tbrho = std::stoi(uopts.options["tbrho"].defval);
	
	std::string name;
	std::string value;
	std::string arg;

	if (argc<3){
		uopts.print(40, 100);
	}
		
	for (int i = 1; i < argc; ++i){		
		arg = argv[i];
		if (arg.compare("-h") == 0 || arg.compare("--help") == 0){
			
			std::cout << whatisthis << std::endl;
			uopts.print(40, 100);
			return 1;
		}
	}
	
	for (int i = 1; i < argc; ++i){	
		arg = argv[i];
		if (i%2 == 1){
			if (arg.substr(0,2).compare("--") == 0){
				name = arg.substr(2);
			}
			
			else if (arg.substr(0,1).compare("-") == 0){
				name = arg.substr(1);
				if (uopts.fullname.count(name)){
					name = uopts.fullname[name];
				}
				
				else {
					std::string erstr = std::string("problem with options : Received flag ") + arg + " which is an invalid shortname. " + name  + " does not correspond to a fullname, ie it is not a key of upots.fullname";
					std::cerr << erstr << std::endl;
					return 1;
				}
			}
			
			else{
				std::cerr << "problem with options : expected a flag (- or --), got a value: " << arg << std::endl;
				return 1;
			}
			
			if (uopts.options.count(name) == 0){
				std::string erstr = std::string("problem with options : invalid fullname : ") + name + ", it is not a key of uopts.options";
				std::cerr << erstr << std::endl;
				return 1;
			} 
			
		}
		else{
			if (arg.substr(0,1).compare("-") == 0){
				std::cerr << "problem with options : expected a value, got a flag  (- or --)" << std::endl;
				return 1;
			}
			
			else{
				

				if (name == "algorithm"){
					if (std::count(algorithms.begin(), algorithms.end(), arg) == 0){
						std::cout << "arg not in raw algorithms list, checking in useralgorithms ... " << std::flush;
						bool sparseversion = false;
						if (arg.substr(arg.size() - 2, 2).compare("-s") == 0){
							std::cout << "alg determined to be sparse ... " << std::flush;
							sparseversion = true;
							arg = arg.substr(0, arg.size() - 2);
						}
						if (std::count(useralgorithms.begin(), useralgorithms.end(), arg) == 0){
							std::cerr << std::string(" problem with options : unrecognised algorithm,  " + arg) << ", please choose one of:\n " << std::endl;
							for (auto & x : useralgorithms){
								std::cerr << x << "  ";
								

								if (std::find(issparsev.begin(), issparsev.end(), x) != issparsev.end()){
									std::cerr << x << "-s";
								}
								
								std::cerr << std::endl;
							}
							for (auto & x : algorithms){
								std::cerr << x << "  ";
							}
							std::cerr << "\n" << std::endl;
							
							
							return 1;
						}
						
						if (sparseversion && (std::find(issparsev.begin(), issparsev.end(), arg) == issparsev.end())){
							std::cerr << "There is no sparse version of " << arg << std::endl;
							return 1;
						}
						
						//arg = mapfromusertoraw[arg];
						if  (sparseversion == true){
							arg = arg + "-s";
						}
						std::cout << "found.\n Mapping " << arg << " to ... " <<  std::flush;
						
						std::cout << mapfromusertoraw[arg] << std::endl;
						arg = mapfromusertoraw[arg];

					}
					algorithm = arg;
					uopts.options["algorithm"].isset = true;
					std::cout << "algorithm set to " << algorithm << std::endl;
				}
				
				else if (name == "nruns"){
					nruns = std::stoi(arg);
					if (nruns > 100){
						std::cerr << std::string("problem with options : runs set to ") + std::to_string(nthreads) + " which exceeds the limit of 100" << std::endl;
						return 1;
					}
					
					else if (nruns < 1){
						std::cerr << std::string("problem with options : nruns set to ") + std::to_string(nruns) + ", it should be greater than equal to 1" << std::endl;
						return 1;
					}
					
					else{
						uopts.options["nruns"].isset = true;
						std::cout << "nruns set to " << nruns << std::endl;
					}
				}
				
				else if (name == "nthreads"){
					nthreads = std::stoi(arg);
					if (nthreads >= 20){
						std::cerr << std::string("problem with options : nthreads set to ") + std::to_string(nthreads) + " which exceeds the limit of 20" << std::endl;
						return 1;
					}
					
					else if (nthreads < 1){
						std::cerr << std::string("problem with options : nthreads set to ") + std::to_string(nthreads) + ", it should be greater than equal to 1" << std::endl;
						return 1;
					}
					
					else{
						uopts.options["nthreads"].isset = true;
						std::cout << "nthreads set to " << nthreads << std::endl;
					}
				}
				
				else if (name == "fprecision"){
					fprecision = std::stoi(arg);
					if ((fprecision == 32) || (fprecision == 64)){
						std::cout << "fprecision set to " << fprecision << std::endl;
						uopts.options["fprecision"].isset = true;
					}
					else{
						std::cerr << "problem with options : fprecision should be 32 or 64, not " << fprecision << std::endl;
						return 1;
					}
				}
				
				
				else if (name == "file_verbosity"){
					file_verbosity = std::stoi(arg);
					if (file_verbosity == 0 || file_verbosity == 2 || file_verbosity == 3){
						std::cout << "file_verbosity set to " << file_verbosity << std::endl;
						uopts.options["file_verbosity"].isset = true;
					}
					else{
						std::cerr << "problem with options : file_verbosity should be 0,2 or 3, not " << file_verbosity  << std::endl;
						return 1;
					}
				}
				
				
				else if (name == "cout_verbosity"){
					cout_verbosity = std::stoi(arg);
					if (cout_verbosity == 0 || cout_verbosity == 1 || cout_verbosity == 2){
						std::cout << "cout_verbosity set to " << cout_verbosity << std::endl;
						uopts.options["cout_verbosity"].isset = true;
					}
					else{
						std::cerr << "problem with options : cout_verbosity should be 0,1 or 2, not " << cout_verbosity  << std::endl;
						return 1;
					}
				}
				
				else if (name == "seed"){
					seed = std::stoi(arg);
					std::cout << "seed set to " << seed << std::endl;
					uopts.options["seed"].isset = true;
				}
				
				else if (name == "ncentroids"){
					
					
					ncentroids = std::stoi(arg);
					
					
					std::cout << "ncentroids set to " << ncentroids << std::endl;
					uopts.options["ncentroids"].isset = true;
				}
				
				
				else if (name == "init0"){
					
					if (std::count(initoptions.begin(), initoptions.end(), arg) == 0){
						std::cerr << std::string("problem with options :  unrecoginsed initialisation option " + arg) << std::endl;
						return 1;
					}
					
					else{
						init0 = arg;
						std::cout << "init0 set to " << arg << std::endl;
						uopts.options["init0"].isset = true;
					}
				}
				
				else if (name == "maxiter"){
					maxiter = std::stoi(arg);
					std::cout << "maxiter set to " << maxiter << std::endl;
					uopts.options["maxiter"].isset = true;
				}

				else if (name == "maxtime"){
					maxtime = std::stof(arg);
					std::cout << "maxtime set to " << maxtime << std::endl;
					uopts.options["maxtime"].isset = true;
				}

				else if (name == "cmserate"){
					cmserate = std::stoi(arg);	
					std::cout << "cmserate set to " << cmserate << std::endl;
					uopts.options["cmserate"].isset = true;
				}
				
				else if (name == "tbrho"){
					std::cout << "setting tbrho ... " << std::flush;
					tbrho = std::stof(arg);
					std::cout << "tbrho set to " << tbrho << std::endl;
					uopts.options["tbrho"].isset = true;
				}
				
				else if (name == "valperiod"){
					valperiod = std::stoi(arg);	
					std::cout << "valperiod set to " << valperiod << std::endl;
					uopts.options["valperiod"].isset = true;
				}
				
				else if (name == "minibatchsize"){
					minibatchsize = std::stoi(arg);	
					std::cout << "minibatchsize set to " << minibatchsize << std::endl;
					uopts.options["minibatchsize"].isset = true;
				}


				
				else if (name == "datainfn" || name == "coutfn" || name == "loutfn" || name == "cinfn" || name == "ind0fn" || name == "ioutfn" || name == "soutfn" || name == "voutfn" || name == "moutfn" || name == "moutdir" || name == "valinfn" || name == "cmsewritefn"){
					
					if ("datainfn" == name){
						datainfn = arg;
						std::cout << "datainfn set to " << arg << std::endl;
						uopts.options["datainfn"].isset = true;
					}
					
					
					else if ("cmsewritefn" == name){
						cmsewritefn = arg;
						std::cout << "cmsewritefn set to " << arg << std::endl;
						uopts.options["cmsewritefn"].isset = true;
					}
					
					
					else if ("coutfn" == name){
						coutfn = arg;
						std::cout << "coutfn set to " << arg << std::endl;
						uopts.options["coutfn"].isset = true;
					}
					
					else if ("loutfn" == name){
						loutfn = arg;
						std::cout << "loutfn set to " << arg << std::endl;
						uopts.options["loutfn"].isset = true;
					}

					else if ("ioutfn" == name){
						ioutfn = arg;
						std::cout << "ioutfn set to " << arg << std::endl;
						uopts.options["ioutfn"].isset = true;
					}

					else if ("soutfn" == name){
						soutfn = arg;
						std::cout << "soutfn set to " << arg << std::endl;
						uopts.options["soutfn"].isset = true;
					}	
					
					else if ("voutfn" == name){
						voutfn = arg;
						std::cout << "voutfn set to " << arg << std::endl;
						uopts.options["voutfn"].isset = true;
					}	

					else if ("moutfn" == name){
						moutfn = arg;
						std::cout << "moutfn set to " << arg << std::endl;
						uopts.options["moutfn"].isset = true;
					}	
					
					else if ("moutdir" == name){
						moutdir = arg;
						std::cout << "moutdir set to " << arg << std::endl;
						uopts.options["moutdir"].isset = true;
					}	
					
					
					else if ("cinfn" == name){
						cinfn = arg;
						std::cout << "cinfn set to " << arg << std::endl;
						uopts.options["cinfn"].isset = true;						
					}
					
					else if ("ind0fn" == name){
						ind0fn = arg;
						std::cout << "ind0fn set to " << arg << std::endl;
						uopts.options["ind0fn"].isset = true;
					}
					
					else if ("valinfn" == name){
						valinfn = arg;
						std::cout << "valinfn set to " << arg << std::endl;
						uopts.options["valinfn"].isset = true;
					}
				}				
			}
		}
	}
	
	std::cout << "all flags have been processed, checking compatibility" << std::endl;
	
	if (uopts.options["minibatchsize"].isset == true && 
	((algorithm.compare("minibatch") != 0) && (algorithm.compare("standardminibatch") != 0) 
	&& (algorithm.compare("gbsimple") != 0) && (algorithm.compare("gbmse") != 0)
	&& (algorithm.compare("gbmse3v1") != 0)
	&& (algorithm.compare("mb3v0") != 0) && (algorithm.compare("elkminibatch") != 0)
	&& (algorithm.compare("sparseminibatch") != 0) && (algorithm.compare("sparsestandardminibatch") != 0)
	&& (algorithm.compare("sparsegbmsesimple") != 0)
	&& (algorithm.compare("sparsegbmse3v1") != 0)
	
	))
	
	{
		
		std::cerr << "minibatchsize set to " << minibatchsize << ", but algorithm is not minibatch (or a variant). Either set algorithm to minibatch or  do not set minibatchsize " << std::endl;
		return 1;
	}
	
	if (!(uopts.options["cmserate"].isset) && uopts.options["valinfn"].isset == true){
		std::cout << "cmse not to be performed and valinfn set, I expect the user wishes to compute validation mses but is not interested in runtimes" << std::endl;
		if (valperiod == 0){
			std::cerr << "problem with options : valinfn is set to " << valinfn << "but valperiod is set to 0. This does not make sense, valperiod should be greater than 0" <<  std::endl;
			return 1;
		}
		
		if (!(cout_verbosity > 1) and !(file_verbosity > 0)){
			std::cerr << "problem with options : cmse not requested, thus : valinfn has been set to " << valinfn << " but verbosity is not compatible : cout_verbosity is " << cout_verbosity  << " and file_verbosity is " << file_verbosity << " . It is required that cout_verbosity > 1 or file_verbosirt > 0" << std::endl;
			return 1;
		} 		
	}
	
	
	if (uopts.options["datainfn"].isset == false){
		std::cerr << "problem with options : datainfn is not an optional parameter, it must be set" << std::endl;
		return 1;
	}
	
	if (uopts.options["ncentroids"].isset == false){
		std::cerr << "problem with options : ncentroids is not an optional parameter, it must be set" << std::endl;
		return 1;
	}
	
	if (uopts.options["cinfn"].isset and uopts.options["ioutfn"].isset){
		std::cerr << "problem with options : cannot set both cinfn and ioutfn" << std::endl;
		return 1;
	}
	
	if ((uopts.options["soutfn"].isset and (file_verbosity != 2 and file_verbosity != 3) )){
		std::cerr << "problem with options : soutfn is set, but file_verbosity is set to " << file_verbosity << ", either unset soutfn or set file_verbosity to 2 (or 3, but then include voutfn as well)" << std::endl;
		return 1;
	}
	
	if ((uopts.options["voutfn"].isset and file_verbosity != 3)){
		std::cerr << "problem with options : voutfn is set, but file_verbosity is set to " << file_verbosity << ", either unset voutfn or set file_verbosity to 3" << std::endl;
		return 1;
	}

	if ((uopts.options["moutfn"].isset and nruns == 1)){
		std::cerr << "problem with options : moutfn is set, nruns is 1, either unset moutfn or set nruns to be greater than 1" << std::endl;
		return 1;
	}
	
	if ((uopts.options["moutdir"].isset and nruns == 1)){
		std::cerr << "problem with options : moutdir is set, nruns is 1, either unset moutdir or set nruns to be greater than 1" << std::endl;
		return 1;
	}
	
	if ((uopts.options["voutfn"].isset and nruns != 1)){
		std::cerr << "problem with options : voutfn is set, nruns is " << nruns << ", either unset voutfn or set nruns to 1" << std::endl;
		return 1;
	}
	
	if ((uopts.options["soutfn"].isset and nruns != 1)){
		std::cerr << "problem with options : soutfn is set, nruns is " << nruns << ", this combination is currently not supported, either unset soutfn or set nruns to 1" << std::endl;
		return 1;
	}
	
		
	if ((uopts.options["voutfn"].isset == false and file_verbosity == 3)){
		std::cerr << "problem with options : voutfn is not set, but file_verbosity is set to 3, either lower file_verbosity to 0 or 2, or set voutfn" << std::endl;
		return 1;
	}
	
	if (nruns != 1 and (uopts.options["cinfn"].isset || uopts.options["ind0fn"].isset )){
		std::cerr << "problem with options : nruns is set to " << nruns << ", but either cinfn or ind0fn is set. For nruns > 1, it is currently not possible to initialise from prespecified centroids or indices. Either reduce nruns to 1, or make sure neither cinfn nor indofn are set. With nruns > 1, init0 may be set. " << std::endl;
		return 1;
	}
	
	int ninitflags = uopts.options["cinfn"].isset + uopts.options["ind0fn"].isset + uopts.options["init0"].isset;
	if (ninitflags > 1){
		std::cerr << "at most one 1 of cinfn, ind0fn and init0 can be set" << std::endl;
		return 1;		
	}
	
	
	if (
	(uopts.options["cmsewritefn"].isset and !(uopts.options["cmserate"].isset)) ||
	(!(uopts.options["cmsewritefn"].isset) and uopts.options["cmserate"].isset))
	{
			std::cerr << "problem with options : either 0 or 2 of cmsewritefn and cmserate need to be set\n";
			return 1;
	}
	
	if (!algisgb(algorithm) && uopts.options["tbrho"].isset){
		std::cerr << "problem with options : The algorithm " << algorithm << " is not a gb algorithm (the frag gb does not appear in its name) and yet tbrho has been set\n" ;
		return 1;		
	}
	
	if (ncentroids == 2){
		if (algorithm.compare("yin-sn") == 0 || algorithm.compare("syin-sn") == 0 || algorithm.compare("syin-ns") == 0  ||algorithm.compare("exp-sn") == 0 || algorithm.compare("exp-ns") == 0){
			std::cerr << "Currently no exp or yin algorithms work with 2 centroids, please select a different algorithm \n" ;
			return 1;
		}
	}
	
	if (ncentroids < 1){
			std::cerr << "ncentroids should be greater than 1 \n" ;
			return 1;
	}


	bool issparse = algissparse(algorithm);
	
	
	std::cout << "flags passed appear to be compatible" << std::endl;	


	srand(seed);
	
	

	bool setseed = uopts.options["seed"].isset;
	

	std::cout << "Calling solvewrite (f/d)" << std::endl;
	if (fprecision == 64){
		cluster::solvewrited(
		algorithm, //1
		issparse, //2
		nruns, //3
		nthreads, //4
		cout_verbosity, //5 
		file_verbosity, //6
		datainfn, //7
		coutfn, //8
		loutfn, //9
		ioutfn, //10
		soutfn, //11
		voutfn, //12
		moutfn, //13
		moutdir, //14
		cinfn, //15
		ind0fn, //16
		init0, //17
		setseed, //18
		seed, //19
		ncentroids, //20
		maxiter, //21
		maxtime, //22
		valinfn, //23
		valperiod, //24
		minibatchsize, //25
		cmsewritefn, //26
		cmserate, //27
		static_cast<double>(1./(tbrho + 0.000001)) //because code works with the inverse of tbrho
		);
	}
	
	else{
		cluster::solvewritef(algorithm, issparse, nruns, nthreads, cout_verbosity, file_verbosity, datainfn, coutfn, loutfn, ioutfn, soutfn, voutfn, moutfn, moutdir, cinfn, ind0fn, init0, setseed, seed, ncentroids, maxiter, maxtime, valinfn, valperiod, minibatchsize, cmsewritefn, cmserate, static_cast<float>(1./(tbrho + 0.000001))); //because code works with the inverse of tbrho
	}
}	
	
