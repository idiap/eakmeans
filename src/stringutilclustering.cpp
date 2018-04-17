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

#include "stringutilclustering.h"
#include <string>
#include <stdexcept>
#include <cmath>

namespace stringutil{
namespace clustering{
	
namespace helper{
std::string getstars(){
	 return "*****************************";
}
}


		
namespace pll{
namespace exact{

//startsummary verbosity 1
std::string getstartsummary_v1(std::string algname, size_t memory_usage, float mse, float val_mse){
	std::string starts = std::string("\n\n") + 
	stringutil::clustering::helper::getstars() + 
	"\nkmeans with algorithm : " + algname + 
	"\t estimated memory required : " + std::to_string(memory_usage) + " bytes " + 
	"\t initial mse : " + std::to_string(mse) +
	"\t initial validation mse : ";
	
	if (val_mse >= 0){ 
		starts = starts + std::to_string(val_mse) + "\n";
	}
	
	else{
		starts = starts + "N/A\n";
	}
	
	return starts;
}

//startsummary verbosity 2
std::string getstartsummary_v2(std::string algname, size_t memory_usage, float mse, float val_mse){
	return getstartsummary_v1(algname, memory_usage, mse, val_mse) + 
	"round\troundchanges\tnXcalcs\tncalcs\ttime\tmse\tvalmse\n";
}



//round summary verbosity 1
std::string getroundsummary_v1(size_t roundchanges){
	return std::to_string(roundchanges) + " ";
}
//round summary verbosity 2
std::string getroundsummary_v2(size_t round, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, float val_mse){
	std::string summary =  std::to_string(round) + "\t" + std::to_string(roundchanges) + "\t" + std::to_string(ncalcs_X) + "\t" + std::to_string(ncalcs) + "\t" + std::to_string(size_t(t_total)) + "\t" + std::to_string(mse);
	if (val_mse > 0){
		summary = summary + "\t" +  std::to_string(val_mse)  + "\n";
	}
	else{
		summary = summary + "\t-\n";
	}
	return summary;
}





//finalsummary verbosity 1
std::string getfinalsummary_v1(size_t round, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total,  float mse, size_t n_empty_clusters, float val_mse){
	std::string finalstring = std::string("\nrounds: ") + std::to_string(round) + std::string("\tround changes : ") + std::to_string(roundchanges)  + "\tdistance calcs in X-update : " + std::to_string(ncalcs_X) + "\ttotal distance calcs : " + std::to_string(ncalcs) +  "\t time : " + std::to_string(size_t(t_total)) + " ms\n"  +
	"final mse : " + std::to_string(mse) + "\t # empty : " 
	+ std::to_string(n_empty_clusters);
	
	if (val_mse > 0){
		finalstring = finalstring + "\t final validation mse : " + std::to_string(val_mse);
	}
	finalstring += "\n" + stringutil::clustering::helper::getstars() + "\n" ;
	
	return finalstring;
}

//finalsummary verbosity 2
std::string getfinalsummary_v2(size_t round, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, size_t n_empty_clusters, float val_mse){
	return getfinalsummary_v1(round, roundchanges, ncalcs_X, ncalcs, t_total, mse, n_empty_clusters, val_mse);

}
}

namespace minibatch{


std::string getstartsummary_v1(std::string algname, size_t memory_usage, float val_mse){
	return std::string("\n\n") + 
	stringutil::clustering::helper::getstars() + 
	"\nkmeans with algorithm : " + algname + 
	"\t estimated memory required : " + std::to_string(memory_usage) + " bytes " + 
	"\t init validation mse : " + std::to_string(val_mse) + 	
	"\n";
}

//startsummary verbosity 2
std::string getstartsummary_v2(std::string algname, size_t memory_usage, float val_mse){
	return getstartsummary_v1(algname, memory_usage, val_mse) + 
	"round \tsubround \troundchanges \tnXcalcs \tncalcs \ttime \tmse \tvalmse \n";
}


//round summary verbosity 1
std::string getroundsummary_v1(size_t subround, size_t nsubrounds, size_t roundchanges){
	if (subround == nsubrounds - 1) { return std::to_string(roundchanges) + " "; }
	else { return "."; }
}

std::string getroundsummary_v2(size_t round, size_t nsubrounds, size_t subround, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, float val_mse){
	
	std::string rcs = "-";
	if (subround == nsubrounds - 1){
		rcs = std::to_string(roundchanges);
	}
	
	std::string summary =  std::to_string(round) + "\t" + std::to_string(subround)  + "\t\t" + rcs + "\t\t" + std::to_string(ncalcs_X) + "\t\t" + std::to_string(ncalcs) + "\t" + std::to_string(size_t(t_total));

	if (mse > 0){
		summary = summary +  "\t" + std::to_string(mse);
	}
	else if (std::fabs(mse + 1) < 0.01) {
		summary = summary +  "\t" + "-";
	}
	
	else{
		throw std::runtime_error("mse value of " + std::to_string(mse) + " is unexpected, caught in minibatch::getroundsummary_v2, bailing ");
	}

	

	if (val_mse > 0){
		summary = summary + "\t" +  std::to_string(val_mse) + "\n";
	}
	else{
		summary = summary + "\t-\n";
	}
	return summary;
}


std::string getfinalsummary_v1(size_t round, size_t nsubrounds, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, size_t n_empty_clusters, float val_mse){
	std::string finalstring = std::string("\nTotal rounds: ") + std::to_string(round/nsubrounds) + std::string("\tMini rounds : ") + std::to_string(round) + std::string("\tround changes : ") + std::to_string(roundchanges)  + "\tdistance calcs in X-update : " + std::to_string(ncalcs_X) + "\ttotal distance calcs : " + std::to_string(ncalcs) +  "\t time : " + std::to_string(size_t(t_total)) + " ms\n"  +
	"final mse : " + std::to_string(mse) + "\t # empty : " 
	+ std::to_string(n_empty_clusters);
	
	if (val_mse > 0){
		finalstring = finalstring + "\t final validation mse : " + std::to_string(val_mse);
	}
	finalstring += "\n" + stringutil::clustering::helper::getstars() + "\n" ;
	
	return finalstring;
}


std::string getfinalsummary_v2(size_t round, size_t nsubrounds, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, size_t n_empty_clusters, float val_mse){
	return getfinalsummary_v1(round, nsubrounds, roundchanges, ncalcs_X, ncalcs, t_total, mse, n_empty_clusters, val_mse);
}

	
	
}


namespace growbatch{


std::string getstartsummary_v1(std::string algname, size_t memory_usage, float val_mse){
	return std::string("\n\n") + 
	stringutil::clustering::helper::getstars() + 
	"\nkmeans with algorithm : " + algname + 
	"\t estimated memory required : " + std::to_string(memory_usage) + " bytes " + 
	"\t init validation mse : " + std::to_string(val_mse) + 	
	"\n";
}

//startsummary verbosity 2
std::string getstartsummary_v2(std::string algname, size_t memory_usage, float val_mse){
	return getstartsummary_v1(algname, memory_usage, val_mse) + 
	"round \tnactive \td_C/d_AB \troundchanges \tnXcalcs \tncalcs \ttime \tmse \tvalmse \n";
}


//round summary verbosity 1
std::string getroundsummary_v1(size_t roundchanges, bool didgrow){
	std::string summary = std::to_string(roundchanges) + " ";
	if (didgrow == true){
		summary = summary + "* ";
	}
	return 	summary;
}

std::string getroundsummary_v2(size_t round, size_t nactive, float d_C__by__d_AB, size_t roundchanges, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, float val_mse){
	
	
	std::string summary =  std::to_string(round) + "\t" + std::to_string(nactive)  + "\t" + std::to_string(d_C__by__d_AB) + "\t" + std::to_string(roundchanges) + "\t" + std::to_string(ncalcs_X) + "\t" + std::to_string(ncalcs) + "\t" + std::to_string(size_t(t_total));

	if (mse > 0){
		summary = summary +  "\t" + std::to_string(mse);
	}
	else{
		summary = summary +  "\t" + "-";
	}

	

	if (val_mse > 0){
		summary = summary + "\t" +  std::to_string(val_mse) + "\n";
	}
	else{
		summary = summary + "\t-\n";
	}
	return summary;
}


std::string getfinalsummary_v1(size_t round, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, size_t n_empty_clusters, float val_mse){
	std::string finalstring = std::string("\nTotal rounds: ") + std::to_string(round) + "\tdistance calcs in X-update : " + std::to_string(ncalcs_X) + "\ttotal distance calcs : " + std::to_string(ncalcs) +  "\t time : " + std::to_string(size_t(t_total)) + " ms\n"  +
	"final mse : " + std::to_string(mse) + "\t # empty : " 
	+ std::to_string(n_empty_clusters);
	
	if (val_mse > 0){
		finalstring = finalstring + "\t final validation mse : " + std::to_string(val_mse);
	}
	finalstring += "\n" + stringutil::clustering::helper::getstars() + "\n" ;
	
	return finalstring;
}


std::string getfinalsummary_v2(size_t round, size_t ncalcs_X, size_t ncalcs, size_t t_total, float mse, size_t n_empty_clusters, float val_mse){
	return getfinalsummary_v1(round, ncalcs_X, ncalcs, t_total, mse, n_empty_clusters, val_mse);
}

	
	
}



}

}
}
