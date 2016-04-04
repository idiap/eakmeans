/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


*/

#ifndef TXTDATASETS_H
#define TXTDATASETS_H

namespace datasets{

static std::string sparse_data_dir("/idiap/temp/jnewling/data/sparsedata/rcv1rcv/");


class TxtDataset{
	public:
		std::string name;
		unsigned nd;
		unsigned dim;
		std::string datapath_dim;
		std::string datapath_dimless;
		TxtDataset(){}; //why do I need this?
		TxtDataset(const std::string & name, unsigned nd, unsigned dim):name(name), nd(nd), dim(dim){
		 datapath_dim = "/idiap/temp/jnewling/data/txtdata/normalised/" + name + "_" + std::to_string(nd) + "_" + std::to_string(dim) + "_cnormed.txt";
		 datapath_dimless = "/idiap/temp/jnewling/data/txtdata/normalised/" + name + "_" + std::to_string(nd) + "_" + std::to_string(dim) + "_cnormed_dimless.txt";
		}
};

class TrainTestDataset{ //TODO replace TrainTestDataset everywhere with TrainTestDataset. 
	public:
	
		std::string name;
				
		std::string datapath_train_dim;
		std::string datapath_train_dimless;
		std::string datapath_test_dim;
		std::string datapath_test_dimless;
		
		TrainTestDataset(){};
		
		//, nd(nd), dim(dim) , unsigned nd, unsigned dim
		//, const std::string & datapath_dim, const std::string & datapath_dimless
		//  unsigned nd;
		//  unsigned dim;

		TrainTestDataset(const std::string & name, 
		std::string rootdir = "/idiap/temp/jnewling/data/sparsedata/trainandtest/"
		): name(name){
			datapath_train_dim = rootdir + name + "_train_withdims.txt";
			datapath_train_dimless = rootdir + name + "_train_dimless.txt";
 			datapath_test_dim = rootdir + name + "_test_withdims.txt";
			datapath_test_dimless = rootdir + name + "_test_dimless.txt";
		}	
};

std::vector<TrainTestDataset> sparse_datasets{
	
	{"truercv"},
	{"truercvos"},
	{"rcv"}, //558700 , 0, "None", "/idiap/temp/jnewling/data/sparsedata/rcv1rcv/all_shuffled.txt"},
	{"nips"}, //1500 , 0, "None", "/idiap/temp/jnewling/data/sparsedata/bagofwords/nips.txt"},
	{"nytimes"},// 299751 , 102661, "None", "/idiap/temp/jnewling/data/sparsedata/bagofwords/nytimes.txt"}
	{"randdim5"},
	{"randdim6"},
	{"infimnist", "/idiap/temp/jnewling/data/densedata/trainandtest/"},
	{"infimnist28by28", "/idiap/temp/jnewling/data/densedata/trainandtest/"}
	
};

std::vector<TxtDataset> txt_datasets {
	
	{"tsn", 200000, 4},
	{"conflongdemo", 164860, 3},
	{"skinseg", 200000, 4},
	{"wcomp", 165630, 15}, 
	{"kegg", 65550, 28},
	
	{"miniboone", 130060, 50},
	{"covtype", 581012, 55},
	{"gassensor", 13910, 128},
	{"uscensus", 2458285, 68}, 
	{"colormoments", 68040,9},
	
	{"ldfpads", 164850, 3},
	{"kddcup98", 95000, 310},
	{"kddcup04bio", 145750, 74},
	{"egeod", 5580, 31099},
	{"mnist50", 60000, 50},

	{"house16H", 22780, 17},
	{"mv", 40760, 11},
	{"europe", 169300, 2},
	{"birch3", 100000, 2},
	{"mnist", 60000, 784},

	{"stl10", 1000000, 108},
	{"random", 1000000, 30},
	{"random2", 1000000, 2},
	{"small", 1000, 2},

};



}


#endif
