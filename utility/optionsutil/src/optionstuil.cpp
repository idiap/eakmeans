/*
EAKMeans is a fast Exact K-means library written in C++ with 
command-line interface, shared library + header files and 
Python bindings

Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

This file is part of EAKMeans.

EAKMeans is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License version 3 as
published by the Free Software Foundation.

EAKMeans is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with EAKMeans. If not, see <http://www.gnu.org/licenses/>.



*/

#include <utility>
#include <string>
#include <map>
#include <iostream>
#include <iomanip>
#include <exception>
#include "optionsutil.h"
namespace optionsutil{

Option::Option(std::string fn, std::string sn, std::string desc, std::string tp, std::string dv): fullname(std::move(fn)), shortname(std::move(sn)), description(std::move(desc)), type(std::move(tp)), defval(std::move(dv)), isset(false) {
	if (defval.compare("") == 0){
		definition = "--" + fullname + " -" + shortname + "  "  + type;
	}
	else{
		definition = "--" + fullname + " -" + shortname + "  "  + type + "  (" + defval + ")  ";
	}
}
		
Option::Option(){
	throw std::logic_error("Default constructor for Option called, this should never happen");
}
		
		
void Option::print(unsigned tab1, unsigned tab2){
	unsigned width = tab2 - tab1;
	unsigned margin = 2;
	
	std::cout << definition;
	if (definition.size() > tab1-margin){
		 std::cout << "\n";
		 std::cout << std::setw(tab1) << " ";
	}
	else{
		std::cout << std::setw(tab1 - definition.size()) << " ";
	}
	
	
	unsigned fragi = 0;
	while(fragi < description.size()/width){
		std::cout << description.substr(fragi*width, width) << "\n";
		std::cout << std::setw(tab1) << " ";
		++fragi;
	}
	std::cout << description.substr(fragi*width) << std::endl;
}

/* To do : is it possible to have the following class using variadic templates:
 * Options<int, float, std::string> options;
 * options.add("name1", anint)
 * options.add("name2", astring)
 * ...
 * ?
 * */
				
void Options::add(Option && o){
	fullname.emplace(std::make_pair(o.shortname, o.fullname));
	options.emplace(std::make_pair(o.fullname, std::move(o)));
}

void Options::add(std::string fn, std::string sn, std::string desc, std::string type, std::string defval){
	fullname.emplace(std::make_pair(sn, fn));
	options.emplace(std::make_pair(fn, Option(fn, std::move(sn), std::move(desc), std::move(type), std::move(defval))));
}

void Options::print(unsigned tab1, unsigned tab2){
	
	std::cout << "\nThe options are, \n";
	
	std::cout << std::endl;
	for (auto & p : options){
		p.second.print(tab1, tab2);
	}
	
	unsigned fragi = 0;
	while(fragi < tail.size()/tab2){
		std::cout << tail.substr(fragi*tab2, tab2) << "\n";
		++fragi;
	}
	std::cout << tail.substr(fragi*tab2) <<"\n" << std::endl;
	
}

}
