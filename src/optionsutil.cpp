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
		 std::cout << " \n";
		 std::cout << std::setw(tab1) << " ";
	}
	else{
		std::cout << std::setw(tab1 - definition.size()) << " ";
	}
	
	
	unsigned fragi = 0;
	//unsigned currenti = 0;
	//std::string nextline("");
	//while (currenti < description.size()){
		//nextline = description.substr(currenti, width);
		//if (nextline.find("\0") != std::string::npos) {
		
		//}
	//}
	
	while(fragi < description.size()/width){
		
		//\033
		
		//if (description.substr(fragi*width, width).find("\\0") != std::string::npos){
			//std::cout << "-------------------------------------------------------------" << std::endl;
		//}
		std::cout << description.substr(fragi*width, width) << " \n";
		std::cout << std::setw(tab1) << " ";
		++fragi;
	}
	std::cout << description.substr(fragi*width) << " \n" << std::endl;
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
	
	//std::cout << "\n-------------------------------------------------\n";
	std::cout << "\nThe options are of the form,\n--full_option_name  -abridged_name  type  (default) \n\n";
	
	
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
