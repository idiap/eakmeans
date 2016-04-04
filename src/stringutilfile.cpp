#include "stringutilfile.h"


#include <fstream>
#include <sstream>
#include <stdexcept>

#include "stringutilbase.h"


#include <cstdlib>

#include <iostream>

namespace stringutilfile{
	

//stolen from http://stackoverflow.com/questions/2844817/how-do-i-check-if-a-c-string-is-an-int
inline bool is_integer(const std::string & s){
	if (s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+'))){
		return false;
	}
	
	char * p;
	
	strtol(s.c_str(), &p, 10);
	
	return (*p == 0);
}
	
bool file_has_2int_header(const std::string & filename){
	std::ifstream dfile(filename, std::ios_base::in);		
	std::string line;
	if (!dfile.is_open()){
		throw std::runtime_error(std::string("The file ") + filename + " probably does not exist. Cannot determine whether the file has the 2 int header or not, as the file does not seem to exist." );
	}
	std::getline(dfile, line);
	auto bob = stringutil::split(line);
	
	/* first determine that it contains 2 nuggets: */
	if (bob.size() != 2){
		//"file does not have 2 frags, it has :  " << bob.size() << " frags " << std::endl;
		return false;
	}
	/* next test that they are indeed integers */	
	if (is_integer(bob[0]) and is_integer(bob[1])){
		return true;
	}

	//"fail due to a non int in header :  " <<  is_integer(bob[0]) << " " << is_integer(bob[1])  << std::endl;
	//"fail due to a non int in header :  |" <<  bob[0] << "|  |" << bob[1]  << "|" <<  std::endl;
	
	return false;
}

}
