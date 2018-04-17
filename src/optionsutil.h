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

#include <string>
#include <map>

namespace optionsutil{

class Option{
	public:
		std::string fullname;
		std::string shortname;
		std::string description;
		std::string type;
		std::string defval;
		bool isset;
		
		Option(std::string fn, std::string sn, std::string desc, std::string tp, std::string dv);
		Option();
		void print(unsigned tab1, unsigned tab2);
		
	private:
		std::string definition;
};

class Options{
	public:
		std::map<std::string, Option> options;
		std::map<std::string, std::string> fullname;
		std::string tail;
		void add(Option && o);
		void add(std::string fn, std::string sn, std::string desc, std::string type, std::string defval);
		void print(unsigned tab1 = 40, unsigned tab2 = 85);
};

}

