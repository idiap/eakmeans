/*
Copyright (c) 2015 Idiap Research Institute, http://www.idiap.ch/
Written by James Newling <jnewling@idiap.ch>

eakmeans is a library for exact and approximate k-means written in C++ and Python. This file is part of eakmeans. eakmeans is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License version 3 as published by the Free Software Foundation. eakmeans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with eakmeans. If not, see <http://www.gnu.org/licenses/>.


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

