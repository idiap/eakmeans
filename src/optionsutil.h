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

