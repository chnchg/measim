/*
  Module multiplexer to initialize the MEAsim system
 */

#ifndef INITIA_HH
#define INITIA_HH
#include "measim.hh"
#include <prm/argu.hh>
class InitModule // initialization modules
{
	std::string id;
	friend class Initia;
protected:
	MT19937 rng;
	prm::Param param;
	size_t seed;
public:
	InitModule(std::string const & id);
	virtual ~InitModule() {}
	virtual void initialize(MEAsim * sys) = 0;
};
class Initia
{
	std::vector<InitModule *> mlist;
	size_t active;
	prm::Argu ag;
public:
	Initia();
	~Initia();
	void set_module(std::string const & name); // make the active module
	std::string const & get_module() const; // return active module name
	bool alter();
	void initialize(MEAsim * sys);
	std::vector<std::string> get_names() const; // return name of all modules
	prm::Param const & get_param() const;
	prm::Param & get_param();
	static bool call_module_line(int, const std::string & mod_line, void * data);
};
#endif
