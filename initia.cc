#include "initia.hh"
#include "mod_init_modules.hh"
#include <cstdlib>

InitModule::InitModule(std::string const & id) :
	id(id)
{
	using namespace prm::tag;
	param.add_var("seed", seed = 0)
		<< Name("seed")
		<< Desc("random seed")
		<< Range(0, 100000)
		<< Step(1, 100)
		<< CmdLine()
		<< Save();
}

Initia::Initia() :
	active(0)
{
	mlist.push_back(new mod_init::MoritaSong);
	mlist.push_back(new mod_init::Minimal);
	mlist.push_back(new mod_init::Random);
	mlist.push_back(new mod_init::SongSmallWorld);
	mlist.push_back(new mod_init::Lattice);
	mlist.push_back(new mod_init::AdjFile);
	mlist.push_back(new mod_init::File);
	mlist.push_back(new mod_init::Rewire);
	mlist.push_back(new mod_init::NK);
	mlist.push_back(new mod_init::AxonSearch);
	mlist.push_back(new mod_init::DistaGau);
	mlist.push_back(new mod_init::DistExp);
	mlist.push_back(new mod_init::DistConst);
	mlist.push_back(new mod_init::Shordist);
	mlist.push_back(new mod_init::Distexpin);
	mlist.push_back(new mod_init::SelfConn);
	mlist.push_back(new mod_init::SelfMoritaSong);
}

Initia::~Initia()
{
	for (std::vector<InitModule *>::iterator i = mlist.begin(); i != mlist.end(); i ++) {
		delete * i;
	}
}

void Initia::set_module(std::string const & name)
{
	for (size_t i = 0; i < mlist.size(); i ++) if (mlist[i]->id == name) {
		active = i;
		ag.clear();
		ag.add(mlist[i]->param);
		return;
	}
	// show help
	bool req_help = name == "help";
	if (! req_help) std::cout << "Unknown module: " << name << '\n';
	std::cout << "\nList of available modules\n\n";
	for (size_t i = 0; i < mlist.size(); i ++) {
		std::cout << '\t' << mlist[i]->id << '\n';
	}
	std::cout << '\n';
	exit(req_help ? 0 : 1);
}

std::string const & Initia::get_module() const
{
	return mlist[active]->id;
}

bool Initia::alter()
{
	return ag.alter();
}

void Initia::initialize(MEAsim * sys)
{
	mlist[active]->initialize(sys);
}

std::vector<std::string> Initia::get_names() const
{
	std::vector<std::string> ns;
	for (size_t i = 0; i < mlist.size(); i ++) ns.push_back(mlist[i]->id);
	return ns;
}

prm::Param const & Initia::get_param() const
{
	return mlist[active]->param;
}

prm::Param & Initia::get_param()
{
	return mlist[active]->param;
}

bool Initia::call_module_line(int, const std::string & mod_line, void * data)
{
	Initia * ini = static_cast<Initia *>(data);
	std::string::size_type k = mod_line.find(':');
	std::string mod_name = mod_line.substr(0, k);
	ini->set_module(mod_name);
	if (k != std::string::npos) {
		arg::SubParser * sp = ini->ag.make_parser();
		sp->set(mod_line.substr(k+1));
		delete sp;
	}
	return true;
}
