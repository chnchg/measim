#ifndef MEASIM_RUN_HH
#define MEASIM_RUN_HH
#include "initia.hh"
#include <sim/runnable.hh>
#include <sim/parser_add.hh>
#include <prm/argu.hh>
#include <sigc.hh>
#include "h5_savedata_file.hh"

class MEAsimRun :
	virtual public Runnable,
	virtual public MEAsim,
	virtual protected ParserAdd
{
protected:
	// running parameter
	double time_step;
	bool dump_fire;
	bool dump_memp;
	bool use_data_file;
	//	bool dump_states;
	double dump_interval;
	std::string out_file;		// specify the file for storing the membrane potentials of neurons and spiking stamps
	H5SimSaveFile sim_data_file;

	// utilities
	Initia initia;

public:
	// overrides
	void init() override;
	void dash() override;
	void initH5OutFile();
	void h5_save(hid_t group) const override;
	void h5_load(hid_t group) override;
	double count() const override;
	void dump(int what, std::ostream & output) const override;
	void addto_parser(arg::Parser & parser) override;
	bool alter() override;

	MEAsimRun();
	Initia & get_initia() {return initia;}
	size_t get_sz() const; // number of neurons
	size_t get_ln() const; // number of synapses
	double get_pot(size_t n) const;
	double get_rst(size_t n) const;
	double get_cal(size_t n) const;
	int get_type(size_t n) const;
	double get_x(size_t n) const;
	double get_y(size_t n) const;
	double get_syny(size_t n) const;
	double get_synz(size_t n) const;
	double get_syns(size_t n) const;
	size_t get_from(size_t l) const;
	size_t get_to(size_t l) const;
	sigc::signal<void, MEAsim *> setup; // Should be removed?
};
#endif
