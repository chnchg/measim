#include "measim_run.hh"
#include <prm/hdf5.hh>
#include <vector>
#include <cassert>

void MEAsimRun::init()
{
	// setup.emit(this);
	initia.initialize(this);
	reset();
}

void MEAsimRun::initH5OutFile()
{
	sim_data_file.set_name(out_file);
	sim_data_file.set_nnum(get_sz());
	sim_data_file.initH5simSaveFile();
}

void MEAsimRun::dash()
{
	static double next_dump_time = 0;
	static bool print_head = true;
	size_t nn = get_nnum();	
	size_t sn = get_snum();
	step(time_step);
	std::cout.precision(9);
	if (use_data_file && (dump_fire || dump_memp)) {
		if (! sim_data_file.get_isOpen()){
			initH5OutFile();
		}
		else if (out_file != sim_data_file.get_name() && sim_data_file.get_isOpen()){
			sim_data_file.done_save();
			initH5OutFile();
		}
	}
	
	if (dump_fire) {
		bool const * fr = get_firing();
		if (use_data_file) sim_data_file.dump_fire(fr, get_age());
		else for (size_t n = 0; n < nn; n ++) if (fr[n]) std::cout << get_age() << '\t' << n << '\n';
	}

	if (dump_memp && get_age() > next_dump_time){
		while (get_age() > next_dump_time ) next_dump_time += dump_interval;
		if (use_data_file) {
			double * pot = new double[get_sz()];
			for (size_t n = 0; n < nn; n++) pot[n] = get_pot(n);
			sim_data_file.dump_mpot(pot);
			delete [] pot;
		}
		else {
			for (size_t n = 0; n < nn; n++) std::cout << get_pot(n) << "\t";
			std::cout << std::endl;
		}
	}
}

void MEAsimRun::h5_save(hid_t group) const
{
	MEAsim::h5_save(group);

	hid_t g = H5Gcreate(group, "init_module", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t stype = H5Tcopy(H5T_C_S1);
	H5Tset_size(stype, H5T_VARIABLE);
	hid_t space = H5Screate(H5S_SCALAR);
	std::string const & mn = initia.get_module();
	hid_t attr = H5Acreate(g, "mod_name", stype, space, H5P_DEFAULT, H5P_DEFAULT);
	char const * s = mn.c_str();
	H5Awrite(attr, stype, & s);
	prm::h5_save(g, initia.get_param());
	H5Aclose(attr);
	H5Sclose(space);
	H5Tclose(stype);
	H5Gclose(g);
}

void MEAsimRun::h5_load(hid_t group)
{
	MEAsim::h5_load(group);

	H5Eset_auto(H5E_DEFAULT, 0, 0); // turn off error reporting
	bool has_init = H5Gget_objinfo(group, "init_module", 0, 0) >= 0;
        H5Eset_auto(H5E_DEFAULT, (H5E_auto_t) H5Eprint, stderr); // turn on default err
	if (! has_init) {
		H5Eclear(H5E_DEFAULT);
		std::cerr << "No init_module information in data file!\n";
		return;
	}
	hid_t g = H5Gopen(group, "init_module", H5P_DEFAULT);
	hid_t stype = H5Tcopy(H5T_C_S1);
	H5Tset_size(stype, H5T_VARIABLE);
	hid_t attr = H5Aopen_name(g, "mod_name");
	char * mn;
	H5Aread(attr, stype, & mn);
	hid_t space = H5Aget_space(attr);
	assert(H5Sget_simple_extent_ndims(space) == 0);
	initia.set_module(mn);
	H5Dvlen_reclaim(stype, space, H5P_DEFAULT, & mn);
	H5Sclose(space);
	H5Aclose(attr);
	prm::h5_load(g, initia.get_param());
	H5Tclose(stype);
	H5Gclose(g);
}

double MEAsimRun::count() const
{
	return get_age();
}

void MEAsimRun::dump(int what, std::ostream & output) const
{
	if (dump_set.check(what, "net_adj")) {
		size_t z = get_nnum();
//		bool * adj = new bool [z * z];
		double * adj = new double [z * z];
		for (size_t i = 0; i < z * z; i ++) adj[i] = 0;
//		for (size_t l = 0; l < get_snum(); l ++) adj[z * get_link(l).from + get_link(l).to] = true;
		for (size_t l = 0; l < get_snum(); l ++) adj[z * get_link(l).from + get_link(l).to] = get_link(l).weight*get_link(l).syn_type;
		for (size_t i = 0; i < z; i ++) {
			for (size_t j = 0; j < z; j ++) {
//				output << ' ' << (adj[i * z + j] ? '1' : '0');
				output << adj[i * z + j] << "\t";
			}
			output << '\n';
		}
		delete [] adj;
	}
	else if (dump_set.check(what, "kdist")) {
		size_t z = get_nnum();
		size_t * k = new size_t [z];
		for (size_t i = 0; i < z; i ++) k[i] = 0;
		for (size_t l = 0; l < get_snum(); l ++) {
			k[get_link(l).from] ++;
			k[get_link(l).to] ++;
		}
		size_t kmax = 0;
		for (size_t i = 0; i < z; i ++) if (k[i] > kmax) kmax = k[i];
		kmax ++;
		size_t * dk = new size_t [kmax];
		for (size_t kk = 0; kk < kmax; kk ++) dk[kk] = 0;
		for (size_t i = 0; i < z; i ++) dk[k[i]] ++;
		size_t tn = z;
		for (size_t kk = 0; kk < kmax; kk ++) {
			output << kk << '\t' << dk[kk] << '\t' << tn << '\n';
			tn -= dk[kk];
		}
		output << kmax << '\t' << 0 << '\t' << 0 << '\n';
		delete [] dk;
		delete [] k;
	}
	else if (dump_set.check(what, "kacc")) {
		size_t z = get_nnum();
		size_t * k = new size_t [z];
		for (size_t i = 0; i < z; i ++) k[i] = 0;
		for (size_t l = 0; l < get_snum(); l ++) {
			k[get_link(l).from] ++;
			k[get_link(l).to] ++;
		}
		size_t kmax = 0;
		for (size_t i = 0; i < z; i ++) if (k[i] > kmax) kmax = k[i];
		kmax ++;
		size_t * dk = new size_t [kmax];
		for (size_t kk = 0; kk < kmax; kk ++) dk[kk] = 0;
		for (size_t i = 0; i < z; i ++) dk[k[i]] ++;
		size_t tn = z;
		if (! dk[0]) output << 0 << '\t' << 1 << '\n';
		for (size_t kk = 0; kk < kmax; kk ++) if (dk[kk]) {
			output << kk << '\t' << double(tn) / z << '\n';
			tn -= dk[kk];
			output << kk << '\t' << double(tn) / z << '\n';
		}
		delete [] dk;
		delete [] k;
	}
	else if (dump_set.check(what, "info")) {
		output << "# \n";
		output << "# init parameters\n";
		output << "# \n";
		output << "module=" << initia.get_module() << '\n';
		initia.get_param().write(output);
	}
}

void MEAsimRun::addto_parser(arg::Parser & parser)
{
	parser.add_opt('o', "outfile")
		.stow(out_file)
		.help("specify the file name for storing data", "filename")
		.show_default();

	arg::Option * io = parser.find("init");  // use existing --init if available
	if (! io) io = & parser.add_opt("init");
	io->store() // null storage
		.call(& Initia::call_module_line, & initia)
		.help("network setup, `help' for usage", "MODULE[:PARAM[,PARAM...]]");
}

bool MEAsimRun::alter()
{
	return initia.alter();
}

MEAsimRun::MEAsimRun()
{
	using namespace prm::tag;
	param.add_var("time_step", time_step = 0.01)
		<< Name("Time Step")
		<< Desc("Time step size for numerical integration (ms)")
		<< Range(0.0001, 1)
		<< Step(0.0001, 0.01)
		<< Save()
		<< CmdLine();
 	param.add_var("dump_interval", dump_interval = 1.0)
		<< Name("dump interval")
		<< Desc("dump interval size for dumping memp or states (ms)")
		<< Range(0.0001, 10)
		<< Step(0.0001, 0.01)
		<< Save()
		<< CmdLine();
        param.add_var("dump_fire", dump_fire = false)
                << Name("dump_fire")
                << Desc("dump firing message to stdout")
                << Control()
//                << Save()
                << CmdLine();

       param.add_var("dump_memp", dump_memp = false)
		   << Name("dump_memp")
		   << Desc("dump membrane potential to stdout")
		   << Control()
		   // << Save()
		   << CmdLine();
	   param.add_var("use_data_file", use_data_file = false)
		   << Name("use_data_file")
		   << Desc("use a file for dumping data")
		   << CmdLine();
       /*
       param.add_var("dump_states", dump_states = false)
                << Name("dump_states")
                << Desc("dump states to stdout")
                << Control()
//                << Save()
                << CmdLine();
       */
	out_file = "simData.h5";
	dump_set.add("net_adj", "adjacency matrix of network");
	dump_set.add("kdist", "degree distribution");
	dump_set.add("kacc", "acc degree distribution for plotting");
}

size_t MEAsimRun::get_sz() const
{
	return get_nnum();
}

size_t MEAsimRun::get_ln() const
{
	return get_snum();
}

double MEAsimRun::get_pot(size_t n) const
{
	return get_neuron(n).V;
}

double MEAsimRun::get_rst(size_t n) const
{
	return get_neuron(n).W;
}

double MEAsimRun::get_cal(size_t n) const
{
	return get_neuron(n).Ca2p_r;
}
int MEAsimRun::get_type(size_t n) const
{
	return get_neuron(n).type;
}
double MEAsimRun::get_x(size_t n) const
{
	return get_node(n).x;
}

double MEAsimRun::get_y(size_t n) const
{
	return get_node(n).y;
}

size_t MEAsimRun::get_from(size_t l) const
{
	return get_link(l).from;
}

size_t MEAsimRun::get_to(size_t l) const
{
	return get_link(l).to;
}
double MEAsimRun::get_syny(size_t l) const
{
	return get_synapse(l).Y;
}
double MEAsimRun::get_synz(size_t l) const
{
	return get_synapse(l).Z;
}
double MEAsimRun::get_syns(size_t l) const
{
	return get_synapse(l).S;
}
