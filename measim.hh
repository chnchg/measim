/*
 * Simple construction of MEAsim's model
 */
#ifndef MEASIM_HH
#define MEASIM_HH
#include <prm/ctrl_tag.hh>
#include <sim/with_param.hh>
#include <sim/mt19937.hh>
#include <sigc.hh>
#if __cplusplus < 201103L
#define override
#define final
#endif
enum TypeID {inhi = -1, exci = 1};   //By Justin, 20141208
class MEAsim :
	virtual protected WithParam
{
	MT19937 rng;
	double age;
	// parameters
	double g_Ca;
	double g_K;
	double g_L;

	double V_Ca;
	double V_K;
	double V_L;

	double V_1;
	double V_2;
	double V_3;
	double V_4;

	double V_r;
	double V_i;    // reversal potential for inhibitory connection, by Justin, 20141208
	double C;
	double theta;
//	double I_bg;
	double I_bgm;
	double I_bgs;

	double tau_d;
	double tau_r;
	double tau_l;
	double tau_s;
	double u_lower;
	double u_upper;

	double xi_ave;
	double eta_max;
	double k_a;
	int mexp;

	double beta;
	double k_r;
	int nexp;
	double I_p;
	double gamma;
	double Ca2p_0;

	double threshold; // firing criteria voltage threshold

	// control
	bool syn_inh;
	// synaptic weight modulation
	double wscaleE; // for excitatory synapses
	double wscaleI; // for inhibitory synapses
	double range;
public:
	struct Neuron
	{
		double V;
		double W;
		double Ca2p_r;
		int type{exci};
	};
	struct Synapse
	{
		double Y;
		double Z;
		double S;
	};
protected:
	size_t nnum;
	size_t snum;
	Neuron * nrns;
	Synapse * syns;
	double * rndu;
	double * rndi;
public:
	// network
	struct Node
	{
		double x;
		double y;
		// auxiliary
		std::vector<size_t> in;
		std::vector<size_t> out;
	};
	struct Link
	{
		double weight;
		size_t from;
		size_t to;
		// auxiliary
		double distance;
		int syn_type;
	};
protected:
	Node * nodes;
	Link * links;
	// workspace
	Neuron * vnrn;
	Synapse * vsyn;
	bool * firing;
protected:
	void make_auxiliary();
	prm::Compound nrn;
	prm::Compound syn;
	prm::Strt * res_nrns(size_t n) {return new prm::Stt<Neuron>(nrns[n]);}
	prm::Strt * res_node(size_t n) {return new prm::Stt<Node>(nodes[n]);}
	prm::Strt * res_syns(size_t s) {return new prm::Stt<Synapse>(syns[s]);}
	prm::Strt * res_link(size_t s) {return new prm::Stt<Link>(links[s]);}
public:
	MEAsim();
	~MEAsim();
	void reset(); // reset the state and age

	// setup
	void set_nnum(size_t nn);
	void set_snum(size_t ns);
	void set_size(size_t nn, size_t ns);
	void place(size_t n, double x, double y);
	// void set_link(size_t idx, size_t n0, size_t n1, double weight);
	void set_link(size_t idx, size_t n0, size_t n1, double weight, TypeID id = exci); // By Justin, 20141208
	void set_nrnType(size_t idx, TypeID id);                                         // By Justin, 20150113
	// integration
	void make_slope(Neuron * vn = 0, Synapse * vs = 0);
	void step_forward(double time, Neuron * vn = 0, Synapse * vs = 0, bool * fr = 0);
	void step(double time);

	// member access
	double get_age() const {return age;}
	size_t get_nnum() const {return nnum;}
	size_t get_snum() const {return snum;}
	double get_wse() const {return wscaleE;}
	double get_wsi() const {return wscaleI;}
	double get_vr() const {return V_r;}
	double get_vi() const {return V_i;}
	bool const * get_firing() const {return firing;}

	// file storage
	void h5_save(hid_t group) const;
	void h5_load(hid_t group);

	// state access
	Neuron const & get_neuron(size_t n) const {return nrns[n];}
	Node const & get_node(size_t n) const {return nodes[n];}
	Synapse const & get_synapse(size_t l) const {return syns[l];}
	Link const & get_link(size_t l) const {return links[l];}
	double get_wmod(size_t l) const {
		double wscale = 0.0;
		(links[l].syn_type==exci)? wscale = wscaleE: wscale = wscaleI;
		if (range) return wscale * links[l].weight * exp(- links[l].distance / range);
		return wscale * links[l].weight;
	}

	// query
	double active_over(size_t n) const {return nrns[n].V > threshold ? nrns[n].V - threshold : 0;}

	// notifications
	sigc::signal<void, size_t> new_nnum;
	sigc::signal<void, size_t> new_snum;
	sigc::signal<void, size_t> fires;
};
#endif
