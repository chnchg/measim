#include "measim.hh"
#include <prm/hdf5.hh>
#include <limits>
#include <cmath>
void MEAsim::make_auxiliary()
{
	for (size_t n = 0; n < nnum; n ++) {
		nodes[n].in.clear();
		nodes[n].out.clear();
	}
	for (size_t l = 0; l < snum; l ++) {
		Link & lk = links[l];
		Node & n0 = nodes[lk.from];
		Node & n1 = nodes[lk.to];
		n1.in.push_back(lk.from);
		n0.out.push_back(lk.to);
		double dx = n0.x - n1.x;
		double dy = n0.y - n1.y;
		lk.distance = sqrt(dx * dx + dy * dy);
	}
}

MEAsim::MEAsim() :
	nnum(0), snum(0)
{
	using namespace prm::tag;
	param.add_var("g_Ca", g_Ca = 1.1)
		<< Name("g_Ca")
		<< Desc("Calcium channel conductance (mS)")
		<< Range(0, 8)
		<< Step(0.01, 0.1)
		<< Save()
		<< CmdLine();
	param.add_var("g_K", g_K = 2)
		<< Name("g_K")
		<< Desc("Potassium channel conductance (mS)")
		<< Range(0, 8)
		<< Step(0.01, 0.1)
		<< Save()
		<< CmdLine();
	param.add_var("g_L", g_L = 0.5)
		<< Name("g_L")
		<< Desc("Leaking channel conductance (mS)")
		<< Range(0, 8)
		<< Step(0.01, 0.1)
		<< Save()
		<< CmdLine();
	param.add_var("V_Ca", V_Ca = 100)
		<< Name("V_Ca")
		<< Desc("Calcium channel reverse potential (mV)")
		<< Range(0, 200)
		<< Step(1, 10)
		<< Save()
		<< CmdLine();
	param.add_var("V_K", V_K = - 70)
		<< Name("V_K")
		<< Desc("Potassium channel reverse potential (mV)")
		<< Range(- 100, - 20)
		<< Step(1, 10)
		<< Save()
		<< CmdLine();
	param.add_var("V_L", V_L = - 65)
		<< Name("V_L")
		<< Desc("Leaking channel reverse potential (mV)")
		<< Range(- 100, 20)
		<< Step(1, 10)
		<< Save()
		<< CmdLine();
	param.add_var("V_1", V_1 = - 1)
		<< Name("V_1")
		<< Desc("Threshold of calcium channel (mV)")
		<< Range(- 20, 20)
		<< Step(1, 10)
		<< Save()
		<< CmdLine();
	param.add_var("V_2", V_2 = 15)
		<< Name("V_2")
		<< Desc("Inverse sensitivity of calcium channel (mV)")
		<< Range(- 10, 30)
		<< Step(1, 10)
		<< Save()
		<< CmdLine();
	param.add_var("V_3", V_3 = 0)
		<< Name("V_3")
		<< Desc("Threshold of potassium channel (mV)")
		<< Range(- 20, 20)
		<< Step(1, 10)
		<< Save()
		<< CmdLine();
	param.add_var("V_4", V_4 = 30)
		<< Name("V_4")
		<< Desc("Inverse sensitivity of potassium channel (mV)")
		<< Range(10, 50)
		<< Step(1, 10)
		<< Save()
		<< CmdLine();
	param.add_var("V_r", V_r = 0)
		<< Name("V_r")
		<< Desc("Reversal potential for excitatory synaptic current (mV)")
		<< Range(-70, 20)
		<< Step(1, 10)
		<< Save()
		<< CmdLine();
	param.add_var("V_i", V_i = -90)
		<< Name("V_i")
		<< Desc("Reversal potential for inhibitory synaptic current (mV)")
		<< Range(-110, 0)
		<< Step(1, 10)
		<< Save()
		<< CmdLine();
	param.add_var("C", C = 1)
		<< Name("C")
		<< Desc("Membrane capacitance density (uF)")
		<< Range(0.01, 10)
		<< Step(0.01, 0.1)
		<< Save()
		<< CmdLine();
	param.add_var("theta", theta = 0.2)
		<< Name("theta")
		<< Desc("response speed of potassium channel (ms)^-1")
		<< Range(0.01, 5)
		<< Step(0.01, 0.1)
		<< Save()
		<< CmdLine();
	param.add_var("I_bgm", I_bgm = 27)
		<< Name("I_bgm")
		<< Desc("Mean value of background input current (uA)")
		<< Range(0, 100)
		<< Step(0.1, 1)
		<< Save()
		<< CmdLine();
	param.add_var("I_bgs", I_bgs = 1)
		<< Name("I_bgs")
		<< Desc("Standard deviation of background input current (uA)")
		<< Range(0, 100)
		<< Step(0.1, 1)
		<< Save()
		<< CmdLine();
	param.add_var("u_lower", u_lower = 0.2)
		<< Name("u_lower")
		<< Desc("Lower bound of the fraction of transmittre release upon action potential")
		<< Range(0, 1)
		<< Step(0.01, 0.1)
		<< Save()
		<< CmdLine();
	param.add_var("u_upper", u_upper = 0.2)
		<< Name("u_upper")
		<< Desc("Upper bound of the fraction of transmittre release upon action potential")
		<< Range(0, 1)
		<< Step(0.01, 0.1)
		<< Save()
		<< CmdLine();
	param.add_var("tau_d", tau_d = 10)
		<< Name("tau_d")
		<< Desc("Decay time of active transmitter (ms)")
		<< Range(0, 50)
		<< Step(1, 10)
		<< Save()
		<< CmdLine();
	param.add_var("tau_r", tau_r = 300)
		<< Name("tau_r")
		<< Desc("Recovery time of inactive transmitter (ms)")
		<< Range(10, 15000)
		<< Step(10, 1000)
		<< Save()
		<< CmdLine();
	param.add_var("tau_l", tau_l = 600)
		<< Name("tau_l")
		<< Desc("Leaking time to super-inactive state (ms)")
		<< Range(10, 15000)
		<< Step(10, 1000)
		<< Save()
		<< CmdLine();
	param.add_var("tau_s", tau_s = 5000) // 5 s
		<< Name("tau_s")
		<< Desc("Recovery time of super-inactive transmitter (ms)")
		<< Range(1000, 50000)
		<< Step(100, 1000)
		<< Save()
		<< CmdLine();
	param.add_var("xi_ave", xi_ave = 0.02)
		<< Name("xi_ave")
		<< Desc("Fractional size of asynchronous release")
		<< Range(0, 0.5)
		<< Step(0.0001, 0.01)
		<< Save()
		<< CmdLine();
	param.add_var("eta_max", eta_max = 0.32)
		<< Name("eta_max")
		<< Desc("Maximum asynchronize release rate (ms)^-1")
		<< Range(0, 2000)
		<< Step(0.01, 1)
		<< Save()
		<< CmdLine();
	param.add_var("k_a", k_a = 0.1)
		<< Name("k_a")
		<< Desc("Activation residual calcium level for asynchronous release (uM)")
		<< Range(0, 1)
		<< Step(0.01, 0.1)
		<< Save()
		<< CmdLine();
	param.add_var("m", mexp = 4)
		<< Name("m")
		<< Desc("Hill exponent of residual calcium induced asynchronous release")
		<< Range(- 10, 10)
		<< Step(1, 3)
		<< Save()
		<< CmdLine();
	param.add_var("beta", beta = 0.005)
		<< Name("beta")
		<< Desc("Maximum calcium extrusion rate (uM/ms)")
		<< Range(0.001, 0.100)
		<< Step(0.0005, 0.01)
		<< Save()
		<< CmdLine();
	param.add_var("k_r", k_r = 0.4)
		<< Name("k_r")
		<< Desc("Calcium extrusion pump activation level (uM)")
		<< Range(0, 1)
		<< Step(0.01, 0.1)
		<< Save()
		<< CmdLine();
	param.add_var("n", nexp = 2)
		<< Name("n")
		<< Desc("Hill exponent for residual calcium extrusion")
		<< Range(- 10, 10)
		<< Step(1, 3)
		<< Save()
		<< CmdLine();
	param.add_var("I_p", I_p = 0.00011)
		<< Name("I_p")
		<< Desc("Passive influx of calcium (uM/ms)")
		<< Range(0, 0.00099)
		<< Step(0.00001, 0.0001)
		<< Save()
		<< CmdLine();
	param.add_var("gamma", gamma = 0.0330)
		<< Name("gamma")
		<< Desc("Release constant for residual calcium")
		<< Range(0, 10)
		<< Step(0.0001, 0.1)
		<< Save()
		<< CmdLine();
	param.add_var("Ca2p_0", Ca2p_0 = 2000)
		<< Name("Ca2p_0")
		<< Desc("extra-synaptic calcium concentration (uM)")
		<< Range(100, 20000)
		<< Step(100, 1000)
		<< Save()
		<< CmdLine();
	param.add_var("threshold", threshold = 10)
		<< Name("threshold")
		<< Desc("voltage threshold for firing (mV)")
		<< Range(- 20, 30)
		<< Step(1, 10)
		<< Save()
		<< CmdLine();
	param.add_var("syn_inh", syn_inh = false)
		<< Name("syn_inh")
		<< Desc("inhibit synchronous release")
		<< Control()
		<< Save()
		<< CmdLine();
	param.add_var("wscalee", wscaleE = 4.0)
		<< Name("wscalee")
		<< Desc("overall scaling factor of excitatory synaptic weights")
		<< Range(0, 100.0)
		<< Step(0.01, 0.1)
		<< Save()
		<< CmdLine();
	param.add_var("wscalei", wscaleI = 4.0)
		<< Name("wscalei")
		<< Desc("overall scaling factor of inhibitory synaptic weights")
		<< Range(0, 100.0)
		<< Step(0.01, 0.1)
		<< Save()
		<< CmdLine();
	param.add_var("range", range = 0)
		<< Name("range")
		<< Desc("decay length of weight")
		<< Range(0, 100)
		<< Step(0.005, 1)
		<< Save()
		<< CmdLine();
	// properties
	param.add_var("age", age = 0)
		<< Name("age")
		<< Desc("Time since last reset")
		<< Save()
		<< CmdLine();

	nrn.set_sizer(new prm::SzrFun<MEAsim>(* this, & MEAsim::set_nnum, & MEAsim::get_nnum));
	syn.set_sizer(new prm::SzrFun<MEAsim>(* this, & MEAsim::set_snum, & MEAsim::get_snum));

	prm::Struct * ss = new prm::Struct;
	ss->add_member("mpot", & Neuron::V)
		<< Name("mpot")
		<< Desc("membrane potential");
	ss->add_member("W", & Neuron::W)
		<< Name("W")
		<< Desc("restoring field");
	ss->add_member("Ca2p_r", & Neuron::Ca2p_r)
		<< Name("Ca2p_r")
		<< Desc("residue calcium concentration");
	ss->add_member("type", & Neuron::type)
		<< Name("type")
		<< Desc("type of the neuron");
	ss->set_resolver(* this, & MEAsim::res_nrns);
	nrn.add_struct(ss);

	ss = new prm::Struct;
	ss->add_member("x", & Node::x)
		<< Name("x")
		<< Desc("x-coordinate");
	ss->add_member("y", & Node::y)
		<< Name("y")
		<< Desc("y-coordinate");
	ss->set_resolver(* this, & MEAsim::res_node);
	nrn.add_struct(ss);

	nrn.add_linear("rndu",rndu)
		<< Name("rndu")
		<< Desc("heterogeneity in u");
	nrn.add_linear("rndi",rndi)
		<< Name("rndi")
		<< Desc("heterogeneity in I");

	ss = new prm::Struct;
	ss->add_member("Y", & Synapse::Y)
		<< Name("Y")
		<< Desc("active transmitter fraction");
	ss->add_member("Z", & Synapse::Z)
		<< Name("Z")
		<< Desc("inactive transmitter fraction");
	ss->add_member("S", & Synapse::S)
		<< Name("S")
		<< Desc("super-inactive transmitter fraction");
	ss->set_resolver(* this, & MEAsim::res_syns);
	syn.add_struct(ss);

	ss = new prm::Struct;
	ss->add_member("weight", & Link::weight)
		<< Name("weight")
		<< Desc("synaptic weight");
	ss->add_member("from", & Link::from)
		<< Name("from")
		<< Desc("from neuron");
	ss->add_member("to", & Link::to)
		<< Name("to")
		<< Desc("to neuron");
	ss->add_member("type", & Link::syn_type)
		<< Name("type")
		<< Desc("type of synapse");
	ss->set_resolver(* this, & MEAsim::res_link);
	syn.add_struct(ss);
}

MEAsim::~MEAsim()
{
	if (nnum) {
		delete [] nrns;
		delete [] nodes;
		delete [] vnrn;
		delete [] firing;
		delete [] rndu;
		delete [] rndi;
	}
	if (snum) {
		delete [] syns;
		delete [] links;
		delete [] vsyn;
	}
}

void MEAsim::reset()
{
	age = 0;
	// find steady state for an isolated neuron
	double V_ss = V_L + g_L * I_bgm;
	double W_ss;
	double ovss;
	do {
		ovss = V_ss;
		W_ss = (1.0 + tanh((V_ss - V_3) / V_4)) / 2;
		double m_inf = (1 + tanh((V_ss - V_1) / V_2)) / 2;
		V_ss = I_bgm + g_Ca * m_inf * V_Ca + g_K * W_ss * V_K + g_L * V_L;
		V_ss /= g_Ca * m_inf + g_K * W_ss + g_L;
		V_ss = (V_ss + 5 * ovss) / 6;
	} while (fabs(V_ss - ovss) > 1e-15);
	double Ca2p_ss = k_r * pow(I_p / (beta - I_p), 1.0 / nexp);
	for (size_t i = 0; i < nnum; i ++) {
		nrns[i].V = V_ss;
		nrns[i].W = W_ss;
		nrns[i].Ca2p_r = Ca2p_ss;
		rndu[i] = rng.uniform();
		rndi[i] = rng.gaussian();
	}

	double Y_ss = 0;
	double Z_ss = 0;
	double S_ss = 0.8;
	for (size_t i = 0; i < snum; i ++) {
		syns[i].Y = Y_ss;
		syns[i].Z = Z_ss;
		syns[i].S = S_ss;
	}
	make_auxiliary();
}

void MEAsim::set_nnum(size_t nn)
{
	if (nnum) {
		delete [] nrns;
		delete [] nodes;
		delete [] vnrn;
		delete [] firing;
		delete [] rndu;
		delete [] rndi;
	}
	nrns = new Neuron [nn];
	nodes = new Node [nn];
	vnrn = new Neuron [nn];
	rndu = new double [nn];
	rndi = new double [nn];
	firing = new bool [nn];
	nnum = nn;
	new_nnum(nn);
}

void MEAsim::set_snum(size_t ns)
{
	if (snum) {
		delete [] syns;
		delete [] links;
		delete [] vsyn;
	}
	syns = new Synapse [ns];
	links = new Link [ns];
	vsyn = new Synapse [ns];
	snum = ns;
	new_snum(ns);
}

void MEAsim::set_size(size_t nn, size_t ns)
{
	set_nnum(nn);
	set_snum(ns);
}

void MEAsim::place(size_t n, double x, double y)
{
	nodes[n].x = x;
	nodes[n].y = y;
}


/*void MEAsim::set_link(size_t idx, size_t n0, size_t n1, double weight)
{
	Link l = {weight, n0, n1, 0, exci};
	links[idx] = l;
	// nodes[n0].out.push_back(idx);
	// nodes[n1].in.push_back(idx);
	}*/

// By Justin, 20141208
void MEAsim::set_link(size_t idx, size_t n0, size_t n1, double weight, TypeID syn_type)
{
	Link l = {weight, n0, n1, 0, syn_type};
	links[idx] = l;
	// nodes[n0].out.push_back(idx);
	// nodes[n1].in.push_back(idx);
}
void MEAsim::set_nrnType(size_t idx, TypeID nrn_type)
{
	nrns[idx].type = nrn_type;
}

void MEAsim::make_slope(Neuron * vn, Synapse * vs)
{
	if (! vn) vn = vnrn;
	if (! vs) vs = vsyn;
	double krn = pow(k_r, nexp);
	// calculate input
	double * exci_inp = new double [nnum];
	double * inhi_inp = new double [nnum];
	for (size_t i = 0; i < nnum; i ++) {
		exci_inp[i] = 0;
		inhi_inp[i] = 0;}
	for (size_t i = 0; i < snum; i ++) {
		Link & l = links[i];
		if (range) {
			if (l.syn_type == exci) exci_inp[l.to] += syns[i].Y * l.weight * exp(- l.distance / range);
			else inhi_inp[l.to] += syns[i].Y * l.weight * exp(- l.distance / range);
		}
		else {
			if (l.syn_type == exci)	exci_inp[l.to] += syns[i].Y * l.weight;
			else inhi_inp[l.to] += syns[i].Y * l.weight;
		}
	}
	for (size_t i = 0; i < nnum; i ++) {
		exci_inp[i] *= wscaleE;
		inhi_inp[i] *= wscaleI;}
	for (size_t i = 0; i < nnum; i ++) {
		Neuron & n = nrns[i];
		double m_inf = (1 + tanh((n.V - V_1) / V_2)) / 2;
		double I_ion = g_Ca * m_inf * (n.V - V_Ca) + g_K * n.W * (n.V - V_K) + g_L * (n.V - V_L);
		vn[i].V = (- I_ion + exci_inp[i] * (V_r - n.V) + inhi_inp[i] * (V_i - n.V) + I_bgm+I_bgs*rndi[i]) / C;
		double tau_W = 1 / cosh((n.V - V_3) / (2 * V_4));
		vn[i].W = theta * ((1 + tanh((n.V - V_3) / V_4)) / 2 - n.W) / tau_W;
		// calcium
		double Cn = pow(n.Ca2p_r, nexp);
		vn[i].Ca2p_r = - beta * Cn / (krn + Cn) + I_p;
	}
	for (size_t i = 0; i < snum; i ++) {
		Synapse & s = syns[i];
		vs[i].Y = - s.Y / tau_d;
		vs[i].Z = s.Y / tau_d - s.Z / tau_r - s.Z / tau_l;
		vs[i].S = s.Z / tau_l - s.S / tau_s;
	}
	delete [] exci_inp;
	delete [] inhi_inp;
}
void MEAsim::step_forward(double time, Neuron * vn, Synapse * vs, bool * fr)
{
	if (! vn) vn = vnrn;
	if (! vs) vs = vsyn;
	if (! fr) fr = firing;
	double * eta = new double [nnum];
	double kam = pow(k_a, mexp);
	for (size_t i = 0; i < nnum; i ++) {
		bool sub = nrns[i].V < threshold;
		nrns[i].V += vn[i].V * time;
		nrns[i].W += vn[i].W * time;
		nrns[i].Ca2p_r += vn[i].Ca2p_r * time;
		fr[i] = sub && nrns[i].V >= threshold;
		double Cm = pow(nrns[i].Ca2p_r, mexp);
		eta[i] = eta_max * Cm / (kam + Cm);
	}
	double teta = 0;
	double * acc_rate = new double [snum];
	for (size_t i = 0; i < snum; i ++) {
		syns[i].Y += vs[i].Y * time;
		syns[i].Z += vs[i].Z * time;
		syns[i].S += vs[i].S * time;
		acc_rate[i] = teta;
		teta += eta[links[i].from];
	}
	// synchronous release synapse
	if (! syn_inh) for (size_t i = 0; i < snum; i ++) if (fr[links[i].from]) {
		Synapse & s = syns[i];
		s.Y += ((u_upper-u_lower)*rndu[links[i].from]+u_lower) * (1 - s.Y - s.Z - s.S);
	}
	// synchronous release calcium
	for (size_t n = 0; n < nnum; n ++) if (fr[n]) {
		double & c = nrns[n].Ca2p_r;
		c += gamma * log(Ca2p_0 / c);
		// emit firing signal
		fires(n);
	}
	// asynchronous release
	double en = exp(- teta * time);
	for (double p = rng.uniform(); p > en; p *= rng.uniform()) {
		// random pick following eta
		double r = rng.uniform() * teta;
		size_t i = 0;
		size_t j = snum;
		// Invariant: r < acc_rate[j] && acc_rate[i] < r
		while (i + 1 < j) {
			size_t k = (i + j) / 2;
			if (r < acc_rate[k]) j = k;
			else i = k;
		}
		Synapse & s = syns[i];
		s.Y += xi_ave * (1 - s.Y - s.Z - s.S);
	}
	delete [] eta;
	delete [] acc_rate;
}
void MEAsim::step(double time)
{
	make_slope();
	step_forward(time);
	age += time;
}

void MEAsim::h5_save(hid_t group) const
{
	rng.h5_save(group);
	prm::h5_save(group, param);
	prm::h5_save(group, nrn, "neurons");
	prm::h5_save(group, syn, "synapses");
}

void MEAsim::h5_load(hid_t group)
{
	rng.h5_load(group);
	prm::h5_load(group, param);
	prm::h5_load(group, nrn, "neurons");
	prm::h5_load(group, syn, "synapses");
	make_auxiliary();
}
