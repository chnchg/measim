#include "mod_init_modules.hh"
#include <sim/h5_file.hh>
#include <fstream>
#include <algorithm>
#include <cmath>
namespace {
	struct Pz
	{
		double x;
		double y;
		Pz(double x, double y) : x(x), y(y) {}
		double r(Pz const & p)
		{
			double dx = p.x - x;
			double dy = p.y - y;
			return sqrt(dx * dx + dy * dy);
		}
	};

	struct Lk
	{
		size_t n0;
		size_t n1;
		TypeID syn_type;
		Lk(size_t n0, size_t n1) : n0(n0), n1(n1), syn_type(exci) {}
		Lk(size_t n0, size_t n1, TypeID syn_type) : n0(n0), n1(n1), syn_type(syn_type) {}
	};
}

using namespace mod_init;

double mod_init::WDev::make_weight()
{
	return wdev == 0 ? weight : weight + wdev * (rng.uniform() - 0.5);
}

double mod_init::WDev::make_gaussian()
{
	if (wdev == 0) return weight;
	double w;
	do w = weight + wdev * rng.gaussian();
	while (w < weight - wdev || w > weight + wdev);
	return w;
}

mod_init::WDev::WDev() :
	InitModule("wdev")
{
	using namespace prm::tag;
	param.add_var("weight", weight = 1)
		<< Name("weight")
		<< Desc("mean synaptic weight")
		<< Range(0, 1000)
		<< Step(0.1, 1)
		<< CmdLine()
		<< Save();
	param.add_var("wdev", wdev = 0.2)
		<< Name("wdev")
		<< Desc("weight spread")
		<< Range(0, 1.0)
		<< Step(0.01, 0.1)
		<< CmdLine()
		<< Save();
}

void mod_init::AddInhiNrn::make_inhi_nrns(int nnrns)
{
	int toal_inhi_nrns = std::ceil(inhiratio * nnrns);
	inhi_nrns.clear();
	inhi_nrns.resize(nnrns, false);
	for (int i = 0; i < toal_inhi_nrns; i ++) {
		int j;
		while (inhi_nrns[j = int(rng.uniform() * nnrns)]);
		inhi_nrns[j] = true;
	}
	return;
}

bool mod_init::AddInhiNrn::is_inhi_nrn(int nrn_idx)
{
	return inhi_nrns[nrn_idx];
}

mod_init::AddInhiNrn::AddInhiNrn():
	InitModule("addinhinrn")
{
	using namespace prm::tag;
	param.add_var("inhiratio", inhiratio = 0)
		<< Name("inhiratio")
		<< Desc("the percentage of inhibitory neurons in the network")
		<< Range(0, 1)
		<< Step(0.01, 0.1)
		<< CmdLine()
		<< Save();

}

void mod_init::Lattice::initialize(MEAsim * sys)
{
	rng.init(seed);
	size_t sz2 = lsz * lsz;
	size_t lsz4 = lsz/2;
	size_t sz4 = lsz4 * lsz4;
	bool * inhi_nrns = new bool[sz2+sz4];
	for (size_t n0 = 0; n0 < sz2 + sz4; n0++){
		if (n0>=sz2) inhi_nrns[n0]=true;
		else inhi_nrns[n0]=false;
	}
	
	std::vector<Pz> pzs;
	double ds1 = 1./double(lsz);
	for (size_t n0 = 0; n0 < lsz; n0 ++) for (size_t n1 = 0; n1 < lsz; n1 ++) pzs.push_back(Pz(ds1*n0,ds1*n1));
	double ds2 = ds1*2;
	for (size_t n0 = 0; n0 < lsz4; n0 ++) for (size_t n1 = 0; n1 < lsz4; n1 ++) pzs.push_back(Pz(ds2*n0+ds1/2,ds2*n1+ds1/2));
	
	std::vector<Lk> lks;
	size_t j;
	for (size_t n0 = 0; n0 < sz2+sz4; n0 ++) for (size_t n1 = n0 + 1; n1 < sz2+sz4; n1 ++) {
		if (pzs[n0].r(pzs[n1]) > range * ds1) continue;
		
		if (!inhi_nrns[n0] && rng.uniform() < rewiredP){
			while (true){
				j = rng.uniform()*(sz2+sz4);
				if (j != n1) break;
			}
		}
		else j=n1;
		if (inhi_nrns[n0])
			lks.push_back((Lk(n0,j,inhi)));
		else
			lks.push_back((Lk(n0,j,exci)));

		if (!inhi_nrns[n1] && rng.uniform() < rewiredP){
			while (true){
				j = rng.uniform()*(sz2+sz4);
				if (j != n0) break;
			}
		}
		else j=n0;
		
		if (inhi_nrns[n1])
			lks.push_back((Lk(n1,j,inhi)));
		else
			lks.push_back((Lk(n1,j,exci)));
	}
	
	sys->set_size(sz2+sz4, lks.size());
	for (size_t n = 0; n < sz2+sz4; n ++) {
		sys->place(n, pzs[n].x, pzs[n].y);
		inhi_nrns[n]?sys->set_nrnType(n,inhi):sys->set_nrnType(n,exci);
	}
	for (size_t i = 0; i < lks.size(); i ++) {
		sys->set_link(i, lks[i].n0, lks[i].n1, make_gaussian(), lks[i].syn_type);
	}
	delete[] inhi_nrns;

}

mod_init::Lattice::Lattice() :
	InitModule("lattice")
{
	using namespace prm::tag;
	param.add_var("lsz", lsz = 8)
		<< Name("lsz")
		<< Desc("linear size")
		<< Range(2, 256)
		<< Step(1, 10)
		<< CmdLine()
		<< Save();
	param.add_var("range", range = 1)
		<< Name("range")
		<< Desc("range of connectivity")
		<< Range(1, 10)
		<< Step(0.125, 2)
		<< CmdLine()
		<< Save();
	param.add_var("rewiredP", rewiredP = 0.0)
		<< Name("rewiredP")
		<< Desc("rewired probability")
		<< Range(0, 1)
		<< Step(0.001, 0.05)
		<< CmdLine()
		<< Save();}


void mod_init::Spaced::place(MEAsim * sys)
{
	sys->set_nnum(nnum);
	double md2 = mdisf / nnum;
	std::vector<Pz> pzs;
 do_over:
	for (size_t n = 0; n < nnum; n ++) {
		double d2;
		double x, y;
		size_t cnt = 0;
		do {
			x = rng.uniform();
			y = rng.uniform();
			d2 = 1;
			for (size_t m = 0; m < n; m ++) {
				double dx = pzs[m].x - x;
				double dy = pzs[m].y - y;
				double nd2 = dx * dx + dy * dy;
				if (nd2 < d2) d2 = nd2;
			}
			cnt ++;
			if (cnt > 100000) {
				pzs.clear();
				std::cerr << "do_over!\n!";
				goto do_over;
			}
		} while (d2 < md2);
		pzs.push_back(Pz(x, y));
	}
	for (size_t i = 0; i < pzs.size(); i ++) {
		sys->place(i, pzs[i].x, pzs[i].y);
	}
}

mod_init::Spaced::Spaced() :
	InitModule("spaced")
{
	using namespace prm::tag;
	param.add_var("nnum", nnum = 64)
		<< Name("nnum")
		<< Desc("number of neurons")
		<< Range(1, 4096)
		<< Step(1, 20)
		<< CmdLine()
		<< Save();	
	param.add_var("mdisf", mdisf = 0.5)
		<< Name("mdisf")
		<< Desc("minial distance factor")
		<< Range(0, 1.0)
		<< Step(0.01, 0.1)
		<< CmdLine()
		<< Save();
}

void mod_init::Random::initialize(MEAsim * sys)
{
	rng.init(seed);
	place(sys);
	make_inhi_nrns(nnum);
	std::vector<Lk> lks;
	if (bidir) for (size_t n0 = 0; n0 < nnum; n0 ++) {
		for (size_t n1 = n0 + 1; n1 < nnum; n1 ++)
		if (rng.uniform() < prob) {
			//lks.push_back(Lk(n0, n1));
			//lks.push_back(Lk(n1, n0));
			is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
			is_inhi_nrn(n1)?lks.push_back(Lk(n1, n0,inhi)):lks.push_back(Lk(n1, n0, exci));
		}
	}
	else for (size_t n0 = 0; n0 < nnum; n0 ++) {
		for (size_t n1 = n0 + 1; n1 < nnum; n1 ++) {
			//if (rng.uniform() < prob) lks.push_back(Lk(n0, n1));
			//if (rng.uniform() < prob) lks.push_back(Lk(n1, n0));
			if (rng.uniform() < prob) is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
			if (rng.uniform() < prob) is_inhi_nrn(n1)?lks.push_back(Lk(n1, n0,inhi)):lks.push_back(Lk(n1, n0, exci));
		}
	}
	// std::cerr << "creating " << lks.size() << " links.\n";
	sys->set_snum(lks.size());
	for (size_t i = 0; i < lks.size(); i ++) {
		sys->set_link(i, lks[i].n0, lks[i].n1, make_weight(), lks[i].syn_type);
	}
	for (size_t n0 = 0; n0 < nnum; n0 ++){
		is_inhi_nrn(n0)?sys->set_nrnType(n0,inhi):sys->set_nrnType(n0,exci);
	}
}

mod_init::Random::Random() :
	InitModule("random")
{
	using namespace prm::tag;
	param.add_var("prob", prob = 0.1)
		<< Name("prob")
		<< Desc("probability of connectivity")
		<< Range(0, 1)
		<< Step(0.001, 0.05)
		<< CmdLine()
		<< Save();
	param.add_var("bidir", bidir = false)
		<< Name("bidir")
		<< Desc("bi-directional connection")
		<< Control()
		<< CmdLine()
		<< Save();
}

void mod_init::SongSmallWorld::initialize(MEAsim * sys)
{
	rng.init(seed);
	place(sys);
	make_inhi_nrns(nnum);
	double Dmax = std::floor((double)nnum/2.0);
	//std::cerr << "Dmax =  " << Dmax << std::endl; 
	double p0 = (double)(knum*2)/(double)(nnum-1);
	//std::cerr << "p0 =  " << p0 << std::endl; 
	double d;
	double pij;
	std::vector<Lk> lks;
	if (bidir) for (int n0 = 0; n0 < nnum; n0 ++) {
		for (int n1 = n0 + 1; n1 < nnum; n1 ++){
			//std::cerr << "n0 =  " << n0 << std::endl; 
			//std::cerr << "n1 =  " << n1 << std::endl;
			//std::cerr << "n0-n1 =  " << n0-n1 << std::endl;
			 
			//std::cerr << "abs(n0-n1) =  " << std::abs((double)(n0-n1)) << std::endl; 
			//std::cerr << "nnum-abs(n0-n1) =  " << nnum-std::abs((double)(n0-n1)) << std::endl; 

			if ((std::abs(n0-n1)) < (nnum-std::abs(n0-n1)))
				d = (double)std::abs(n0-n1)/Dmax;
			else
				d = (double)(nnum-std::abs(n0-n1))/Dmax;
			if ((p0-d)>=0)
				pij = rewiredP*p0 + (1-rewiredP);
			else
				pij = rewiredP*p0;
			//std::cerr << "d =  " << d << std::endl; 
			//std::cerr << "pij =  " << pij << std::endl; 

			if (rng.uniform() < pij) {
				//lks.push_back(Lk(n0, n1));
				//lks.push_back(Lk(n1, n0));
				is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
				is_inhi_nrn(n1)?lks.push_back(Lk(n1, n0,inhi)):lks.push_back(Lk(n1, n0, exci));
			}
		}
	}
	else for (int n0 = 0; n0 < nnum; n0 ++) {
		for (int n1 = n0 + 1; n1 < nnum; n1 ++) {
			//std::cerr << "n0 =  " << n0 << std::endl; 
			//std::cerr << "n1 =  " << n1 << std::endl;
			//std::cerr << "n0-n1 =  " << n0-n1 << std::endl;
			//std::cerr << "abs(n0-n1) =  " << std::abs((double)(n0-n1)) << std::endl; 
			//std::cerr << "nnum-abs(n0-n1) =  " << nnum-std::abs((double)(n0-n1)) << std::endl; 

			if ((std::abs(n0-n1)) < (nnum-std::abs(n0-n1)))
				d = (double)std::abs(n0-n1)/Dmax;
			else
				d = (double)(nnum-std::abs(n0-n1))/Dmax;
			if ((p0-d)>=0)
				pij = rewiredP*p0 + (1-rewiredP);
			else
				pij = rewiredP*p0;
			//std::cerr << "d =  " << d << std::endl; 
			//std::cerr << "pij =  " << pij << std::endl; 
		
			//if (rng.uniform() < prob) lks.push_back(Lk(n0, n1));
			//if (rng.uniform() < prob) lks.push_back(Lk(n1, n0));
			if (rng.uniform() < pij) is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
			if (rng.uniform() < pij) is_inhi_nrn(n1)?lks.push_back(Lk(n1, n0,inhi)):lks.push_back(Lk(n1, n0, exci));
		}
	}
	// std::cerr << "creating " << lks.size() << " links.\n";
	sys->set_snum(lks.size());
	for (size_t i = 0; i < lks.size(); i ++) {
		sys->set_link(i, lks[i].n0, lks[i].n1, make_weight(), lks[i].syn_type);
	}
	for (size_t n0 = 0; n0 < nnum; n0 ++){
		is_inhi_nrn(n0)?sys->set_nrnType(n0,inhi):sys->set_nrnType(n0,exci);
	}
}

mod_init::SongSmallWorld::SongSmallWorld() :
	InitModule("songsw")
{
	using namespace prm::tag;
	param.add_var("knum", knum = 4)
		<< Name("knum")
		<< Desc("number of connected neighboring neurons")
		<< Range(1, 4095)
		<< Step(1, 20)
		<< CmdLine()
		<< Save();
	param.add_var("rewiredP", rewiredP = 0.0)
		<< Name("rewiredP")
		<< Desc("rewired probability")
		<< Range(0, 1)
		<< Step(0.001, 0.05)
		<< CmdLine()
		<< Save();
	param.add_var("bidir", bidir = false)
		<< Name("bidir")
		<< Desc("bi-directional connection")
		<< Control()
		<< CmdLine()
		<< Save();
}

namespace {
	struct DS
	{
		double r2;
		size_t nn;
	};
	bool cds(DS a, DS b) {return a.r2 < b.r2;}
}

void mod_init::Minimal::initialize(MEAsim * sys)
{
	rng.init(seed);
	place(sys);
	std::vector<DS> dss;
	for (size_t n0 = 0; n0 < nnum; n0 ++)
	for (size_t n1 = n0 + 1; n1 < nnum; n1 ++) {
		double dx = sys->get_node(n0).x - sys->get_node(n1).x;
		double dy = sys->get_node(n0).y - sys->get_node(n1).y;
		double r2 = dx * dx + dy * dy;
		DS d = {r2, n0 * nnum + n1};
		dss.push_back(d);
	}
	std::sort(dss.begin(), dss.end(), cds);
	std::vector<size_t> clr(nnum);
	for (size_t n = 0; n < nnum; n ++) clr[n] = n;
	std::vector<DS>::iterator i = dss.begin();
	std::vector<Lk> lks;
	for (size_t nclr = nnum; nclr > 1; i ++) {
		lks.push_back(Lk(i->nn / nnum, i->nn % nnum));
		lks.push_back(Lk(i->nn % nnum, i->nn / nnum));
		size_t c0 = clr[i->nn / nnum];
		size_t c1 = clr[i->nn % nnum];
		if (c0 == c1) continue;
		for (size_t n = 0; n < nnum; n ++) if (clr[n] == c1) clr[n] = c0;
		nclr --;
	}
	size_t rem = dss.size() - lks.size();
	for (size_t ov = 0; ov < over * rem; ov ++) {
		lks.push_back(Lk(i->nn / nnum, i->nn % nnum));
		lks.push_back(Lk(i->nn % nnum, i->nn / nnum));
		i ++;
	}
	sys->set_snum(lks.size());
	for (size_t i = 0; i < lks.size(); i ++) {
		sys->set_link(i, lks[i].n0, lks[i].n1, make_weight());
	}
}

mod_init::Minimal::Minimal() :
	InitModule("minimal")
{
	using namespace prm::tag;
	param.add_var("over", over = 0)
		<< Name("over")
		<< Desc("over connect fraction")
		<< Range(0, 1)
		<< Step(0.01, 0.1)
		<< CmdLine()
		<< Save();
}


void mod_init::MoritaSong::initialize(MEAsim * sys)
{
	rng.init(seed);
	sys->set_nnum(nnum);
	make_inhi_nrns(nnum);
	double theta;
	if (gamma > 2) theta = pow((gamma - 2) / (gamma - 1), 2) * al / nnum;
	else theta = pow(al * (2 - gamma) / (nnum * log(nnum)), 1 / (gamma - 1));
	std::vector<Pz> pzs;
	std::vector<double> a;
	double dmax = 0;
	for (size_t n = 0; n < nnum; n ++) {
		Pz p(rng.uniform(), rng.uniform());
		for (size_t m = 0; m < pzs.size(); m ++) {
			double d = p.r(pzs[m]);
			if (d > dmax) dmax = d;
		}
		pzs.push_back(p);
		a.push_back(pow((n + 1.0) / nnum, 1 / (1 - gamma)));
		sys->place(n, p.x, p.y);
	}
	std::vector<Lk> lks;
	for (size_t i = 0; i < nnum; i ++) for (size_t j = i + 1; j < nnum; j ++) {
		double d = pzs[i].r(pzs[j]);
		double q = pow(2 * d, delta) / (a[i] * a[j]);
		if (q < theta) {
			//lks.push_back(Lk(i, j));
			//lks.push_back(Lk(j, i));
			is_inhi_nrn(i)?lks.push_back(Lk(i, j,inhi)):lks.push_back(Lk(i, j, exci));
			is_inhi_nrn(j)?lks.push_back(Lk(j, i,inhi)):lks.push_back(Lk(j, i, exci));
		}
	}
	sys->set_snum(lks.size());
	for (size_t i = 0; i < lks.size(); i ++) {
		sys->set_link(i, lks[i].n0, lks[i].n1, make_weight(), lks[i].syn_type);
	}
	for (size_t n0 = 0; n0 < nnum; n0 ++){
		is_inhi_nrn(n0)?sys->set_nrnType(n0,inhi):sys->set_nrnType(n0,exci);
	}

}

mod_init::MoritaSong::MoritaSong() :
	InitModule("morita-song")
{
	using namespace prm::tag;
	weight = 1; // change default weight
	wdev = 0.2;
	param.add_var("nnum", nnum = 100)
		<< Name("nnum")
		<< Desc("number of neurons")
		<< Range(4, 16384)
		<< Step(1, 20)
		<< CmdLine()
		<< Save();	
	param.add_var("gamma", gamma = 3)
		<< Name("gamma")
		<< Desc("gamma parameter")
		<< Range(1.001, 10)
		<< Step(0.001, 0.1)
		<< CmdLine()
		<< Save();
	param.add_var("delta", delta = 4.7)
		<< Name("delta")
		<< Desc("d, dimension")
		<< Range(0.1, 10)
		<< Step(0.1, 1)
		<< CmdLine()
		<< Save();
	param.add_var("al", al = 8)
		<< Name("al")
		<< Desc("average number of linking per node")
		<< Range(1, 20)
		<< Step(0.1, 1)
		<< CmdLine()
		<< Save();
}

void mod_init::AdjFile::initialize(MEAsim * sys)
{
	std::ifstream f(file.c_str());
	std::vector<int> vs;
	int ns = 0;
	while (true) {
		int v;
		f >> v;
		if (! f.good()) break;
		vs.push_back(v);
		if (v) ns ++;
	}
	size_t nn = 0;
	while (nn * nn < vs.size()) nn ++;
	if (nn * nn != vs.size()) return;
	sys->set_nnum(nn);
	for (size_t n = 0; n < nn; n ++) sys->place(n, rng.uniform(), rng.uniform());
	sys->set_snum(ns);
	size_t l = 0;
	for (size_t i = 0; i < nn * nn; i ++) if (vs[i]) sys->set_link(l ++, i / nn, i % nn, make_gaussian());
}

mod_init::AdjFile::AdjFile() :
	InitModule("adj-file")
{
	using namespace prm::tag;
	weight = 4; // change default weight
	wdev = 0.2;
	param.add_var("file", file)
		<< Name("file")
		<< Desc("file to read network from")
		<< Word("FILE")
		<< prm::tag::File()
		<< Control()
		<< CmdLine()
		<< Save();
}

namespace {
	struct Link
	{
		size_t n0;
		size_t n1;
	};
}

void mod_init::File::initialize(MEAsim * sys)
{
	H5File f(file);
	hid_t g = f.open_load();
	sys->h5_load(g);

	if (wiring == "") return;

	std::vector<Link> links;
	std::ifstream lf(wiring.c_str());
	while (! lf.eof()) {
		Link l;
		lf >> l.n0;
		lf >> l.n1;
		if (lf.good()) links.push_back(l);
	}
	size_t snum = sys->get_snum();
	if (links.size() != snum) throw;
	
	for (size_t l = 0; l < snum; l ++) {
		MEAsim::Link lk = sys->get_link(l);
		sys->set_link(l, links[l].n0, links[l].n1, lk.weight);
	}
}

mod_init::File::File() :
	InitModule("file")
{
	using namespace prm::tag;
	param.add_var("file", file)
		<< Name("file")
		<< Desc("file to read network from")
		<< Word("FILE")
		<< prm::tag::File()
		<< Control()
		<< CmdLine()
		<< Save();
	param.add_var("wiring", wiring = "")
		<< Name("wiring")
		<< Desc("wiring replacement file")
		<< Word("FILE")
		<< prm::tag::File()
		<< Control()
		<< CmdLine()
		<< Save();
}

void mod_init::Rewire::initialize(MEAsim * sys)
{
	File::initialize(sys);
	size_t snum = sys->get_snum();
	size_t nnum = sys->get_nnum();

	bool * rwd = new bool[snum];
	for (size_t l = 0; l < snum; l ++) rwd[l] = false;

	size_t rwn = 0;
	if (keepk) while (rwn + 1 < snum * frac) {
		size_t lk1 = rng.uniform() * snum;
		size_t lk2 = rng.uniform() * snum;
		MEAsim::Link l1 = sys->get_link(lk1);
		MEAsim::Link l2 = sys->get_link(lk2);
		if (l1.from == l2.from || l1.from == l2.to || l1.to == l2.from || l1.to == l2.to) continue;
		size_t lk1r = snum;
		size_t lk2r = snum;
		double l1rw = 0;
		double l2rw = 0;
		for (size_t lk = 0; lk < snum; lk ++) {
			MEAsim::Link const & l = sys->get_link(lk);
			if (l.from == l1.to && l.to == l1.from) {
				lk1r = lk;
				l1rw = l.weight;
			}
			if (l.from == l2.to && l.to == l2.from) {
				lk2r = lk;
				l2rw = l.weight;
			}
		}
		sys->set_link(lk1, l1.from, l2.to, l1.weight);
		if (! rwd[lk1]) {
			rwd[lk1] = true;
			rwn ++;
		}
		sys->set_link(lk2, l2.from, l1.to, l2.weight);
		if (! rwd[lk2]) {
			rwd[lk2] = true;
			rwn ++;
		}
		if (lk1r < snum) {
			sys->set_link(lk1r, l2.to, l1.from, l1rw);
			if (! rwd[lk1r]) {
				rwd[lk1r] = true;
				rwn ++;
			}
		}
		if (lk2r < snum) {
			sys->set_link(lk2r, l1.to, l2.from, l2rw);
			if (! rwd[lk2r]) {
				rwd[lk1r] = true;
				rwn ++;
			}
		}
	}
	else while (rwn + 1 < snum * frac) {
		size_t lk = rng.uniform() * snum;
		MEAsim::Link l = sys->get_link(lk);
		size_t lkr = snum;
		double lrw = 0;
		size_t lton = 0;
		bool dup = false;
		size_t nto = rng.uniform() * nnum;
		for (size_t i = 0; i < snum; i ++) {
			MEAsim::Link const & ll = sys->get_link(i);
			if (ll.from == l.to && ll.to == l.from) {
				lkr = i;
				lrw = ll.weight;
			}
			if (ll.to == l.to) lton ++;
			if (ll.from == l.from && ll.to == nto) dup = true;
		}
		if (lton < 2 || dup == true) continue;
		sys->set_link(lk, l.from, nto, l.weight);
		if (! rwd[lk]) {
			rwd[lk] = true;
			rwn ++;
		}
		if (lkr < snum) {
			sys->set_link(lkr, nto, l.to, lrw);
			if (! rwd[lkr]) {
				rwd[lkr] = true;
				rwn ++;
			}
		}
	}

	delete [] rwd;
}

mod_init::Rewire::Rewire() :
	InitModule("rewire")
{
	using namespace prm::tag;
	param.add_var("frac", frac = 0.1)
		<< Name("frac")
		<< Desc("fraction of rewired links")
		<< Range(0, 1)
		<< Step(0.001, 0.05)
		<< Control()
		<< CmdLine()
		<< Save();
	param.add_var("keepk", keepk = true)
		<< Name("keepk")
		<< Desc("keep degrees for nodes")
		<< Control()
		<< CmdLine()
		<< Save();
}

void mod_init::NK::initialize(MEAsim * sys)
{
	rng.init(seed);
	place(sys);

	sys->set_snum(nnum * degree);
	size_t sn = 0;
	for (size_t n = 0; n < nnum; n ++) {
		std::vector<size_t> n0s;
		for (size_t m = 0; m < nnum; m ++) if (m != n) n0s.push_back(m);
		for (size_t i = 0; i < degree; i ++) {
			size_t d = rng.uniform() * n0s.size();
			sys->set_link(sn, n0s[d], n, make_weight());
			sn ++;
			n0s.erase(n0s.begin() + d);
		}
	}

}

mod_init::NK::NK() :
	InitModule("nk")
{
	using namespace prm::tag;
	param.add_var("degree", degree = 8)
		<< Name("degree")
		<< Desc("input degree")
		<< Range(1, 100)
		<< Step(1, 10)
		<< CmdLine()
		<< Save();
}

double mod_init::AxonSearch::make_axlen()
{
	if (ldevf == 0) return 1/sqrt(nnum) * axlenf;
	double l;
	do l = 1/sqrt(nnum) * axlenf + ldevf * 1/sqrt(nnum) * axlenf * rng.gaussian();
	while (l < 0);
	return l;
}

bool mod_init::AxonSearch::iden_connection(double r, double l, double xs, double ys, double xt, double yt, double vex, double vey)
{
	double lnn = sqrt(pow(xt-xs,2)+pow(yt-ys,2)); // distance between two neurons
	double alpha1 = acos(((xt-xs)*vex+(yt-ys)*vey)/l/lnn);
	if (alpha1 >= pi/2)
	  if (lnn < r) return true;
	  else return false;
	else{
	        if (lnn * sin(alpha1) > r) return false;
	        else{
	               if (lnn * cos(alpha1) <= l) return true;
	               else if (sqrt(pow(lnn * cos(alpha1)-l,2)+pow(lnn * sin(alpha1),2)) < r) return true;
	        }
	}

}

void mod_init::AxonSearch::initialize(MEAsim * sys)
{
	rng.init(seed);
	place(sys);
	make_inhi_nrns(nnum);
	double denR = 1/sqrt(nnum) * dradf; // radius of dendrite region
	double axlen = make_axlen(); // axon length
	std::vector<Lk> lks;
	for (size_t n0 = 0; n0 < nnum; n0 ++) {
		double axlen = make_axlen(); // assign axon length
		double theta = make_growAng(); // assign the direction of axon
		MEAsim::Node node_inf = sys->get_node(n0);
		double xs = node_inf.x;
		double ys = node_inf.y;
		double vex = axlen * cos(theta);
		double vey = axlen * sin(theta);
		for (size_t n1 = 0; n1 < nnum; n1 ++) {
			if (n0==n1) continue;
			MEAsim::Node node_inf = sys->get_node(n1);
			double xt = node_inf.x;
			double yt = node_inf.y;
			if (iden_connection(denR,axlen,xs,ys,xt,yt,vex,vey)) is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
			if (xs+vex > 1 || xs+vex < 0 || ys+vey > 1 || ys+vey < 0){
			         if (xs+vex > 1) xt = 2 - xt;
			         if (xs+vex < 0) xt = -xt;
			         if (ys+vey > 1) yt = 2 - yt;
			         if (ys+vey < 0) yt = -yt;
			         if (iden_connection(denR,axlen,xs,ys,xt,yt,vex,vey)) is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
			}
		}
	}
	// std::cerr << "creating " << lks.size() << " links.\n";
	sys->set_snum(lks.size());
	for (size_t i = 0; i < lks.size(); i ++) {
		sys->set_link(i, lks[i].n0, lks[i].n1, make_weight(), lks[i].syn_type);
	}
	for (size_t n0 = 0; n0 < nnum; n0 ++){
		is_inhi_nrn(n0)?sys->set_nrnType(n0,inhi):sys->set_nrnType(n0,exci);
	}
}

mod_init::AxonSearch::AxonSearch() :
	InitModule("axonsearch")
{
	using namespace prm::tag;
	param.add_var("dradf", dradf = 0.5)
		<< Name("dradf")
		<< Desc("factor for radius of dendrite region")
		<< Range(0, 1)
		<< Step(0.001, 0.05)
		<< CmdLine()
		<< Save();
	param.add_var("axlenf", axlenf = 3)
		<< Name("axlenf")
		<< Desc("factor for mean length of axon")
	        << Range(1, 20)
		<< Step(0.01, 0.5)
		<< CmdLine()
		<< Save();
	param.add_var("ldevf", ldevf = 0)
		<< Name("ldevf")
		<< Desc("factor for axon length spread")
		<< Range(0, 0.3)
		<< Step(0.001, 0.05)
		<< CmdLine()
		<< Save();
}

void mod_init::DistaGau::initialize(MEAsim * sys)
{
	rng.init(seed);
	place(sys);
	make_inhi_nrns(nnum);
	std::vector<Lk> lks;
	if (bidir) for (size_t n0 = 0; n0 < nnum; n0 ++) {
		for (size_t n1 = n0 + 1; n1 < nnum; n1 ++){
			double xs = sys->get_node(n0).x;
			double ys = sys->get_node(n0).y;			
			double xt = sys->get_node(n1).x;
			double yt = sys->get_node(n1).y;
			double lnn = sqrt((xt-xs)*(xt-xs)+(yt-ys)*(yt-ys));
			double gau_std = adisf/sqrt(nnum); 
			if (abs(rng.gaussian()) > lnn/gau_std) {
				is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
				is_inhi_nrn(n1)?lks.push_back(Lk(n1, n0,inhi)):lks.push_back(Lk(n1, n0, exci));
			}
		}
	}
	else for (size_t n0 = 0; n0 < nnum; n0 ++) {
		for (size_t n1 = n0 + 1; n1 < nnum; n1 ++) {
			double xs = sys->get_node(n0).x;
			double ys = sys->get_node(n0).y;			
			double xt = sys->get_node(n1).x;
			double yt = sys->get_node(n1).y;
			double lnn = sqrt((xt-xs)*(xt-xs)+(yt-ys)*(yt-ys));
			double gau_std = adisf/sqrt(nnum); 
			if (abs(rng.gaussian()) > lnn/gau_std) is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
			if (abs(rng.gaussian()) > lnn/gau_std) is_inhi_nrn(n1)?lks.push_back(Lk(n1, n0,inhi)):lks.push_back(Lk(n1, n0, exci));
		}
	}
	// std::cerr << "creating " << lks.size() << " links.\n";
	sys->set_snum(lks.size());
	for (size_t i = 0; i < lks.size(); i ++) {
		sys->set_link(i, lks[i].n0, lks[i].n1, make_weight(), lks[i].syn_type);
	}
	for (size_t n0 = 0; n0 < nnum; n0 ++){
		is_inhi_nrn(n0)?sys->set_nrnType(n0,inhi):sys->set_nrnType(n0,exci);
	}
}
mod_init::DistaGau::DistaGau() :
	InitModule("distagau")
{
	using namespace prm::tag;
	param.add_var("adisf", adisf = 3)
		<< Name("adisf")
		<< Desc("factor for standard deviation of gaussian probability function")
	        << Range(0, 40)
		<< Step(0.01, 0.5)
		<< CmdLine()
		<< Save();
	param.add_var("bidir", bidir = false)
		<< Name("bidir")
		<< Desc("bi-directional connection")
		<< Control()
		<< CmdLine()
		<< Save();
}
void mod_init::DistExp::initialize(MEAsim * sys)
{
	rng.init(seed);
	place(sys);
	make_inhi_nrns(nnum);
	std::vector<Lk> lks;
	if (bidir) for (size_t n0 = 0; n0 < nnum; n0 ++) {
		for (size_t n1 = n0 + 1; n1 < nnum; n1 ++){
			double xs = sys->get_node(n0).x;
			double ys = sys->get_node(n0).y;			
			double xt = sys->get_node(n1).x;
			double yt = sys->get_node(n1).y;
			double lnn = sqrt((xt-xs)*(xt-xs)+(yt-ys)*(yt-ys));
			double adisc = adisf/sqrt(nnum); 
			if (rng.uniform() < exp(-lnn/adisc)) {
				is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
				is_inhi_nrn(n1)?lks.push_back(Lk(n1, n0,inhi)):lks.push_back(Lk(n1, n0, exci));
			}
		}
	}
	else for (size_t n0 = 0; n0 < nnum; n0 ++) {
		for (size_t n1 = n0 + 1; n1 < nnum; n1 ++) {
			double xs = sys->get_node(n0).x;
			double ys = sys->get_node(n0).y;			
			double xt = sys->get_node(n1).x;
			double yt = sys->get_node(n1).y;
			double lnn = sqrt((xt-xs)*(xt-xs)+(yt-ys)*(yt-ys));
			double adisc = adisf/sqrt(nnum); 
			if (rng.uniform() < exp(-lnn/adisc)) is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
			if (rng.uniform() < exp(-lnn/adisc)) is_inhi_nrn(n1)?lks.push_back(Lk(n1, n0,inhi)):lks.push_back(Lk(n1, n0, exci));
		}
	}
	// std::cerr << "creating " << lks.size() << " links.\n";
	sys->set_snum(lks.size());
	for (size_t i = 0; i < lks.size(); i ++) {
		sys->set_link(i, lks[i].n0, lks[i].n1, make_weight(), lks[i].syn_type);
	}
	for (size_t n0 = 0; n0 < nnum; n0 ++){
		is_inhi_nrn(n0)?sys->set_nrnType(n0,inhi):sys->set_nrnType(n0,exci);
	}
}
mod_init::DistExp::DistExp() :
	InitModule("distexp")
{
	using namespace prm::tag;
	param.add_var("adisf", adisf = 3)
		<< Name("adisf")
		<< Desc("factor for space constant of the exponential probability function")
	        << Range(0, 40)
		<< Step(0.01, 0.5)
		<< CmdLine()
		<< Save();
	param.add_var("bidir", bidir = false)
		<< Name("bidir")
		<< Desc("bi-directional connection")
		<< Control()
		<< CmdLine()
		<< Save();
}
void mod_init::DistConst::initialize(MEAsim * sys)
{
	rng.init(seed);
	place(sys);
	make_inhi_nrns(nnum);
	std::vector<Lk> lks;
	for (size_t n0 = 0; n0 < nnum; n0 ++) {
		for (size_t n1 = n0 + 1; n1 < nnum; n1 ++){
			double xs = sys->get_node(n0).x;
			double ys = sys->get_node(n0).y;			
			double xt = sys->get_node(n1).x;
			double yt = sys->get_node(n1).y;
			double lnn = sqrt((xt-xs)*(xt-xs)+(yt-ys)*(yt-ys));
			double adisc = adisf/sqrt(nnum); 
			if (lnn < adisc) {
				is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
				is_inhi_nrn(n1)?lks.push_back(Lk(n1, n0,inhi)):lks.push_back(Lk(n1, n0, exci));
			}
		}
	}

	// std::cerr << "creating " << lks.size() << " links.\n";
	sys->set_snum(lks.size());
	for (size_t i = 0; i < lks.size(); i ++) {
		sys->set_link(i, lks[i].n0, lks[i].n1, make_weight(), lks[i].syn_type);
	}
	for (size_t n0 = 0; n0 < nnum; n0 ++){
		is_inhi_nrn(n0)?sys->set_nrnType(n0,inhi):sys->set_nrnType(n0,exci);
	}
}
mod_init::DistConst::DistConst() :
	InitModule("distconst")
{
	using namespace prm::tag;
	param.add_var("adisf", adisf = 3)
		<< Name("adisf")
		<< Desc("factor for space constant")
	        << Range(0, 40)
		<< Step(0.01, 0.5)
		<< CmdLine()
		<< Save();
}


void mod_init::Shordist::initialize(MEAsim * sys)
{
	rng.init(seed);
	place(sys);
	make_inhi_nrns(nnum);
	std::vector<Lk> lks;

	for (size_t n0 = 0; n0 < nnum; n0 ++) {
		chosed_list.clear();
		double xs = sys->get_node(n0).x;
		double ys = sys->get_node(n0).y;			
		for (size_t sn = 0; sn < outdeg; sn ++){
			double short_dist = 2.0;
			size_t s = 0;
			for (size_t n1 = 0; n1 < nnum; n1 ++) {
				if (n0==n1 || isChosed(n1)) continue;
				double xt = sys->get_node(n1).x;
				double yt = sys->get_node(n1).y;
				double lnn = sqrt((xt-xs)*(xt-xs)+(yt-ys)*(yt-ys));
				if (lnn < short_dist) {short_dist = lnn; s = n1;}
			}
			chosed_list.push_back(s);
			is_inhi_nrn(n0)?lks.push_back(Lk(n0, s,inhi)):lks.push_back(Lk(n0, s, exci));
		}
	}
	// std::cerr << "creating " << lks.size() << " links.\n";
	sys->set_snum(lks.size());
	for (size_t i = 0; i < lks.size(); i ++) {
		sys->set_link(i, lks[i].n0, lks[i].n1, make_weight(), lks[i].syn_type);
	}
	for (size_t n0 = 0; n0 < nnum; n0 ++){
		is_inhi_nrn(n0)?sys->set_nrnType(n0,inhi):sys->set_nrnType(n0,exci);
	}
}
bool mod_init::Shordist::isChosed(size_t n)
{
    std::vector<size_t>::iterator biter = chosed_list.begin();
    while (biter != chosed_list.end()) {
        if (*biter == n) return true;
        ++biter;
    }

    return false;
}

mod_init::Shordist::Shordist() :
	InitModule("shordist")
{
	using namespace prm::tag;
	param.add_var("outdeg", outdeg = 4)
		<< Name("outdeg")
		<< Desc("number of out-degree connections")
	        << Range(4, 1024)
		<< Step(1, 4)
		<< CmdLine()
		<< Save();
	param.add_var("bidir", bidir = false)
		<< Name("bidir")
		<< Desc("bi-directional connection")
		<< Control()
		<< CmdLine()
		<< Save();
}

void mod_init::Distexpin::initialize(MEAsim * sys)
{
	rng.init(seed);
	place(sys);
	make_inhi_nrns(nnum);
	std::vector<Lk> lks;
	if (bidir) for (size_t n0 = 0; n0 < nnum; n0 ++) {
		for (size_t n1 = n0 + 1; n1 < nnum; n1 ++){
			double xs = sys->get_node(n0).x;
			double ys = sys->get_node(n0).y;			
			double xt = sys->get_node(n1).x;
			double yt = sys->get_node(n1).y;
			double lnn = sqrt((xt-xs)*(xt-xs)+(yt-ys)*(yt-ys));
			double adisc = adisf/sqrt(nnum);
			double maxdisc = maxdisf/sqrt(nnum);
			if ( lnn < maxdisc && rng.uniform() < exp(-lnn/adisc)) {
				is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
				is_inhi_nrn(n1)?lks.push_back(Lk(n1, n0,inhi)):lks.push_back(Lk(n1, n0, exci));
			}
		}
	}
	else for (size_t n0 = 0; n0 < nnum; n0 ++) {
		for (size_t n1 = n0 + 1; n1 < nnum; n1 ++) {
			double xs = sys->get_node(n0).x;
			double ys = sys->get_node(n0).y;			
			double xt = sys->get_node(n1).x;
			double yt = sys->get_node(n1).y;
			double lnn = sqrt((xt-xs)*(xt-xs)+(yt-ys)*(yt-ys));
			double adisc = adisf/sqrt(nnum); 
			double maxdisc = maxdisf/sqrt(nnum);
			if (lnn < maxdisc && rng.uniform() < exp(-lnn/adisc)) is_inhi_nrn(n0)?lks.push_back(Lk(n0, n1,inhi)):lks.push_back(Lk(n0, n1, exci));
			if (lnn < maxdisc && rng.uniform() < exp(-lnn/adisc)) is_inhi_nrn(n1)?lks.push_back(Lk(n1, n0,inhi)):lks.push_back(Lk(n1, n0, exci));
		}
	}
	// std::cerr << "creating " << lks.size() << " links.\n";
	sys->set_snum(lks.size());
	for (size_t i = 0; i < lks.size(); i ++) {
		sys->set_link(i, lks[i].n0, lks[i].n1, make_weight(), lks[i].syn_type);
	}
	for (size_t n0 = 0; n0 < nnum; n0 ++){
		is_inhi_nrn(n0)?sys->set_nrnType(n0,inhi):sys->set_nrnType(n0,exci);
	}
}
mod_init::Distexpin::Distexpin() :
	InitModule("distexpin")
{
	using namespace prm::tag;
	param.add_var("adisf", adisf = 3)
		<< Name("adisf")
		<< Desc("factor for space constant of the exponential probability function")
	        << Range(0, 40)
		<< Step(0.01, 0.5)
		<< CmdLine()
		<< Save();
	param.add_var("maxdisf", maxdisf = 6)
		<< Name("maxdisf")
		<< Desc("factor for max space range of possible connection")
	        << Range(0, 40)
		<< Step(0.01, 0.5)
		<< CmdLine()
		<< Save();
	param.add_var("bidir", bidir = false)
		<< Name("bidir")
		<< Desc("bi-directional connection")
		<< Control()
		<< CmdLine()
		<< Save();
}

void mod_init::SelfConn::initialize(MEAsim * sys)
{
	rng.init(seed);
	place(sys);
	std::vector<Lk> lks;
	if (bidir) for (size_t n0 = 0; n0 < nnum; n0 ++) {
		for (size_t n1 = n0 + 1; n1 < nnum; n1 ++)
		if (rng.uniform() < prob) {
			lks.push_back(Lk(n0, n1));
			lks.push_back(Lk(n1, n0));
		}
	}
	else for (size_t n0 = 0; n0 < nnum; n0 ++) {
		for (size_t n1 = n0 + 1; n1 < nnum; n1 ++) {
			if (rng.uniform() < prob) lks.push_back(Lk(n0, n1));
			if (rng.uniform() < prob) lks.push_back(Lk(n1, n0));
		}
	}
	for (size_t n0 = 0; n0 < nnum; n0 ++) {
		lks.push_back(Lk(n0, n0));
	}

	// std::cerr << "creating " << lks.size() << " links.\n";
	sys->set_snum(lks.size());
	for (size_t i = 0; i < lks.size(); i ++) {
		if (lks[i].n0 == lks[i].n1)	sys->set_link(i, lks[i].n0, lks[i].n1, wself);
		else sys->set_link(i, lks[i].n0, lks[i].n1, make_weight());
	}
}

mod_init::SelfConn::SelfConn() :
	InitModule("selfconn")
{
	using namespace prm::tag;
	param.add_var("prob", prob = 0.1)
		<< Name("prob")
		<< Desc("probability of connectivity")
		<< Range(0, 1)
		<< Step(0.001, 0.05)
		<< CmdLine()
		<< Save();
	param.add_var("bidir", bidir = false)
		<< Name("bidir")
		<< Desc("bi-directional connection")
		<< Control()
		<< CmdLine()
		<< Save();	
	param.add_var("self weight", wself = 1)
		<< Name("self weight")
		<< Desc("mean synaptic weight of self-connection")
		<< Range(0, 1000)
		<< Step(0.1, 1)
		<< CmdLine()
		<< Save();
}

void mod_init::SelfMoritaSong::initialize(MEAsim * sys)
{
	rng.init(seed);
	sys->set_nnum(nnum);
	double theta;
	if (gamma > 2) theta = pow((gamma - 2) / (gamma - 1), 2) * al / nnum;
	else theta = pow(al * (2 - gamma) / (nnum * log(nnum)), 1 / (gamma - 1));
	std::vector<Pz> pzs;
	std::vector<double> a;
	double dmax = 0;
	for (size_t n = 0; n < nnum; n ++) {
		Pz p(rng.uniform(), rng.uniform());
		for (size_t m = 0; m < pzs.size(); m ++) {
			double d = p.r(pzs[m]);
			if (d > dmax) dmax = d;
		}
		pzs.push_back(p);
		a.push_back(pow((n + 1.0) / nnum, 1 / (1 - gamma)));
		sys->place(n, p.x, p.y);
	}
	std::vector<Lk> lks;
	for (size_t i = 0; i < nnum; i ++) for (size_t j = i + 1; j < nnum; j ++) {
		double d = pzs[i].r(pzs[j]);
		double q = pow(2 * d, delta) / (a[i] * a[j]);
		if (q < theta) {
			lks.push_back(Lk(i, j));
			lks.push_back(Lk(j, i));
		}
	}
	for (size_t n0 = 0; n0 < nnum; n0 ++) {
		lks.push_back(Lk(n0, n0));
	}

	sys->set_snum(lks.size());
	for (size_t i = 0; i < lks.size(); i ++) {
		if (lks[i].n0 == lks[i].n1)	sys->set_link(i, lks[i].n0, lks[i].n1, wself);
		else sys->set_link(i, lks[i].n0, lks[i].n1,  make_gaussian());
	}
}

mod_init::SelfMoritaSong::SelfMoritaSong() :
	InitModule("self-morita-song")
{
	using namespace prm::tag;
	weight = 1; // change default weight
	wdev = 0.2;
	param.add_var("nnum", nnum = 100)
		<< Name("nnum")
		<< Desc("number of neurons")
		<< Range(4, 16384)
		<< Step(1, 20)
		<< CmdLine()
		<< Save();	
	param.add_var("gamma", gamma = 3)
		<< Name("gamma")
		<< Desc("gamma parameter")
		<< Range(1.001, 10)
		<< Step(0.001, 0.1)
		<< CmdLine()
		<< Save();
	param.add_var("delta", delta = 4.7)
		<< Name("delta")
		<< Desc("d, dimension")
		<< Range(0.1, 10)
		<< Step(0.1, 1)
		<< CmdLine()
		<< Save();
	param.add_var("al", al = 8)
		<< Name("al")
		<< Desc("average number of linking per node")
		<< Range(1, 20)
		<< Step(0.1, 1)
		<< CmdLine()
		<< Save();
	param.add_var("self weight", wself = 1)
		<< Name("self weight")
		<< Desc("mean synaptic weight of self-connection")
		<< Range(0, 1000)
		<< Step(0.1, 1)
		<< CmdLine()
		<< Save();
}
