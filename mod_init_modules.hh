#include "initia.hh"
#include <vector>
#include <cmath>
#define pi M_PI

namespace mod_init {
	class WDev : // with a weight distribution
		public virtual InitModule
	{
	protected:
		double weight;
		double wdev;
		double make_weight();
		double make_gaussian();
	public:
		WDev();
	};

	class AddInhiNrn:             // By Justin, 20141208
		public virtual InitModule
	{
	protected:
		double inhiratio;
		std::vector<bool> inhi_nrns;
		void make_inhi_nrns(int nnrns);
		bool is_inhi_nrn(int nrn_idx);
	public:
		AddInhiNrn();
	};

	class Lattice :
		public virtual WDev
	{
		void initialize(MEAsim * sys) override;
		// parameters
		size_t lsz;
		double range;
		double rewiredP;
	public:
		Lattice();
	};

	class Spaced : // randomly spaced reurons
		public virtual InitModule
	{
	protected:
		// parameters
		size_t nnum;
		double mdisf;

		void place(MEAsim * sys);
	public:
		Spaced();
	};

	class Random :
		public virtual Spaced,
		public virtual WDev,
		public virtual AddInhiNrn	
	{
		void initialize(MEAsim * sys) override;
		// parameter
		double prob;
		bool bidir;
	public:
		Random();
	};
	class SongSmallWorld :
		public virtual Spaced,
		public virtual WDev,
		public virtual AddInhiNrn	
	{
		void initialize(MEAsim * sys) override;
		// parameter
		double rewiredP;
		size_t knum;
		bool bidir;
	public:
		SongSmallWorld();
	};
	class SelfConn :
		public virtual Spaced,
		public virtual WDev	
		{
		void initialize(MEAsim * sys) override;
		// parameter
		double prob;
		bool bidir;
		double wself;
	public:
		SelfConn();
	};
	class Minimal :
		public virtual Spaced,
		public virtual WDev
	{
		void initialize(MEAsim * sys) override;
		double over;
	public:
		Minimal();
	};

	class MoritaSong :
		public virtual WDev,
		public virtual AddInhiNrn
	{
		void initialize(MEAsim * sys) override;
		size_t nnum;
		double gamma;
		double delta;
		double al;
	public:
		MoritaSong();
	};
	class SelfMoritaSong :
		public virtual WDev
	{
		void initialize(MEAsim * sys) override;
		size_t nnum;
		double gamma;
		double delta;
		double al;
		double wself;		
	public:
		SelfMoritaSong();
	};
	class AdjFile :
		public virtual WDev
	{
		void initialize(MEAsim * sys) override;
		std::string file;
	public:
		AdjFile();
	};

	class File : // straight from file, allow command line parameter change
		public virtual InitModule
	{
	protected:
		void initialize(MEAsim * sys) override;
		std::string file;
		std::string wiring;
	public:
		File();
	};

	class Rewire : // load file and rewire the network
		public virtual File
	{
		void initialize(MEAsim * sys) override;
		double frac;
		bool keepk; // whether to keep degree for nodes
	public:
		Rewire();
	};

	class NK : // NK network
		public virtual Spaced,
		public virtual WDev
	{
		void initialize(MEAsim * sys) override;
		// parameter
		size_t degree;
	public:
		NK();
	};
	class AxonSearch : // RL axon searching network, Justin, 20141201
		public virtual Spaced,
		public virtual WDev,
		public virtual AddInhiNrn
	{
		void initialize(MEAsim * sys) override;
		double make_axlen();
		double make_growAng(){return rng.uniform() * 2 * pi;}
	        virtual bool iden_connection(double r, double l, double xs, double ys, double xt, double yt, double vex, double vey);
		// parameters
		double dradf; // factor for radius of dendrite region, radius = dradf times 1/sqrt(nnum)
		double axlenf; // factor for mean length of axon, length = axlenf times 1/sqrt(nnum)
		double ldevf; // factor for axon length spread, dev = ldevf times axon's length 
	public:
		AxonSearch();
	};
	class DistaGau : // distance-dependent connection with a Gaussian probability function, Justin, 20141204
		public virtual Spaced,
		public virtual WDev,
		public virtual AddInhiNrn
	{
		void initialize(MEAsim * sys) override;
		// parameters
		double adisf; // factor for standard deviation of the Gaussian probability function, std= adisf times 1/sqrt(nnum)
		bool bidir;
	public:
		 DistaGau();
	};
	class DistExp: // distance-dependent connection with a exponential probability function, Justin, 20141204
		public virtual Spaced,
		public virtual WDev,
		public virtual AddInhiNrn
	{
		void initialize(MEAsim * sys) override;
		// parameters
		double adisf; // factor for space constant of the exponential probability function, adisc= adisf times 1/sqrt(nnum)
		bool bidir;
	public:
		DistExp();
	};
	class DistConst: // distance-dependent connection with a constant distance, Justin, 20151001
		public virtual Spaced,
		public virtual WDev,
		public virtual AddInhiNrn
	{
		void initialize(MEAsim * sys) override;
		// parameters
		double adisf; // factor for space constant, dist = adisf times 1/sqrt(nnum)
	public:
		DistConst();
	};
	class Shordist: // top indeg other neurons with shortest distance from it are connected, Justin, 20150121
		public virtual Spaced,
		public virtual WDev,
		public virtual AddInhiNrn
	{
		void initialize(MEAsim * sys) override;
		bool isChosed(size_t n);
		std::vector<size_t> chosed_list;
		// parameters
		size_t outdeg; // the number of in-degree for each neuron 
		bool bidir;
	public:
		Shordist();
	};
	class Distexpin: // distance-dependent connection with a exponential probability function, Justin, 20141204
		public virtual Spaced,
		public virtual WDev,
		public virtual AddInhiNrn
	{
		void initialize(MEAsim * sys) override;
		// parameters
		double adisf; // factor for space constant of the exponential probability function, adisc= adisf times 1/sqrt(nnum)
		double maxdisf; // for max space range of possible connection, maxdisc= maxdisf times 1/sqrt(nnum)
		bool bidir;
	public:
		Distexpin();
	};

}
