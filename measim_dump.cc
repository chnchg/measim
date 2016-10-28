#include "measim_dump.hh"

void MEAsimDump::dash()
{
	double a0 = floor(get_age() / resolution);
	MEAsimRun::dash();
	double a1 = floor(get_age() / resolution);
	if (a0 == a1) return;
	std::cout << get_age(); 

	size_t nn = get_nnum();
	double aV = 0;
	double aW = 0;
	double aCa2p_r = 0;
	for (size_t i = 0; i < nn; i ++) {
		Neuron const & n = get_neuron(i);
		aV += n.V;
		aW += n.W;
		aCa2p_r += n.Ca2p_r;
	}
	aV /= nn;
	aW /= nn;
	aCa2p_r /= nn;
	std::cout << '\t' << aV;
	std::cout << '\t' << aW;
	std::cout << '\t' << aCa2p_r;

	size_t sn = get_snum();
	double aY = 0;
	double aZ = 0;
	double aS = 0;
	for (size_t i = 0; i < sn; i ++) {
		Synapse const & s = get_synapse(i);
		aY += s.Y;
		aZ += s.Z;
		aS += s.S;
	}
	aY /= sn;
	aZ /= sn;
	aS /= sn;
	std::cout << '\t' << aY;
	std::cout << '\t' << aZ;
	std::cout << '\t' << aS;

	std::cout << '\n';
}

MEAsimDump::MEAsimDump()
{
	description = "Dump average values of various variables:\n\n"
		"\tage  V  W  Ca2_p  Y  Z  S";
	using namespace prm::tag;
	param.add_var("resolution", resolution = 1)
		<< Name("resolution")
		<< Desc("time resolution")
		<< CmdLine();
}
