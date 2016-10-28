#ifndef VIEW_MISC_HH
#define VIEW_MISC_HH
#include "measim.hh"
#include <gltk.hh>
#include <vector>
extern double const pi;
GLuint make_circle(size_t n = 16);
struct ModW
{
	size_t l;
	double w;
	ModW(size_t l, double w) : l(l), w(w) {}
};
extern std::vector<ModW> mwlist;
void make_mwlist(MEAsim const * vol);
#endif
