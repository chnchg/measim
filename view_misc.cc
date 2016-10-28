#include "view_misc.hh"
#include <algorithm>

double const pi = 3.141592653589793238462;

GLuint make_circle(size_t n)
{
	GLuint circ = glGenLists(1);
	glNewList(circ, GL_COMPILE);
	glBegin(GL_POLYGON);
	for (size_t i = 0; i < n; i ++) {
		double t = i * 2 * pi / n;
		glVertex2f(sin(t), cos(t));
	}
	glEnd();
	glEndList();
	return circ;
}

namespace {
	bool clw(ModW a, ModW b) {return a.w < b.w;}
}

std::vector<ModW> mwlist; // sorted weight list

void make_mwlist(MEAsim const * vol) {
	if (mwlist.size()) return;
	double mxw = 0;
	size_t ln = vol->get_snum();
	for (size_t l = 0; l < ln; l ++) {
		double w = vol->get_wmod(l);
		mwlist.push_back(ModW(l, w));
		if (w > mxw) mxw = w;
	}
	std::sort(mwlist.begin(), mwlist.end(), clw);
	for (size_t l = 0; l < ln; l ++) mwlist[l].w /= mxw;
}
