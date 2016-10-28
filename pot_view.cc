#include "pot_view.hh"
#include "net_view.hh"
#include "view_misc.hh"
#include <cmath>

struct CTrans
{
	double lower;
	double min;
	double ratio;
	bool has_zero;
	CTrans(double min, double max, double lower, double upper) :
		lower(lower), min(min),
		ratio((upper - lower) / (max - min)),
		has_zero(min * max <= 0)
	{}
	double operator()(double v) const
	{
		return lower + (v - min) * ratio;
	}
};

void PotView::on_draw()
{
	double const vmin = - 33;
	double const vmax = 17;
	double const wmin = 0.14;
	double const wmax = 0.54;

	static GLuint circ = make_circle(8);
	double ww = get_width();
	double hh = get_height();

	CTrans xt(wmin, wmax, 0, ww);
	CTrans yt(vmin, vmax, hh, 0);

	struct Point
	{
		double x;
		double y;
	} * pt = new Point [vol->get_sz()];
	for (size_t n = 0; n < vol->get_sz(); n ++) {
		double v = vol->get_pot(n);
		double w = vol->get_rst(n);
		pt[n].x = xt(w);
		pt[n].y = yt(v);
	}

	glClearColor(1, 1, 1, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// axis at zeros
	glColor3f(0.6, 0.6, 0.6);
	glBegin(GL_LINES);
	if (xt.has_zero) {
		glVertex2f(xt(0), 0);
		glVertex2f(xt(0), hh);
	}
	if (yt.has_zero) {
		glVertex2f(0, yt(0));
		glVertex2f(ww, yt(0));
	}
	glEnd();

	// x ticks
	size_t tticks = 0; // total number of ticks
	bool divf = true; // current scale divided by div10?
	double div10 = 1000;
	double tlen = 8;
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
	while (tticks < 8) {
		int b = ceil(divf ? wmin / div10 : wmin * div10);
		while (b <= (divf ? wmax / div10 : wmax * div10)) {
			double x = xt(divf ? b * div10 : b / div10);
			glVertex2f(x, 0);
			glVertex2f(x, tlen);
			glVertex2f(x, hh - tlen);
			glVertex2f(x, hh);
			tticks ++;
			b ++;
		}
		if (tticks) tlen /= 2;
		if (divf) {
			if (div10 >= 10) {
				div10 /= 10;
				continue;
			}
			divf = false;
		}
		div10 *= 10;
	}
	glEnd();

	// y ticks
	tticks = 0;
	divf = true;
	div10 = 1000;
	tlen = 8;
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
	while (tticks < 8) {
		int b = ceil(divf ? vmin / div10 : vmin * div10);
		while (b <= (divf ? vmax / div10 : vmax * div10)) {
			double y = yt(divf ? b * div10 : b / div10);
			glVertex2f(0, y);
			glVertex2f(tlen, y);
			glVertex2f(ww - tlen, y);
			glVertex2f(ww, y);
			tticks ++;
			b ++;
		}
		if (tticks > 2) tlen /= 2;
		if (divf) {
			if (div10 >= 10) {
				div10 /= 10;
				continue;
			}
			divf = false;
		}
		div10 *= 10;
	}
	glEnd();

	// threshold
	glColor3f(0, 0.7, 0);
	double thl = yt(10); // threshold == 10 mv
	glBegin(GL_LINES);
	glVertex2f(0, thl);
	glVertex2f(ww, thl);
	glEnd();

	if (mwlist.size()) {
		glBegin(GL_LINES);
		for (std::vector<ModW>::iterator i = mwlist.begin(); i != mwlist.end(); i ++) {
			size_t f = vol->get_from(i->l);
			size_t t = vol->get_to(i->l);
			double w = i->w;
			glColor3f(1 - w / 2, 1 - w / 8, 1 - w / 8);
			glVertex2f(pt[f].x, pt[f].y);
			glVertex2f(pt[t].x, pt[t].y);
		}
		glEnd();
	}
	else {
		glColor3f(0, 1, 1);
		glBegin(GL_LINES);
		for (size_t l = 0; l < vol->get_ln(); l ++) {
			size_t f = vol->get_from(l);
			size_t t = vol->get_to(l);
			glVertex2f(pt[f].x, pt[f].y);
			glVertex2f(pt[t].x, pt[t].y);
		}
		glEnd();
	}

	double r = 0.003 * (ww + hh);
	for (size_t n = 0; n < vol->get_sz(); n ++) {
		double c = vol->get_cal(n);
		glColor3f(c, 0, 1 - c);
		glPushMatrix();
		glTranslatef(pt[n].x, pt[n].y, 0);
		glScalef(r, r, r);
		glCallList(circ);
		glPopMatrix();
	}
	delete [] pt;
}
