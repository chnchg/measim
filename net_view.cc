#include "net_view.hh"
#include "view_misc.hh"

void NetView::on_draw()
{
	static GLuint circ = make_circle(10);
	double ww = get_width();
	double bw = ww / 200;
	ww -= 3 * bw;
	double hh = get_height();
	size_t sz = vol->get_nnum();
	double r = 0.07 * (ww + hh) / sqrt(sz);
	struct Point
	{
		double x;
		double y;
	} * pt = new Point [sz];
	for (size_t n = 0; n < sz; n ++) {
		pt[n].x = r + (ww - 2 * r) * vol->get_node(n).x;
		pt[n].y = r + (hh - 2 * r) * vol->get_node(n).y;
	}
	make_mwlist(vol);
	glClearColor(1, 1, 1, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// activity trail
	double c = trail.size() / 100;
	for (std::deque<NetView::TrailPoint>::iterator i = trail.begin(); i != trail.end(); i ++) {
		double k = 0.6 + 0.4 * c;
		glColor3f(k, k, c);
		glPushMatrix();
		glTranslatef(i->x, i->y, 0);
		glScalef(i->r, i->r, i->r);
		glCallList(circ);
		glPopMatrix();
		c -= 0.01;
	}

	// links
	glBegin(GL_LINES);
	for (std::vector<ModW>::iterator i = mwlist.begin(); i != mwlist.end(); i ++) {
		size_t f = vol->get_link(i->l).from;
		size_t t = vol->get_link(i->l).to;
		double w = i->w;
		// double s = vol->get_synapse(i->l).S * 2;
		if (vol->get_link(i->l).syn_type == exci) glColor3f(1 - w, 1 - w / 4, 1 - w);
		else glColor3f(1, 1 - w / 4, 1 - w);
		glVertex2f(pt[f].x, pt[f].y);
		glVertex2f(pt[t].x, pt[t].y);
	}
	glEnd();

	// neurons
	double ax = 0;
	double ay = 0;
	double at = 0;
	for (size_t n = 0; n < sz; n ++) {
		MEAsim::Neuron const & nn = vol->get_neuron(n);
		double c = nn.Ca2p_r;
		double a = vol->active_over(n);
		at += a;
		ax += pt[n].x * a;
		ay += pt[n].y * a;
		glColor3f(c, 0, 1 - c);
		glPushMatrix();
		glTranslatef(pt[n].x, pt[n].y, 0);
		glScalef(r, r, r);
		glCallList(circ);
		glPopMatrix();
	}
	ax /= at;
	ay /= at;
	NetView::TrailPoint p = {ax, ay, 10 * r * at / sz};
	trail.push_back(p);
	while (trail.size() > 100) trail.pop_front();

	double tts = 0;
	double ttx = 0;
	double tty = 0;
	glBegin(GL_TRIANGLES);
	for (std::vector<ModW>::iterator i = mwlist.begin(); i != mwlist.end(); i ++) {
		size_t f = vol->get_link(i->l).from;
		size_t t = vol->get_link(i->l).to;
		// double w = i->w;
		double dx = pt[t].x - pt[f].x;
		double dy = pt[t].y - pt[f].y;
		double rr = sqrt(dx * dx + dy * dy);
		double ex = dx / rr;
		double ey = dy / rr;
		MEAsim::Synapse const & sn = vol->get_synapse(i->l);
		double s = sn.S;
		tts += s;
		tty += sn.Y;
		ttx += 1.0 - s - sn.Y - sn.Z;
		// glColor3f(s, 1, 1 - s);
		if (vol->get_link(i->l).syn_type == exci) glColor3f(1 - s, 1 - s / 4, 1 - s);
		else glColor3f(1, 1 - s / 4, 1 - s);

		glVertex2f(pt[t].x - r * ex, pt[t].y - r * ey);
		glVertex2f(pt[t].x - r * (4 * ex + ey) / 2, pt[t].y - r * (4 * ey - ex) / 2);
		glVertex2f(pt[t].x - r * (4 * ex - ey) / 2, pt[t].y - r * (4 * ey + ex) / 2);
	}
	glEnd();
	tts /= mwlist.size();
	tty /= mwlist.size();
	ttx /= mwlist.size();

	glColor3f(tts, 1, 1 - tts);
	glBegin(GL_QUADS);
	glVertex2f(ww, hh);
	glVertex2f(ww, hh * (1 - tts));
	glVertex2f(ww + bw, hh * (1 - tts));
	glVertex2f(ww + bw, hh);
	glEnd();

	ww += bw;
	glColor3f(1.0, 0.1, 0.1);
	glBegin(GL_QUADS);
	glVertex2f(ww, hh);
	glVertex2f(ww, hh * (1 - tty));
	glVertex2f(ww + bw, hh * (1 - tty));
	glVertex2f(ww + bw, hh);
	glEnd();

	ww += bw;
	glColor3f(0.5, 0.5, 1);
	glBegin(GL_QUADS);
	glVertex2f(ww, hh);
	glVertex2f(ww, hh * (1 - ttx));
	glVertex2f(ww + bw, hh * (1 - ttx));
	glVertex2f(ww + bw, hh);
	glEnd();

	delete [] pt;
}
