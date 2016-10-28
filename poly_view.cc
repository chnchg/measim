#include "poly_view.hh"
#include <fstream>
void PolyView::on_draw()
{
	double ww = get_width();
	double hh = get_height();
	size_t sz = vol->get_nnum();
	double start = vol->get_age() - range;
	double dh = 0.5 * hh / sz;

	glClearColor(1, 1, 1, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glColor3f(0, 0, 0);
	glBegin(GL_LINES);
	glVertex2f(0, hh / 2);
	glVertex2f(ww, hh / 2);
	glEnd();

	glColor3f(1, 0, 0);
	glBegin(GL_LINES);
	for (std::deque<Event>::iterator i = events.begin(); i != events.end(); i ++) {
		if (i->i >= sz) continue;
		double x = (i->time - start) * ww / range;
		double y = (i->i * dh);
		glVertex2f(x, y);
		glVertex2f(x, y + dh);
	}
	glEnd();

	glColor3f(0.5, 0.5, 0);
	glBegin(GL_LINE_STRIP);
	double tm = start; // window end
	double wd = 5.0; // window size
	size_t np = 0; // points in window
	double hs = 0.5 * hh / sz;
	std::deque<Event>::iterator b = events.begin();
	std::deque<Event>::iterator c = b;
	glVertex2f(0, hh);
	double tp = tm;
	size_t pp = 0;
	while (c != events.end()) {
		double pn = np;
		if (b->time <= c->time - wd) {
			tm = b->time + wd;
			b ++;
			np --;
		}
		else {
			tm = c->time;
			c ++;
			np ++;
		}
		if ((tm - tp) * ww / range > 1) {
			if (pn != pp) glVertex2f((tp - start) * ww / range, hh - pn * hs);
			glVertex2f((tm - start) * ww / range, hh - pn * hs);
			glVertex2f((tm - start) * ww / range, hh - np * hs);
			tp = tm;
			pp = np;
		}
	}
	while (b != events.end() && start + range >= b->time + wd) {
		double pn = np;
		tm = b->time + wd;
		b ++;
		np --;
		if ((tm - tp) * ww / range > 1) {
			if (pn != pp) glVertex2f((tp - start) * ww / range, hh - pn * hs);
			glVertex2f((tm - start) * ww / range, hh - pn * hs);
			glVertex2f((tm - start) * ww / range, hh - np * hs);
			tp = tm;
			pp = np;
		}
	}
	tm = start + range;
	if (pp != np) glVertex2f((tp - start) * ww / range, hh - np * hs);
	glVertex2f((tm - start) * ww / range, hh - np * hs);
	glEnd();
	glColor3f(0.3, 0.3, 0);
	std::ostringstream os;
	os << "range: " << range << "ms";
	glRasterPos2i(4, 14);
	for (size_t i = 0; i < os.str().size(); i ++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, os.str()[i]);
	os.str("");
	os.precision(1);
	os << std::fixed << "age: " << start + range << "ms";
	glRasterPos2i(4, 28);
	for (size_t i = 0; i < os.str().size(); i ++) glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, os.str()[i]);
}

void PolyView::on_menu_select(int item)
{
	using namespace gltk;
	if (item == 0) { // change range
		if (! rdlg) {
			rdlg = new Window;
			Box * bx = new Box(ORIENTATION_VERTICAL);
			rdlg->add(bx);
			Box * hb = new Box(ORIENTATION_HORIZONTAL);
			bx->pack_start(hb);
			hb->pack_start(new Label("Range:"));
			hb->pack_start(range_entry);
			hb->show();
			hb = new Box(ORIENTATION_HORIZONTAL);
			bx->pack_start(hb);
			Button * b = new Button("Cancel");
			hb->pack_start(b);
			b->signal_clicked().connect(sigc::mem_fun(rdlg, & Window::hide));
			b = new Button("OK");
			hb->pack_start(b);
			b->signal_clicked().connect(sigc::mem_fun(this, & PolyView::range_ok));
		}
		range_entry.set_value(range);
		rdlg->show();
	}
	else if (item == 1) { // save data
		if (! fdlg) {
			fdlg = new Window;
			Box * bx = new Box(ORIENTATION_VERTICAL);
			fdlg->add(bx);
			FileBrowser * bsr = new FileBrowser;
			bx->pack_start(bsr);
			bsr->signal_selected().connect(sigc::mem_fun(this, & PolyView::file_select));
			Box * hb = new Box(ORIENTATION_HORIZONTAL);
			bx->pack_start(hb);
			hb->pack_start(new Label("Filename:"));
			hb->pack_start(fn_entry);
			Button * b = new Button("Cancel");
			bx->pack_start(b);
			b->signal_clicked().connect(sigc::mem_fun(fdlg, & Window::hide));
			b = new Button("OK");
			bx->pack_start(b);
			b->signal_clicked().connect(sigc::mem_fun(this, & PolyView::file_ok));
		}
		fdlg->show();
	}
}

void PolyView::add(size_t i)
{
	double t = vol->get_age();
	if (events.size() && events.back().time > t) events.clear();
	events.push_back(Event(i, t));
	double start = t - range;
	while (events.front().time < start) events.pop_front();
}

void PolyView::change_range()
{
	range_entry.set_value(range);
}

void PolyView::range_ok()
{
	range = range_entry.get_value();
	rdlg->hide();
	queue_draw();
}

void PolyView::file_select(std::string name)
{
	fn_entry.set_text(name);
}

void PolyView::file_ok()
{
	std::string fn = fn_entry.get_text();
	if (! fn.size()) return;
	std::ofstream f(fn.c_str());
	size_t sz = vol->get_nnum();
	double start = vol->get_age() - range;
	f.precision(4);
	f << std::fixed;
	for (std::deque<Event>::iterator i = events.begin(); i != events.end(); i ++) {
		if (i->i >= sz || i->time <= start) continue;
		// f << i->time - start << '\t' << i->i << '\n';
		f << i->time << '\t' << i->i << '\n';
	}
	fdlg->hide();
}

PolyView::PolyView(MEAsim * vr) :
	vol(vr),
	range(5000),
	rdlg(0),
	fdlg(0)
{
	vr->fires.connect(sigc::mem_fun(this, & PolyView::add));
	menu.push_back("Change range");
	menu.push_back("Export data");
	// menu.push_back("Export SVG");
}

PolyView::~PolyView()
{
	delete rdlg;
	delete fdlg;
}
