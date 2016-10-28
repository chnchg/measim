#ifndef POLY_VIEW_HH
#define POLY_VIEW_HH
#include "measim_run.hh"
#include <s_entry.hh>
#include <gltk.hh>
#include <deque>

class PolyView :
	public virtual gltk::Widget,
	public virtual gltk::Popup
{
	MEAsim const * vol;
	double range;
	void on_draw();
	void on_menu_select(int item);
	struct Event
	{
		size_t i;
		double time;
		Event(size_t i, double time) : i(i), time(time) {}
	};
	std::deque<Event> events;
	void add(size_t i);
	gltk::Window * rdlg; // range dialog
	SEntry<double> range_entry;
	void change_range();
	void range_ok();
	gltk::Window * fdlg; // file dialog
	gltk::Entry fn_entry;
	void file_select(std::string name);
	void file_ok();
public:
	PolyView(MEAsim * vol);
	~PolyView();
	void clear() {events.clear();}
};

#endif
