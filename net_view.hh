#ifndef NET_VIEW_HH
#define NET_VIEW_HH
#include "measim_run.hh"
#include <gltk.hh>
#include <deque>
class NetView :
	public virtual gltk::Widget
{
	MEAsim const * vol;
	void on_draw() /* final */;
	struct TrailPoint
	{
		double x;
		double y;
		double r;
	};
	std::deque<TrailPoint> trail;
public:
	NetView(MEAsim * vr) : vol(vr) {}
};
#endif
