#ifndef POT_VIEW_HH
#define POT_VIEW_HH
#include "measim_run.hh"
#include <gltk.hh>
class PotView :
	public virtual gltk::Widget
{
	MEAsimRun * vol;
	void on_draw() final;
public:
	PotView(MEAsimRun * vol) : vol(vol) {}
};
#endif
