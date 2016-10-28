/*
 *  Dump trace of various variables
 */

#include "measim_run.hh"
#include <deque>
class MEAsimDump :
	virtual public MEAsimRun
{

	double resolution;
	void dash();
public:
	MEAsimDump();
};
