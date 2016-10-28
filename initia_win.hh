#ifndef INITIA_WIN_HH
#define INITIA_WIN_HH
#include "initia.hh"
#include <gltk.hh>
#include <prm/gui.hh>
class InitiaWin :
	public gltk::Window
{
	Initia & ini;
	gltk::ComboBoxText cbt; // module selector
	gltk::Box setup_box;
	prm::Gltk gk;
	gltk::Widget * active_setup;
	void change_mod();
	void update_setup();
public:
	InitiaWin(Initia & initia);
	sigc::signal<void> select_ok;
	void sync_mod(); // follow change in initia
};
#endif
