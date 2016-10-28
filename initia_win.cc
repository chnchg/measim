#include "initia_win.hh"

void InitiaWin::change_mod()
{
	std::string id = cbt.get_active_text();
	ini.set_module(cbt.get_active_text());
	update_setup();
}

void InitiaWin::update_setup()
{
	if (active_setup) {
		setup_box.remove(active_setup);
		active_setup = 0;
		gk.clear();
	}
	gk.add(ini.get_param());
	active_setup = & gk.widget();
	setup_box.pack_start(active_setup);
	active_setup->show();
	resize();
}

InitiaWin::InitiaWin(Initia & initia) :
	ini(initia),
	active_setup(0)
{
	std::vector<std::string> nl = initia.get_names();
	for (std::vector<std::string>::iterator i = nl.begin(); i != nl.end(); i ++) cbt.append(* i);
	sync_mod();
	using namespace gltk;
	Box * bx = new Box(ORIENTATION_VERTICAL);
	add(bx);
	bx->show();
	bx->pack_start(cbt);
	cbt.show();
	cbt.signal_changed().connect(sigc::mem_fun(this, & InitiaWin::change_mod));
	bx->pack_start(setup_box);
	setup_box.show();
	Box * hb = new Box(ORIENTATION_HORIZONTAL);
	bx->pack_start(hb);
	hb->show();
	Button * b = new Button("Cancel");
	b->signal_clicked().connect(sigc::mem_fun(this, & Window::hide));
	hb->pack_start(b);
	b->show();
	b = new Button("OK");
	b->signal_clicked().connect(sigc::mem_fun(this, & Window::hide));
	b->signal_clicked().connect(select_ok.make_slot());
	hb->pack_start(b);
	b->show();
}

void InitiaWin::sync_mod()
{
	cbt.set_active_text(ini.get_module());
	update_setup();
}
