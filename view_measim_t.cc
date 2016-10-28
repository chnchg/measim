// where the main function is...

#include "measim_run.hh"
#include "sim/gkviewt.hh"
#include "pot_view.hh"
#include "net_view.hh"
#include "view_misc.hh"
#include "poly_view.hh"
#include "initia_win.hh"

struct Update
{
	gltk::Window * w;
	double next;
	double interval;
	Update(gltk::Window * w, double interval) : w(w), next(0), interval(interval) {}
};

std::vector<Update> update_list;
Runnable * sys;

void update_reset()
{
	double cnt = sys->count();
	for (std::vector<Update>::iterator i = update_list.begin(); i != update_list.end(); i ++) {
		i->w->refresh();
		i->next = cnt + i->interval;
	}
}

void update_views()
{
	double cnt = sys->count();
	for (std::vector<Update>::iterator i = update_list.begin(); i != update_list.end(); i ++) {
		if (cnt < i->next) continue;
		i->w->refresh();
		i->next += i->interval;
	}
}

int main(int argc, char ** argv, char ** envp)
{
	gkViewT view(argc, argv, envp);
	MEAsimRun vol;
	sys = & vol;
	view.set_system(sys);
	gltk::Window * w;
	view.sys_dashed.connect(sigc::ptr_fun(update_views));
	view.set_title("View MEAsim");

	PotView pv(& vol);
	w = new gltk::Window;
	w->add(pv);
	pv.show();
	w->set_title("W-V space");
	w->show();
	update_list.push_back(Update(w, 1.0));

	NetView nv(& vol);
	w = new gltk::Window;
	w->add(nv);
	nv.show();
	w->set_title("network view");
	w->show();
	update_list.push_back(Update(w, 1.0));
	view.tk_changed("range").connect(sigc::mem_fun(mwlist, & std::vector<ModW>::clear));
	view.tk_changed("range").connect(sigc::mem_fun(* w, & gltk::Window::refresh));
	view.sys_changed.connect(sigc::mem_fun(mwlist, & std::vector<ModW>::clear));
	view.sys_changed.connect(sigc::ptr_fun(& update_reset));

	PolyView poly(& vol);
	w = new gltk::Window;
	w->add(poly);
	poly.show();
	w->set_default_size(800, 320);
	w->set_title("firing polygram");
	w->show();
	update_list.push_back(Update(w, 8.0));

	// initia window
	InitiaWin ini(vol.get_initia());
	view.add_action("Init", sigc::mem_fun(ini, & InitiaWin::show));
	view.sys_loaded.connect(sigc::mem_fun(ini, & InitiaWin::sync_mod));
	ini.select_ok.connect(sigc::mem_fun(sys, & Runnable::init));
	ini.select_ok.connect(sigc::mem_fun(poly, & PolyView::clear));
	ini.select_ok.connect(view.sys_changed.make_slot());

	return view.execute();
}

/*
 * Following is copied from
 *     http://stackoverflow.com/questions/9555375/qt-cmake-vc-produce-a-command-prompt
 */
#ifdef _WIN32
extern "C" {
#include <shellapi.h>
}
class Win32CommandLineConverter {
private:
	std::unique_ptr<char*[]> argv_;
	std::vector<std::unique_ptr<char[]>> storage_;
public:
	Win32CommandLineConverter()
	{
		LPWSTR cmd_line = GetCommandLineW();
		int argc;
		LPWSTR* w_argv = CommandLineToArgvW(cmd_line, &argc);
		argv_ = std::unique_ptr<char*[]>(new char*[argc]);
		storage_.reserve(argc);
		for(int i=0; i<argc; ++i) {
			storage_.push_back(ConvertWArg(w_argv[i]));
			argv_[i] = storage_.back().get();
		}
		LocalFree(w_argv);
	}
	int argc() const {return static_cast<int>(storage_.size());}
	char ** argv() const {return argv_.get();}
	static std::unique_ptr<char[]> ConvertWArg(LPWSTR w_arg)
	{
		int size = WideCharToMultiByte(CP_UTF8, 0, w_arg, -1, nullptr, 0, nullptr, nullptr);
		std::unique_ptr<char[]> ret(new char[size]);
		WideCharToMultiByte(CP_UTF8, 0, w_arg, -1, ret.get(), size, nullptr, nullptr);
		return ret;
	}
};
int CALLBACK WinMain(HINSTANCE /* hInstance */, HINSTANCE /* hPrevInstance */, LPSTR /* lpCmdLine */, int /* nCmdShow */)
{
	Win32CommandLineConverter cmd_line;
	return main(cmd_line.argc(), cmd_line.argv(), 0);
}
#endif
