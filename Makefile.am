# measim/Makefile.am
AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4
CLEANFILES = *~ */*~
SUBDIRS = cst
LDADD = libvol.a cst/libcst.a $(HDF5_LDFLAGS)
AM_CXXFLAGS = -I$(srcdir)/cst $(PROG_CFLAGS)
include cst/make_main.am

noinst_LIBRARIES = libvol.a
libvol_a_SOURCES = \
	measim.cc measim.hh measim_run.cc measim_run.hh h5_savedata_file.cc h5_savedata_file.hh\
	initia.cc initia.hh mod_init_modules.cc mod_init_modules.hh

bin_PROGRAMS =
bin_PROGRAMS += rvol
nodist_rvol_SOURCES = measim_run-runner.cc

bin_PROGRAMS += vdmp
vdmp_SOURCES = measim_dump.hh measim_dump.cc
nodist_vdmp_SOURCES = measim_dump-runner.cc

if HAVE_GLUT
noinst_LIBRARIES += libgui.a
libgui_a_SOURCES = \
	pot_view.cc pot_view.hh net_view.cc net_view.hh \
	view_misc.cc view_misc.hh poly_view.cc poly_view.hh \
	initia_win.cc initia_win.hh

bin_PROGRAMS += vvol
vvol_SOURCES = view_measim.cc
vvol_LDADD = libgui.a $(LDADD) $(GLUT_LIBS) $(GL_LIBS)

# threaded version
bin_PROGRAMS += vvlt
vvlt_SOURCES = view_measim_t.cc
vvlt_LDADD = libgui.a $(LDADD) $(GLUT_LIBS) $(GL_LIBS)
endif
