#ifndef H5SAVEDATAFILE_HH
#define H5SAVEDATAFILE_HH 1
#include <hdf5.h>
#include <string>



class H5SimSaveFile
{
	typedef struct {
		double time;
		size_t    nrn;
//		fstamp(double n1, size_t n2):time(n1),nrn(n2){}
	}fstamp;
	std::string name;
	bool is_open;
	hid_t file, group;
	hid_t mpotD,mpotS, fstampD, fstampS;
	size_t curSteps,curSpkNum,nnum;
public:
	H5SimSaveFile(std::string name = "") : name(name),nnum(0),is_open(false),curSteps(0),curSpkNum(0) {}
	H5SimSaveFile(std::string name, size_t n ) : name(name),nnum(n),is_open(false),curSteps(0),curSpkNum(0) {}
	~H5SimSaveFile();
	void open_save(std::string path = "/");
	void done_save();
	void close();
	void set_name(std::string n) {name=n;}
	void set_nnum(size_t n){nnum = n;}
	void initH5simSaveFile();
	void dump_fire(bool const * fire, double curAge);
	void dump_mpot(double * pot);
	size_t get_curSteps(){return curSteps;}
	size_t get_curSpkNum(){return curSpkNum;}
	std::string get_name(){return name;}
	bool get_isOpen(){return is_open;}
};
#endif