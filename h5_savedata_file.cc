#include "h5_savedata_file.hh"
#include <vector>
#include <iostream>
#include <cerrno>
#include <cstring>
extern "C" {
#include <fcntl.h>
#if defined (_WIN32)
#include <io.h>
#define F_OK (0)
#define R_OK (2)
#else
#include <unistd.h>
#endif
}

#define mpotRank 2
#define fstampRank 1
#define LENGTH 100000

H5SimSaveFile::~H5SimSaveFile()
{
	done_save();	
}

void H5SimSaveFile::initH5simSaveFile()
{
	done_save();
	curSteps = 0;
	curSpkNum = 0;	
	open_save();
}
void H5SimSaveFile::close()
{
	if (is_open){
		is_open = false;
		H5Sclose(mpotS);
		H5Sclose(fstampS);
		H5Dclose(mpotD);
		H5Dclose(fstampD);
		H5Gclose(group);
		H5Fclose(file);
	}
}

void H5SimSaveFile::open_save(std::string path)
{
	close();
	std::string tmpfile = name + ".new~";
	file = H5Fcreate(tmpfile.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file < 0) {
		std::cerr << "fail to open hdf5 file: " << tmpfile << '\n';
		throw;
	}
	group = H5Gopen(file, path.c_str(), H5P_DEFAULT);

	hsize_t dimv[2] = {LENGTH,nnum};
	hsize_t dimf[] = {LENGTH};
	hsize_t mdim1[2] = {1,nnum};
	hsize_t mdim2[] = {nnum};
	hsize_t maxdimv[2] = {H5S_UNLIMITED,nnum};
	hsize_t maxdimf[] = {H5S_UNLIMITED};
	hsize_t chunk_dimv[2] ={1, nnum};
	hsize_t chunk_dimf[] ={nnum};
	herr_t status;

	mpotS = H5Screate_simple(mpotRank, dimv, maxdimv);
	hid_t mpotP = H5Pcreate (H5P_DATASET_CREATE);
	status = H5Pset_layout(mpotP, H5D_CHUNKED);
	status = H5Pset_chunk( mpotP, mpotRank, chunk_dimv);
	mpotD = H5Dcreate2(file, "/mpot", H5T_NATIVE_DOUBLE, mpotS, H5P_DEFAULT, mpotP, H5P_DEFAULT);
	H5Pclose(mpotP);
	
	fstampS = H5Screate_simple(fstampRank, dimf, maxdimf);
	hid_t fstampT = H5Tcreate (H5T_COMPOUND, sizeof(fstamp));
	H5Tinsert(fstampT, "time", HOFFSET(fstamp, time), H5T_NATIVE_DOUBLE);
	H5Tinsert(fstampT, "neuron", HOFFSET(fstamp, nrn), H5T_NATIVE_ULONG);
	hid_t fstampP = H5Pcreate (H5P_DATASET_CREATE);
	status = H5Pset_chunk( fstampP, fstampRank, chunk_dimf);	
	fstampD = H5Dcreate2(file, "/fstamp", fstampT, fstampS, H5P_DEFAULT, fstampP, H5P_DEFAULT);
	H5Pclose(fstampP);
	H5Tclose(fstampT);

	is_open = true;
}

void H5SimSaveFile::done_save()
{
	if (is_open){
		hsize_t size1[2] = {curSteps,nnum};
		hsize_t size2[] = {curSpkNum};
		herr_t status;
		status = H5Dset_extent(mpotD, size1);
		status = H5Dset_extent(fstampD, size2);
		
		close();
		if (access(name.c_str(), F_OK) == 0) { // old data file exists?
			std::string bak = name + "~";
			int ret = rename(name.c_str(), bak.c_str()); // backup old file
			if (ret && access(bak.c_str(), F_OK) == 0) { // error renaming file, and old backup exists
				unlink(bak.c_str()); // try removing old backup first...
				rename(name.c_str(), bak.c_str()); // do again...
			}
		}
		rename((name + ".new~").c_str(), name.c_str());
	}
}

void H5SimSaveFile::dump_mpot(double * pot)
{
	hid_t dataspce, filespace;
	herr_t status;
	hsize_t count2D[2];              /* size of the hyperslab in the file */
	hsize_t offset2D[2];             /* hyperslab offset in the file */
	
	// resize the dataset

//	if (curSteps%LENGTH==0){
	if ((curSteps+1)*nnum* H5Tget_size(H5Dget_type(mpotD)) > H5Dget_storage_size(mpotD)){
		hsize_t size1[2] = {curSteps+LENGTH,nnum};
		status = H5Dset_extent(mpotD, size1);
	}

	filespace = H5Dget_space(mpotD);
	// Select a hyperslab for file sapce
	offset2D[0] = curSteps;
	offset2D[1] = 0;
	count2D[0] = 1;
	count2D[1] = nnum;
	
	status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset2D, NULL,count2D, NULL);
	if (status<0)
		std::cout << "Can't select region!(mpot)\n";
	hsize_t dims2[2]={1,nnum};
	dataspce = H5Screate_simple(mpotRank, dims2, NULL);
	
	// save membrane potential
	status = H5Dwrite(mpotD, H5T_NATIVE_DOUBLE, dataspce, filespace,H5P_DEFAULT, pot);	
	curSteps += 1;
	H5Sclose(dataspce);
}

void H5SimSaveFile::dump_fire(bool const * fire, double curAge)
{
	hid_t dataspce, filespace;
	herr_t status;
	hsize_t count1D[1];              /* size of the hyperslab in the file */
	hsize_t offset1D[1];             /* hyperslab offset in the file */
	

	
	//resize the dataset
	if ((curSpkNum+nnum)*2* H5Tget_size(H5Dget_type(fstampD)) > H5Dget_storage_size(fstampD)){
		hsize_t size2[] = {curSpkNum+LENGTH};
		status = H5Dset_extent(fstampD, size2);
	}
	
	// extract spikes
	std::vector<fstamp> spks;
	fstamp tempfs;
	for(size_t i = 0; i<nnum;i++){
		if (fire[i]){
			spks.push_back(fstamp{curAge,i});
		}
	}
	fstamp * fstampptr = new fstamp[spks.size()];
	for(size_t i = 0; i<spks.size();i++){
		fstampptr[i].time =spks.at(i).time;
		fstampptr[i].nrn =spks.at(i).nrn;
	}

	filespace = H5Dget_space(fstampD);
	// Select a hyperslab
	offset1D[0] = curSpkNum;
	count1D[0] = spks.size();
	status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset1D, NULL,count1D, NULL);
	if (status<0)
		std::cout << "Can't select region!(fstamp)\n";
	hsize_t dims1[]={spks.size()};
	dataspce = H5Screate_simple(fstampRank, dims1, NULL);
	
	// save spike train
	status = H5Dwrite(fstampD, H5Dget_type(fstampD), dataspce, filespace,H5P_DEFAULT, fstampptr);
	curSpkNum += spks.size();
	H5Sclose(dataspce);
	delete[] fstampptr;
}
