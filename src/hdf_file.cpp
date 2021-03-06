#include <sstream>
#include <cassert>
#include <mutex>
#include <exception>
#include <fstream>
#include <sstream>
#include "boost/uuid/uuid.hpp"
#include "boost/uuid/uuid_generators.hpp"
#include "hdf_file.hpp"
#include "hdf5.h"
#include "smv.hpp"


herr_t IterateTrajectories(hid_t group_id, const char* group_name,
  const H5L_info_t* info, void* op_data) {
  int* chosen_idx=static_cast<int*>(op_data);

  std::string sname{group_name};
  if (sname.substr(0,10)==std::string("trajectory")) {
    if (sname.size()>10) {
      std::istringstream convert{sname.substr(10)};
      int val{0};
      try {
        convert>>val;
        *chosen_idx=std::max(*chosen_idx, val+1);
      } catch (std::exception& e) {
        BOOST_LOG_TRIVIAL(warning)<<"Could not convert "<<sname;
      }
    } else {
      *chosen_idx=std::max(*chosen_idx, 1);
    }
  }
  return 0;
}


using EndType=std::tuple<std::vector<double>,std::vector<int64_t>>;
// Take just the last value from a dataset.
herr_t WalkFindEnd(hid_t group_id, const char* group_name,
  const H5L_info_t* info, void* op_data) {
  auto end=static_cast<EndType*>(op_data);

  hid_t dsg_id=H5Gopen(group_id, group_name, H5P_DEFAULT);
  hid_t ds_id=H5Dopen(dsg_id, "seirtotaltimes", H5P_DEFAULT);
  if (ds_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not open seirtotal for "<<group_name;
    return -1;
  }
  hid_t space_id=H5Dget_space(ds_id);
  int ndims=H5Sget_simple_extent_ndims(space_id);
  hsize_t dims[ndims];
  hsize_t maxdims[ndims];
  ndims=H5Sget_simple_extent_dims(space_id, dims, maxdims);
  if (ndims<0) {
    BOOST_LOG_TRIVIAL(error)<<"Couldn't get extent of dataset?";
    return -2;
  }

  hsize_t one{1};
  hid_t mem_space=H5Screate_simple(1, &one, NULL);
  hid_t sub_file=H5Scopy(space_id);
  hsize_t offset(dims[0]-1);
  herr_t offset_status=H5Sselect_elements(sub_file, H5S_SELECT_SET, 1, &offset);
  if (offset_status<0) {
    BOOST_LOG_TRIVIAL(error)<<"Couldn't offset dataspace";
    return -3;
  }

  double value;
  H5Dread(ds_id, H5T_NATIVE_DOUBLE, mem_space, sub_file, H5P_DEFAULT, &value);
  std::get<0>(*end).push_back(value);
  std::get<1>(*end).push_back(dims[0]);

  H5Sclose(sub_file);
  H5Sclose(mem_space);
  H5Sclose(space_id);
  H5Dclose(ds_id);
  H5Gclose(dsg_id);
  return 0;
}


herr_t TrajectoryNames(hid_t group_id, const char* relative_name,
  const H5L_info_t* info, void* op_data) {
  typedef std::vector<std::string> Names;
  Names* name=static_cast<Names*>(op_data);
  name->push_back(std::string{relative_name});
  return 0;
}


HDFFile::HDFFile(const std::string& filename)
: filename_(filename), open_(false), file_id_(0), trajectory_group_(0),
  comp_cnt_{5}, image_group_(0)
{}

HDFFile::~HDFFile() {}

bool HDFFile::Open(bool truncate) {
  if (truncate) {
    file_id_=H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
        H5P_DEFAULT);
  } else {
    std::ifstream file_exists(filename_.c_str());
    bool exists=file_exists.good();
    file_exists.close();

    if (exists) {
      file_id_=H5Fopen(filename_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
      BOOST_LOG_TRIVIAL(debug)<<"File exists: "<<filename_;
    } else {
      BOOST_LOG_TRIVIAL(debug)<<"Creating file: "<<filename_;
      file_id_=H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
          H5P_DEFAULT);
    }
  }
  if (file_id_<0) return false;

  open_=true;

  int chosen_idx=0;

  if (!truncate) {
    // Is there already a trajectory in this file?
    herr_t find_status=H5Literate_by_name(file_id_, "/", H5_INDEX_NAME,
      H5_ITER_INC, NULL, IterateTrajectories, &chosen_idx, H5P_DEFAULT);
    BOOST_LOG_TRIVIAL(debug)<<"HDFFile::Open chosen index "<<chosen_idx;
    if (find_status<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not iterate over groups in file.";
    }
  }
  std::stringstream trajname;
  trajname<<"/trajectory";
  if (chosen_idx>0) {
    trajname << chosen_idx;
  }
  BOOST_LOG_TRIVIAL(info)<<"Writing to file directory "<<trajname.str();
  trajectory_group_=H5Gcreate(file_id_, trajname.str().c_str(), H5P_DEFAULT,
    H5P_DEFAULT, H5P_DEFAULT);

  WriteUUIDTo(trajectory_group_);

  return true;
}

bool HDFFile::OpenRead(bool readwrite) {
  std::ifstream file_exists(filename_.c_str());
  bool exists=file_exists.good();
  file_exists.close();
  if (!exists) {
    BOOST_LOG_TRIVIAL(error)<<"File doesn't exist to open "<<filename_;
    return false;
  }

  if (readwrite) {
    BOOST_LOG_TRIVIAL(debug)<<"Open readwrite: "<<filename_;
    file_id_=H5Fopen(filename_.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  } else {
    BOOST_LOG_TRIVIAL(debug)<<"Open read-only: "<<filename_;
    file_id_=H5Fopen(filename_.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  }

  if (file_id_<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not open file for reading "<<filename_;
    return false;
  }
  if (H5Lexists(file_id_, "/trajectory", H5P_DEFAULT)!=0) {
    BOOST_LOG_TRIVIAL(debug)<<"Opening /trajectory";
    trajectory_group_=H5Gopen(file_id_, "/trajectory", H5P_DEFAULT);
  }

  open_=true;
  return open_;
}

bool HDFFile::WriteExecutableData(const std::map<std::string,std::string>& compile,
  const boost::program_options::basic_parsed_options<char>& options,
  const std::vector<int64_t>& siri) const {
  std::unique_lock<std::mutex> only_me(single_writer_);

  {
    hsize_t adims=1;
    hid_t dspace_id=H5Screate_simple(1, &adims, NULL);

    for (const auto& kv : compile) {
      hid_t strtype=H5Tcopy(H5T_C_S1);
      herr_t strstatus=H5Tset_size(strtype, kv.second.size());
      if (strstatus<0) {
        BOOST_LOG_TRIVIAL(error)
          <<"Could not create string for executable data.";
          return false;
      }

      hid_t attr0_id=H5Acreate2(trajectory_group_, kv.first.c_str() , strtype,
        dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      herr_t atstatus=H5Awrite(attr0_id, strtype, kv.second.c_str());
      if (atstatus<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not write attribute "<<kv.first;
        return false;
      }
      H5Tclose(strtype);
      H5Aclose(attr0_id);
    }
    H5Sclose(dspace_id);
  }

  {
    hsize_t sdims=siri.size();
    hid_t sirspace_id=H5Screate_simple(1, &sdims, NULL);

    hid_t attr1_id=H5Acreate2(trajectory_group_, "Initial Values",
      H5T_STD_I64LE, sirspace_id, H5P_DEFAULT, H5P_DEFAULT);
    herr_t at1status=H5Awrite(attr1_id, H5T_NATIVE_LONG, &siri[0]);
    if (at1status<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write attribute Initial Values";
      return false;
    }
    H5Sclose(sirspace_id);
    H5Aclose(attr1_id);
  }

  // We save the command-line options used to call the program.
  std::stringstream optstring;
  optstring << "<options>";
  for (auto& opt : options.options) {
    optstring << "<option><name>" << opt.string_key << "</name>";
    optstring << "<values>";
    for (auto& v : opt.value) {
      optstring << "<value>" << v << "</value>";
    }
    optstring <<"</values></option>";
  }
  optstring << "</options>";
  
  {
    hsize_t odims=1;
    hid_t ospace_id=H5Screate_simple(1, &odims, NULL);

    hid_t ostrtype=H5Tcopy(H5T_C_S1);
    herr_t ostrstatus=H5Tset_size(ostrtype, optstring.str().size());
    if (ostrstatus<0) {
      BOOST_LOG_TRIVIAL(error)
        <<"Could not create string for executable data.";
        return false;
    }

    hid_t oattr_id=H5Acreate2(trajectory_group_, "Options", ostrtype,
      ospace_id, H5P_DEFAULT, H5P_DEFAULT);
    herr_t ostatus=H5Awrite(oattr_id, ostrtype, optstring.str().c_str());
    if (ostatus<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write attribute Options";
      return false;
    }
    H5Sclose(ospace_id);
    H5Tclose(ostrtype);
    H5Aclose(oattr_id);
  }

  return true;
}


bool HDFFile::WriteEnsembleVariables(int64_t elapsed_ns) const {
  hsize_t odims=1;
  hid_t ospace_id=H5Screate_simple(1, &odims, NULL);
  hid_t attr0_id=H5Acreate2(trajectory_group_, "totalelapsedwall",
    H5T_STD_I64LE, ospace_id, H5P_DEFAULT, H5P_DEFAULT);
  if (attr0_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not write walltime";
    H5Sclose(ospace_id);    
    return false;
  }
  herr_t wstatus=H5Awrite(attr0_id, H5T_STD_I64LE, &elapsed_ns);
  if (wstatus<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not write walltime";
    return false;
  }

  H5Aclose(attr0_id);
  H5Sclose(ospace_id);
  return wstatus>=0;
}


bool HDFFile::Close() {
  if (open_) {
    if (trajectory_group_>0) {
      herr_t group_status=H5Gclose(trajectory_group_);
      if (group_status<0) {
        BOOST_LOG_TRIVIAL(warning)<<"Could not close HDF5 group "<<group_status;
        return false;
      }
    }
    if (image_group_>0) {
      herr_t group_status=H5Gclose(image_group_);
      if (group_status<0) {
        BOOST_LOG_TRIVIAL(warning)<<"Could not close HDF5 image group "
            <<group_status;
        return false;
      }
    }
    herr_t status=H5Fclose(file_id_);
    if (status<0) {
      BOOST_LOG_TRIVIAL(warning)<<"Could not close HDF5 file "<<status;
      return false;
    }
    open_=false;
  }
  return true;
}


std::vector<std::string> HDFFile::Trajectories() const {
  std::vector<std::string> name;
  herr_t find_status=H5Literate(trajectory_group_, H5_INDEX_NAME,
    H5_ITER_NATIVE, NULL, TrajectoryNames, &name);
  if (find_status<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not iterate over groups in file.";
  }
  return name;
}

/*! How many events and what end times for each trajectory in a file.
 */
std::tuple<std::vector<double>,std::vector<int64_t>> HDFFile::EndTimes() const {
  H5G_info_t group_info;
  herr_t gi_status=H5Gget_info(trajectory_group_, &group_info);
  assert(gi_status>=0);

  EndType data;
  std::get<0>(data).reserve(group_info.nlinks);
  std::get<1>(data).reserve(group_info.nlinks);

  herr_t iter_status=H5Literate(trajectory_group_, H5_INDEX_NAME,
    H5_ITER_NATIVE, NULL, WalkFindEnd, &data);
  if (iter_status<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not iterate over trajectories to find "
        <<"end times "<<iter_status;
  }

  return data;
}


std::vector<int64_t> HDFFile::InitialValues() const {
  hid_t attr1_id=H5Aopen_by_name(file_id_, "/trajectory",
      "Initial Values", H5P_DEFAULT, H5P_DEFAULT);
  auto siri=std::vector<int64_t>(comp_cnt_, 0);

  herr_t at1status=H5Aread(attr1_id, H5T_NATIVE_LONG, &siri[0]);
  if (at1status<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not read attribute Initial Values";
    return siri;
  }
  H5Aclose(attr1_id);
  return siri;
}

std::vector<int64_t> HDFFile::LoadInitialPen(
    const std::string dataset_name) const {
  std::stringstream initial_name_str;
  initial_name_str<<"/trajectory/"<<dataset_name<<"/initialpencount";
  const char* initial_name=initial_name_str.str().c_str();

  hid_t ds_id=H5Dopen(file_id_, initial_name, H5P_DEFAULT);
  if (ds_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not read dataset from file: "
        <<initial_name_str.str();
    return std::vector<int64_t>{0};
  }
  hid_t space_id=H5Dget_space(ds_id);
  int ndims=H5Sget_simple_extent_ndims(space_id);
  if (ndims!=2) {
    BOOST_LOG_TRIVIAL(error)<<"initial values dataset has wrong number of "
        <<"dimensions? "<<ndims;
  }
  hsize_t dims[ndims];
  hsize_t maxdims[ndims];
  ndims=H5Sget_simple_extent_dims(space_id, dims, maxdims);
  if (ndims<0) {
    BOOST_LOG_TRIVIAL(error)<<"Couldn't get extent of dataset?";
  }
  std::vector<int64_t> initial(dims[0]*dims[1], 0);
  herr_t read_status=H5Dread(ds_id, H5T_STD_I64LE, space_id, H5S_ALL,
      H5P_DEFAULT, &initial[0]);

  H5Sclose(space_id);
  H5Dclose(ds_id);
  return initial;
}

std::array<int64_t, 2> HDFFile::EventsInFile() const {
  int64_t trajectory_cnt{0};
  int64_t event_cnt{0};
  auto trajectories=Trajectories();
  for (auto traj_name : trajectories) {
    std::stringstream initial_name_str;
    initial_name_str<<"/trajectory/"<<traj_name<<"/trajectory";
    const char* initial_name=initial_name_str.str().c_str();

    hid_t ds_id=H5Dopen(file_id_, initial_name, H5P_DEFAULT);
    if (ds_id<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not read dataset from file: "
          <<initial_name_str.str()<<" error "<<ds_id;
      return {0, 0};
    }

    hid_t space_id=H5Dget_space(ds_id);
    int ndims=H5Sget_simple_extent_ndims(space_id);
    if (ndims!=1) {
      BOOST_LOG_TRIVIAL(error)<<"initial values dataset has wrong number of "
          <<"dimensions? "<<ndims;
      return {0, 0};
    }
    hsize_t dims[ndims];
    hsize_t maxdims[ndims];
    ndims=H5Sget_simple_extent_dims(space_id, dims, maxdims);
    if (ndims<0) {
      BOOST_LOG_TRIVIAL(error)<<"Couldn't get extent of dataset?";
      return {0, 0};
    }
    trajectory_cnt+=1;
    event_cnt+=dims[0];

    H5Dclose(space_id);
    H5Dclose(ds_id);
  }
  return {trajectory_cnt, event_cnt};
}


bool HDFFile::LoadTrajectoryCounts(const std::string& name,
    std::vector<int64_t>& seirc, int64_t& cnt) const {
  std::stringstream name_str;
  name_str<<name<<"/seirtotal";
  BOOST_LOG_TRIVIAL(debug)<<"Reading "<<name_str.str();
  hid_t ds_id=H5Dopen(trajectory_group_, name_str.str().c_str(), H5P_DEFAULT);
  if (ds_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not read dataset "<<name
      <<" error "<<ds_id;
    return false;
  }
  hid_t space_id=H5Dget_space(ds_id);
  int ndims=H5Sget_simple_extent_ndims(space_id);
  if (ndims!=2) {
    BOOST_LOG_TRIVIAL(error)<<"initial values dataset has wrong number of "
        <<"dimensions? "<<ndims;
    return false;
  }
  hsize_t dims[ndims];
  hsize_t maxdims[ndims];
  ndims=H5Sget_simple_extent_dims(space_id, dims, maxdims);
  if (ndims<0) {
    BOOST_LOG_TRIVIAL(error)<<"Couldn't get extent of dataset?";
    return false;
  }
  if (seirc.size()<dims[0]*dims[1]) {
    BOOST_LOG_TRIVIAL(debug)<<"LoadTrajectoryCounts: Resizing array";
    seirc.resize(dims[0]*dims[1]);
  }
  cnt=dims[0];
  herr_t read_status=H5Dread(ds_id, H5T_STD_I64LE, space_id, H5S_ALL,
      H5P_DEFAULT, seirc.data());
  H5Sclose(space_id);
  H5Dclose(ds_id);
  return true;
}

bool HDFFile::LoadTrajectoryTimes(const std::string& name,
    std::vector<double>& time, int64_t& cnt) const {
  std::stringstream name_str;
  name_str<<name<<"/seirtotaltimes";
  BOOST_LOG_TRIVIAL(debug)<<"Reading "<<name_str.str();
  hid_t ds_id=H5Dopen(trajectory_group_, name_str.str().c_str(), H5P_DEFAULT);
  if (ds_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not read dataset "<<name;
    return false;
  }
  hid_t space_id=H5Dget_space(ds_id);
  int ndims=H5Sget_simple_extent_ndims(space_id);
  if (ndims!=1) {
    BOOST_LOG_TRIVIAL(error)<<"initial values dataset has wrong number of "
        <<"dimensions? "<<ndims;
    return false;
  }
  hsize_t dims[ndims];
  hsize_t maxdims[ndims];
  ndims=H5Sget_simple_extent_dims(space_id, dims, maxdims);
  if (ndims<0) {
    BOOST_LOG_TRIVIAL(error)<<"Couldn't get extent of dataset?";
    return false;
  }
  if (time.size()<dims[0]) {
    BOOST_LOG_TRIVIAL(debug)<<"LoadTrajectoryTimes: resizing array";
    time.resize(dims[0]);
  }
  cnt=dims[0];
  herr_t read_status=H5Dread(ds_id, H5T_NATIVE_DOUBLE, space_id, H5S_ALL,
      H5P_DEFAULT, time.data());
  H5Sclose(space_id);
  H5Dclose(ds_id);
  return true;
}


std::vector<TrajectoryEntry> HDFFile::LoadTrajectoryFromPens(
    const std::string dataset_name) const {
  std::stringstream initial_name_str;
  initial_name_str<<"/trajectory/"<<dataset_name<<"/trajectory";
  const char* initial_name=initial_name_str.str().c_str();

  hid_t ds_id=H5Dopen(file_id_, initial_name, H5P_DEFAULT);
  if (ds_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not read dataset from file: "
        <<initial_name_str.str();
    return std::vector<TrajectoryEntry>(0);
  }
  hid_t space_id=H5Dget_space(ds_id);
  int ndims=H5Sget_simple_extent_ndims(space_id);
  if (ndims!=1) {
    BOOST_LOG_TRIVIAL(error)<<"initial values dataset has wrong number of "
        <<"dimensions? "<<ndims;
  }
  hsize_t dims[ndims];
  hsize_t maxdims[ndims];
  ndims=H5Sget_simple_extent_dims(space_id, dims, maxdims);
  if (ndims<0) {
    BOOST_LOG_TRIVIAL(error)<<"Couldn't get extent of dataset?";
  }
  std::vector<PenTrajectory> events{dims[0]};
  hid_t type_id=H5Dget_type(ds_id);
  herr_t read_status=H5Dread(ds_id, type_id, space_id, H5S_ALL,
      H5P_DEFAULT, &events[0]);

  std::vector<TrajectoryEntry> trajectory{dims[0]+1};
  auto initial_pen=LoadInitialPen(dataset_name);
  std::vector<int64_t> initial{comp_cnt_};
  initial.assign(comp_cnt_, 0);
  for (size_t ip=0; ip<initial_pen.size(); ++ip) {
    for (size_t cidx=0; cidx<comp_cnt_; ++cidx) {
      initial[cidx]+=initial_pen[ip+cidx];
    }
    ip+=comp_cnt_;
  }

  trajectory[0].s=initial[0];
  trajectory[0].e=initial[1];
  trajectory[0].i=initial[2];
  trajectory[0].r=initial[3];
  trajectory[0].t=0.;
  for (size_t i=0; i<dims[0]; ++i) {
    switch (events[i].transition) {
      case 0:
        initial[0]-=1;
        initial[1]+=1;
        break;
      case 1:
        initial[1]-=1;
        initial[2]+=1;
        break;
      case 2:
        initial[2]-=1;
        initial[3]+=1;
        initial[4]-=1;
        break;
      case 3:
        initial[4]+=1;
        break;
      default:
        BOOST_LOG_TRIVIAL(error)<<"Unknown transition type."
            <<events[i].transition;
        break;
    }
    trajectory[i+1].s=initial[0];
    trajectory[i+1].e=initial[1];
    trajectory[i+1].i=initial[2];
    trajectory[i+1].r=initial[3];
    trajectory[i+1].c=initial[4];
    trajectory[i+1].t=events[i].time;
  }

  H5Dclose(type_id);
  H5Sclose(space_id);
  H5Dclose(ds_id);
  return trajectory;
}



bool HDFFile::Save2DPDF(const std::vector<double>& interpolant,
    const std::vector<double>& x, const std::vector<double>& y,
    std::string name) {
  ImageGroup();

  std::stringstream xstr;
  xstr<<name<<"x";
  Write1DFloat(image_group_, x, xstr.str());
  std::stringstream ystr;
  ystr<<name<<"y";
  Write1DFloat(image_group_, y, ystr.str());

  hsize_t dims[2];
  dims[0]=y.size();
  dims[1]=x.size();
  assert(interpolant.size()==dims[0]*dims[1]);
  BOOST_LOG_TRIVIAL(debug)<<"Saving ensemble with dims "
      <<dims[0]<<" "<<dims[1];
  hid_t dspace=H5Screate_simple(2, dims, NULL);
  hid_t ds_id=H5Dcreate(image_group_, name.c_str(), H5T_IEEE_F64LE,
      dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  herr_t interpolant_write=H5Dwrite(ds_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &interpolant[0]);
  H5Dclose(ds_id);
  H5Sclose(dspace);
  return true;
}

bool HDFFile::Save1DArray(const std::vector<double>& x,
    const std::string& name) {
  ImageGroup();
  return this->Write1DFloat(image_group_, x, name);
}

bool HDFFile::Save1DArray(const std::vector<int64_t>& x,
    const std::string& name) {
  ImageGroup();
  return this->Write1DInt(image_group_, x, name);
}

bool HDFFile::Save2DArray(const std::vector<int64_t>& x, int cols, 
      const std::string& name) {
  ImageGroup();
  return this->Write2DInt(image_group_, x, cols, name);
}

bool HDFFile::WriteUUIDTo(hid_t group) const {
  auto gen=boost::uuids::random_generator();
  boost::uuids::uuid tag(gen());
  std::stringstream uuidstring;
  uuidstring << tag;

  hsize_t odims=1;
  hid_t ospace_id=H5Screate_simple(1, &odims, NULL);

  hid_t strtype=H5Tcopy(H5T_C_S1);
  BOOST_LOG_TRIVIAL(debug)<<"uuidstring size " <<uuidstring.str().size();
  herr_t strstatus=H5Tset_size(strtype, uuidstring.str().size());
  if (strstatus<0) {
    BOOST_LOG_TRIVIAL(error)
      <<"Could not create string for executable data.";
      return false;
  }

  hid_t attr0_id=H5Acreate2(group, "uuid", strtype,
    ospace_id, H5P_DEFAULT, H5P_DEFAULT);
  herr_t atstatus=H5Awrite(attr0_id, strtype, uuidstring.str().c_str());
  if (atstatus<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not write attribute "<<"uuid";
    return false;
  }
  H5Tclose(strtype);
  H5Aclose(attr0_id);
  H5Sclose(ospace_id);
  return true;
}

bool HDFFile::WriteImageAttribute(std::string attrib, std::string name) {
  ImageGroup();
  hsize_t adims=1;
  hid_t space_id=H5Screate_simple(1, &adims, NULL);
  hid_t strtype=H5Tcopy(H5T_C_S1);
  herr_t strstatus=H5Tset_size(strtype, attrib.size());
  if (strstatus<0) {
    BOOST_LOG_TRIVIAL(error)
      <<"Could not create string for image attribute.";
      return false;
  }
  hid_t attr0_id=H5Acreate2(image_group_, name.c_str() , strtype,
    space_id, H5P_DEFAULT, H5P_DEFAULT);
  herr_t atstatus=H5Awrite(attr0_id, strtype, attrib.c_str());
  if (atstatus<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not write attribute "<<name;
    return false;
  }
  H5Tclose(strtype);
  H5Sclose(space_id);
  return true;
}

bool HDFFile::Write1DFloat(hid_t group, const std::vector<double>& x,
    const std::string& name) {
  hsize_t dims[1];
  dims[0]=x.size();
  hid_t x_dspace=H5Screate_simple(1, dims, NULL);
  hid_t x_id=H5Dcreate(group, name.c_str(), H5T_IEEE_F64LE,
      x_dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  herr_t x_write=H5Dwrite(x_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &x[0]);
  H5Dclose(x_id);
  H5Sclose(x_dspace);
  return true;
}


bool HDFFile::Write1DInt(hid_t group, const std::vector<int64_t>& x,
    const std::string& name) {
  hsize_t dims[1];
  dims[0]=x.size();
  hid_t x_dspace=H5Screate_simple(1, dims, NULL);
  hid_t x_id=H5Dcreate(group, name.c_str(), H5T_STD_I64LE,
      x_dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  herr_t x_write=H5Dwrite(x_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &x[0]);
  H5Dclose(x_id);
  H5Sclose(x_dspace);
  return true;
}


bool HDFFile::Write2DInt(hid_t group, const std::vector<int64_t>& x,
    int cols, const std::string& name) {
  hsize_t dims[2];
  dims[0]=x.size()/cols;
  dims[1]=cols;
  assert(x.size()%cols==0);
  hid_t x_dspace=H5Screate_simple(2, dims, NULL);
  hid_t x_id=H5Dcreate(group, name.c_str(), H5T_STD_I64LE,
      x_dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  herr_t x_write=H5Dwrite(x_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL,
      H5P_DEFAULT, &x[0]);
  H5Dclose(x_id);
  H5Sclose(x_dspace);
  return true;
}

hid_t HDFFile::ImageGroup() {
  if (image_group_>0) {
    return image_group_;
  }
  hid_t root=H5Gopen(file_id_, "/", H5P_DEFAULT);
  htri_t bImages=H5Lexists(root, "images", H5P_DEFAULT);
  if (bImages<0) {
    BOOST_LOG_TRIVIAL(error), "Could not find if /images exists in file.";
    return false;
  }
  hid_t image_group;
  if (bImages==0) {
    std::stringstream trajname;
    trajname<<"/images";
    BOOST_LOG_TRIVIAL(info)<<"Creating HDF5 directory "<<trajname.str();
    image_group_=H5Gcreate(file_id_, trajname.str().c_str(), H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);
  } else {
    image_group=H5Gopen(file_id_, "/images", H5P_DEFAULT);
  }
  H5Gclose(root);
  return image_group_;
}