#ifndef _HDF_FILE_H_
#define _HDF_FILE_H_ 1

#include <memory>
#include <string>
#include <tuple>
#include <vector>
#include <map>
#include <array>
#include "boost/program_options.hpp"
#include "boost/uuid/uuid.hpp"
#include "boost/uuid/uuid_generators.hpp"
#include "boost/uuid/uuid_io.hpp"
#include "trajectory.hpp"
#include "parameter.hpp"
#include "hdf5.h"
#include "smv.hpp"

class HDFFile {
  hid_t file_id_;
  hid_t trajectory_group_;
  std::string filename_;
  bool open_;
  unsigned int comp_cnt_;
  mutable std::mutex single_writer_;
 public:
  using TrajectoryType=std::vector<TrajectoryEntry>;
  explicit HDFFile(const std::string& filename);
  HDFFile(const HDFFile& o);
  HDFFile& operator=(const HDFFile&)=delete;
  ~HDFFile();
  bool Open(bool truncate=true);
  bool OpenRead(bool readwrite=true);
  bool Close();

  template<typename Params>
  bool SaveTrajectory(const Params& params,
    int seed, int idx, const TrajectoryType& trajectory) const;

  template<typename Params>
  bool SavePenTrajectory(const Params& params,
    int seed, int idx, const std::vector<PenTrajectory>& trajectory,
    const std::vector<TrajectoryEntry>& initial, int64_t nanosecs) const;

  bool WriteExecutableData(const std::map<std::string,std::string>& compile,
    const boost::program_options::basic_parsed_options<char>& cmdline,
    const std::vector<int64_t>& initial_values) const;

  bool WriteEnsembleVariables(int64_t elapsed_ns) const;

  bool Save2DPDF(const std::vector<double>& interpolant,
    const std::vector<double>& x, const std::vector<double>& y,
    std::string name) const;

  // Reading.
  std::vector<std::string> Trajectories() const;
  std::tuple<std::vector<double>,std::vector<int64_t>> EndTimes() const;
  std::vector<int64_t> InitialValues() const;
  std::vector<int64_t> LoadInitialPen(const std::string dataset_name) const;
  std::array<int64_t,2> EventsInFile() const;
  bool LoadTrajectoryCounts(const std::string& name,
    std::vector<int64_t>& seirc, int64_t& cnt) const;
  bool LoadTrajectoryTimes(const std::string& name,
    std::vector<double>& time, int64_t& cnt) const;
  TrajectoryType LoadTrajectoryFromPens(const std::string dataset_name) const;
 private:
  template<typename TrajType>
  bool SaveTotalTimes(int seed, int idx,
    const TrajType& trajectory) const;

  bool WriteUUIDTo(hid_t group) const;
  bool Write1DFloat(hid_t group, const std::vector<double>& x,
    const std::string& name) const;
};


template<typename TrajType>
bool HDFFile::SaveTotalTimes(int seed, int idx,
    const TrajType& trajectory) const {
  // No mutex!
  hsize_t dims[1];
  dims[0]=trajectory.size();
  hid_t dataspace_id=H5Screate_simple(1, dims, NULL);

  std::stringstream dset_name;
  dset_name << "dset" << seed << "-" << idx << "/seirtotaltimes";
  hid_t dataset_id=H5Dcreate2(trajectory_group_, dset_name.str().c_str(),
    H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not make HD5 dataset "<<dataset_id;
    herr_t space_status=H5Sclose(dataspace_id);
    return false;
  }

  if (trajectory.size()>0) {
    std::vector<double> vals(dims[0]);
    for (size_t i=0; i<dims[0]; ++i) {
      vals[i]=trajectory[i].t;
    }
    herr_t write_status=H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &vals[0]);

    if (write_status<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write HD5 dataset "<<write_status;
      herr_t space_status=H5Sclose(dataspace_id);
      herr_t close_status=H5Dclose(dataset_id);
      return false;
    }
  } else {
    BOOST_LOG_TRIVIAL(warning)<<"There was no trajectory to write.";
  }

  herr_t close_status=H5Dclose(dataset_id);
  herr_t space_status=H5Sclose(dataspace_id);
  return true;
}


template<typename Params>
bool HDFFile::SaveTrajectory(const Params& params,
    int seed, int idx, const TrajectoryType& trajectory) const {
  std::unique_lock<std::mutex> only_me(single_writer_);
  assert(open_);
  hsize_t dims[2];
  dims[0]=trajectory.size();
  dims[1]=comp_cnt_;
  hid_t dataspace_id=H5Screate_simple(2, dims, NULL);

  std::stringstream dset_name;
  dset_name << "dset" << seed << "-" << idx << "/seirtotal";
  hid_t dataset_id=H5Dcreate2(trajectory_group_, dset_name.str().c_str(),
    H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not make HD5 dataset "<<dataset_id;
    herr_t space_status=H5Sclose(dataspace_id);
    return false;
  }

  if (trajectory.size()>0) {
    int comp_cnt_=5;
    std::vector<int64_t> vals(dims[0]*dims[1]);
    for (size_t i=0; i<dims[0]; ++i) {
      vals[comp_cnt_*i+0]=trajectory[i].s;
      vals[comp_cnt_*i+1]=trajectory[i].e;
      vals[comp_cnt_*i+2]=trajectory[i].i;
      vals[comp_cnt_*i+3]=trajectory[i].r;
      vals[comp_cnt_*i+4]=trajectory[i].c;
    }
    herr_t write_status=H5Dwrite(dataset_id, H5T_STD_I64LE,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &vals[0]);

    if (write_status<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write HD5 dataset "<<write_status;
      herr_t space_status=H5Sclose(dataspace_id);
      herr_t close_status=H5Dclose(dataset_id);
      return false;
    }
  } else {
    BOOST_LOG_TRIVIAL(warning)<<"There was no trajectory to write.";
  }
  // Now write dataset attributes.
  hsize_t adims=1;
  hid_t dspace_id=H5Screate_simple(1, &adims, NULL);
  for (auto& p : params) {
    hid_t attr0_id=H5Acreate2(dataset_id, p.second.name.c_str(), H5T_IEEE_F64LE,
      dspace_id, H5P_DEFAULT, H5P_DEFAULT);
    herr_t atstatus=H5Awrite(attr0_id, H5T_NATIVE_DOUBLE, &p.second.value);
    if (atstatus<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write attribute "<<p.second.name;
    }
    H5Aclose(attr0_id);
  }
  H5Sclose(dspace_id);

  herr_t close_status=H5Dclose(dataset_id);
  herr_t space_status=H5Sclose(dataspace_id);

  this->SaveTotalTimes(seed, idx, trajectory);
  return true;
}



template<typename Params>
bool HDFFile::SavePenTrajectory(const Params& params,
    int seed, int idx, const std::vector<PenTrajectory>& trajectory,
    const std::vector<TrajectoryEntry>& initial, int64_t nanosecs) const {
  std::unique_lock<std::mutex> only_me(single_writer_);
  BOOST_LOG_TRIVIAL(debug)<<"Writing pen trajectory "<<open_;
  assert(open_);
  hsize_t dims[1];
  dims[0]=trajectory.size();
  hid_t dataspace_id=H5Screate_simple(1, dims, NULL);
  if (dataspace_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not create the dataspace";
  }

  // Disk storage types are defined exactly.
  hid_t trajectory_type=H5Tcreate(H5T_COMPOUND, sizeof(PenTrajectory));
  if (trajectory_type<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not create the dataset format datatype";
  }
  herr_t h5tres=H5Tinsert(trajectory_type, "individual",
      HOFFSET(PenTrajectory, individual), H5T_STD_I64LE);
  if (h5tres<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not insert datatype into dset trajectory";
  }
  H5Tinsert(trajectory_type, "pen", HOFFSET(PenTrajectory, pen),
      H5T_STD_I64LE);
  H5Tinsert(trajectory_type, "transition", HOFFSET(PenTrajectory, transition),
      H5T_STD_I64LE);
  H5Tinsert(trajectory_type, "time", HOFFSET(PenTrajectory, time),
      H5T_IEEE_F64LE);
  if (trajectory_type<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not make HD5 type "<<trajectory_type;
    herr_t space_status=H5Sclose(dataspace_id);
    return false;
  }

  // When writing, ask library to translate from native type (H5T_NATIVE_LONG)
  // to disk storage type (H5T_STD_I64LE) if necessary.
  hid_t write_trajectory_type=H5Tcreate(H5T_COMPOUND, sizeof(PenTrajectory));
  if (write_trajectory_type<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not create trajectory type for writing";
  }
  herr_t h5tres2=H5Tinsert(write_trajectory_type, "individual",
      HOFFSET(PenTrajectory, individual), H5T_NATIVE_LONG);
  if (h5tres2<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not insert datatype into dset trajectory";
  }
  H5Tinsert(write_trajectory_type, "pen", HOFFSET(PenTrajectory, pen),
      H5T_NATIVE_LONG);
  H5Tinsert(write_trajectory_type, "transition", HOFFSET(PenTrajectory,
      transition), H5T_NATIVE_LONG);
  H5Tinsert(write_trajectory_type, "time", HOFFSET(PenTrajectory, time),
      H5T_NATIVE_DOUBLE);
  if (write_trajectory_type<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not make HD5 native type "
      <<write_trajectory_type;
    herr_t trajt_status=H5Tclose(trajectory_type);
    herr_t space_status=H5Sclose(dataspace_id);
    return false;
  }

  BOOST_LOG_TRIVIAL(debug)<<"Creating group";
  std::stringstream dset_group_name;
  dset_group_name << "dset" << seed << "-" << idx;
  hid_t dataset_group_id=H5Gcreate(trajectory_group_,
    dset_group_name.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_group_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not create trajectory dataset group";
  }
  BOOST_LOG_TRIVIAL(debug)<<"Creating trajectory dataset";
  hid_t dataset_id=H5Dcreate2(dataset_group_id, "trajectory",
    write_trajectory_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not make HD5 dataset "<<dataset_id;
    herr_t dg_status=H5Gclose(dataset_group_id);
    herr_t trajt_status=H5Tclose(trajectory_type);
    herr_t wtrajt_status=H5Tclose(write_trajectory_type);
    herr_t space_status=H5Sclose(dataspace_id);
    return false;
  }

  BOOST_LOG_TRIVIAL(debug)<<"Writing trajectory dataset "<<trajectory.size();
  if (trajectory.size()>0) {
    herr_t write_status=H5Dwrite(dataset_id, write_trajectory_type,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &trajectory[0]);

    if (write_status<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write HD5 dataset "<<write_status;
      herr_t trajt_status=H5Tclose(trajectory_type);
      herr_t dg_status=H5Gclose(dataset_group_id);
      herr_t wtrajt_status=H5Tclose(write_trajectory_type);
      herr_t space_status=H5Sclose(dataspace_id);
      herr_t close_status=H5Dclose(dataset_id);
      return false;
    }
  } else {
    BOOST_LOG_TRIVIAL(warning)<<"There was no trajectory to write.";
  }
  // Now write dataset attributes.
  BOOST_LOG_TRIVIAL(debug)<<"Create parameter attributes";
  hsize_t adims=1;
  hid_t dspace_id=H5Screate_simple(1, &adims, NULL);
  for (auto& p : params) {
    hid_t attr0_id=H5Acreate2(dataset_group_id, p.second.name.c_str(),
      H5T_IEEE_F64LE, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
    if (attr0_id<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not create attribute for "<<p.second.name;
    }
    herr_t atstatus=H5Awrite(attr0_id, H5T_NATIVE_DOUBLE, &p.second.value);
    if (atstatus<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write attribute "<<p.second.name;
    }
    H5Aclose(attr0_id);
  }

  {
    hid_t attr0_id=H5Acreate2(dataset_group_id, "elapsed_time", H5T_STD_I64LE,
      dspace_id, H5P_DEFAULT, H5P_DEFAULT);
    if (attr0_id<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not create attribute for elapsed time";
    }
    herr_t atstatus=H5Awrite(attr0_id, H5T_NATIVE_DOUBLE, &nanosecs);
    if (atstatus<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write attribute elapsed time";
    }
    H5Aclose(attr0_id);
  }

  WriteUUIDTo(dataset_group_id);

  H5Sclose(dspace_id);

  // This is the set of initial pen values.
  std::vector<hsize_t> idims{initial.size(),comp_cnt_};
  std::vector<int64_t> initial_data(initial.size()*comp_cnt_);
  size_t eidx{0};
  for (auto entry : initial) {
    initial_data[comp_cnt_*eidx]=initial[eidx].s;
    initial_data[comp_cnt_*eidx+1]=initial[eidx].e;
    initial_data[comp_cnt_*eidx+2]=initial[eidx].i;
    initial_data[comp_cnt_*eidx+3]=initial[eidx].r;
    initial_data[comp_cnt_*eidx+4]=initial[eidx].c;
    ++eidx;
  }
  BOOST_LOG_TRIVIAL(debug)<<"Create initial value dataspace";
  hid_t pen_dspace_id=H5Screate_simple(2, &idims[0], NULL);
  if (pen_dspace_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not create simple dspace for initial vals";
  }
  BOOST_LOG_TRIVIAL(debug)<<"Create initial value dataset";
  hid_t pen_dset_id=H5Dcreate2(dataset_group_id, "initialpencount",
    H5T_STD_I64LE, pen_dspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (pen_dset_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not create dataset for initial vals.";
  }
  BOOST_LOG_TRIVIAL(debug)<<"Writing initial value dataset";
  herr_t pends_write=H5Dwrite(pen_dset_id, H5T_STD_I64LE, H5S_ALL, H5S_ALL,
    H5P_DEFAULT, &initial_data[0]);
  if (pends_write<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not write dataset for initial vals.";
  }

  herr_t pen_dset_status=H5Dclose(pen_dset_id);
  herr_t idspace_status=H5Sclose(pen_dspace_id);

  herr_t wtrajt_status=H5Tclose(write_trajectory_type);
  herr_t trajt_status=H5Tclose(trajectory_type);
  herr_t close_status=H5Dclose(dataset_id);
  herr_t space_status=H5Sclose(dataspace_id);
  herr_t dg_status=H5Gclose(dataset_group_id);
  return true;
}


#endif

