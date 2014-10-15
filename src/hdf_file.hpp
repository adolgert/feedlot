#ifndef _HDF_FILE_H_
#define _HDF_FILE_H_ 1

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <array>
#include "boost/program_options.hpp"
#include "trajectory.hpp"
#include "parameter.hpp"
#include "hdf5.h"
#include "smv.hpp"

class HDFFile {
  hid_t file_id_;
  hid_t trajectory_group_;
  std::string filename_;
  bool open_;
  mutable std::mutex single_writer_;
 public:
  using TrajectoryType=std::vector<TrajectoryEntry>;
  explicit HDFFile(const std::string& filename);
  HDFFile(const HDFFile& o);
  HDFFile& operator=(const HDFFile&)=delete;
  ~HDFFile();
  bool Open(bool truncate=true);
  bool Close();
  template<typename Params>
  bool SaveTrajectory(const Params& params,
    int seed, int idx, const TrajectoryType& trajectory) const;
  template<typename Params>
  bool SavePenTrajectory(const Params& params,
    int seed, int idx, const std::vector<PenTrajectory>& trajectory) const;
  bool WriteExecutableData(const std::map<std::string,std::string>& compile,
    const boost::program_options::basic_parsed_options<char>& cmdline,
    const std::vector<int64_t>& initial_values) const;
};


template<typename Params>
bool HDFFile::SaveTrajectory(const Params& params,
    int seed, int idx, const TrajectoryType& trajectory) const {
  std::unique_lock<std::mutex> only_me(single_writer_);
  assert(open_);
  hsize_t dims[1];
  dims[0]=trajectory.size();
  hid_t dataspace_id=H5Screate_simple(1, dims, NULL);

  // Disk storage types are defined exactly.
  hid_t trajectory_type=H5Tcreate(H5T_COMPOUND, sizeof(TrajectoryEntry));
  H5Tinsert(trajectory_type, "s", HOFFSET(TrajectoryEntry, s),
      H5T_STD_I64LE);
  H5Tinsert(trajectory_type, "i", HOFFSET(TrajectoryEntry, i),
      H5T_STD_I64LE);
  H5Tinsert(trajectory_type, "r", HOFFSET(TrajectoryEntry, r),
      H5T_STD_I64LE);
  H5Tinsert(trajectory_type, "t", HOFFSET(TrajectoryEntry, t),
      H5T_IEEE_F64LE);
  if (trajectory_type<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not make HD5 type "<<trajectory_type;
    herr_t space_status=H5Sclose(dataspace_id);
    return false;
  }

  // When writing, ask library to translate from native type (H5T_NATIVE_LONG)
  // to disk storage type (H5T_STD_I64LE) if necessary.
  hid_t write_trajectory_type=H5Tcreate(H5T_COMPOUND, sizeof(TrajectoryEntry));
  H5Tinsert(write_trajectory_type, "s", HOFFSET(TrajectoryEntry, s),
      H5T_NATIVE_LONG);
  H5Tinsert(write_trajectory_type, "i", HOFFSET(TrajectoryEntry, i),
      H5T_NATIVE_LONG);
  H5Tinsert(write_trajectory_type, "r", HOFFSET(TrajectoryEntry, r),
      H5T_NATIVE_LONG);
  H5Tinsert(write_trajectory_type, "t", HOFFSET(TrajectoryEntry, t),
      H5T_NATIVE_DOUBLE);
  if (write_trajectory_type<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not make HD5 native type "
      <<write_trajectory_type;
    herr_t trajt_status=H5Tclose(trajectory_type);
    herr_t space_status=H5Sclose(dataspace_id);
    return false;
  }

  std::stringstream dset_name;
  dset_name << "dset" << seed << "-" << idx;
  hid_t dataset_id=H5Dcreate2(trajectory_group_, dset_name.str().c_str(),
    write_trajectory_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not make HD5 dataset "<<dataset_id;
    herr_t trajt_status=H5Tclose(trajectory_type);
    herr_t wtrajt_status=H5Tclose(write_trajectory_type);
    herr_t space_status=H5Sclose(dataspace_id);
    return false;
  }

  if (trajectory.size()>0) {
    herr_t write_status=H5Dwrite(dataset_id, write_trajectory_type,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &trajectory[0]);

    if (write_status<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write HD5 dataset "<<write_status;
      herr_t trajt_status=H5Tclose(trajectory_type);
      herr_t wtrajt_status=H5Tclose(write_trajectory_type);
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
    hid_t attr0_id=H5Acreate2(dataset_id, p.name.c_str(), H5T_IEEE_F64LE,
      dspace_id, H5P_DEFAULT, H5P_DEFAULT);
    herr_t atstatus=H5Awrite(attr0_id, H5T_NATIVE_DOUBLE, &p.value);
    if (atstatus<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write attribute "<<p.name;
    }
    H5Aclose(attr0_id);
  }
  H5Sclose(dspace_id);

  herr_t wtrajt_status=H5Tclose(write_trajectory_type);
  herr_t trajt_status=H5Tclose(trajectory_type);
  herr_t close_status=H5Dclose(dataset_id);
  herr_t space_status=H5Sclose(dataspace_id);
  return true;
}



template<typename Params>
bool HDFFile::SavePenTrajectory(const Params& params,
    int seed, int idx, const std::vector<PenTrajectory>& trajectory) const {
  std::unique_lock<std::mutex> only_me(single_writer_);
  assert(open_);
  hsize_t dims[1];
  dims[0]=trajectory.size();
  hid_t dataspace_id=H5Screate_simple(1, dims, NULL);

  // Disk storage types are defined exactly.
  hid_t trajectory_type=H5Tcreate(H5T_COMPOUND, sizeof(PenTrajectory));
  H5Tinsert(trajectory_type, "individual", HOFFSET(PenTrajectory, individual),
      H5T_STD_I64LE);
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
  H5Tinsert(write_trajectory_type, "individual", HOFFSET(PenTrajectory,
      individual), H5T_NATIVE_LONG);
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

  std::stringstream dset_name;
  dset_name << "dset" << seed << "-" << idx;
  hid_t dataset_id=H5Dcreate2(trajectory_group_, dset_name.str().c_str(),
    write_trajectory_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset_id<0) {
    BOOST_LOG_TRIVIAL(error)<<"Could not make HD5 dataset "<<dataset_id;
    herr_t trajt_status=H5Tclose(trajectory_type);
    herr_t wtrajt_status=H5Tclose(write_trajectory_type);
    herr_t space_status=H5Sclose(dataspace_id);
    return false;
  }

  if (trajectory.size()>0) {
    herr_t write_status=H5Dwrite(dataset_id, write_trajectory_type,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &trajectory[0]);

    if (write_status<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write HD5 dataset "<<write_status;
      herr_t trajt_status=H5Tclose(trajectory_type);
      herr_t wtrajt_status=H5Tclose(write_trajectory_type);
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
    hid_t attr0_id=H5Acreate2(dataset_id, p.name.c_str(), H5T_IEEE_F64LE,
      dspace_id, H5P_DEFAULT, H5P_DEFAULT);
    herr_t atstatus=H5Awrite(attr0_id, H5T_NATIVE_DOUBLE, &p.value);
    if (atstatus<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not write attribute "<<p.name;
    }
    H5Aclose(attr0_id);
  }
  H5Sclose(dspace_id);

  herr_t wtrajt_status=H5Tclose(write_trajectory_type);
  herr_t trajt_status=H5Tclose(trajectory_type);
  herr_t close_status=H5Dclose(dataset_id);
  herr_t space_status=H5Sclose(dataspace_id);
  return true;
}

#endif

