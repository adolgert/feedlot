#ifndef _SIR_EXP_HPP_
#define _SIR_EXP_HPP_ 1
#include <map>
#include <random>
#include "boost/random/mersenne_twister.hpp"
#include "mt19937.hpp"
#include "rider_enums.hpp"
#include "parameter.hpp"
#include "trajectory.hpp"
#include "model_options.hpp"

using RandGen=afidd::rng::mt19937;
//using RandGen=boost::mt19937;

int64_t SEIR_run(double end_time, const std::vector<int64_t>& seir_cnt,
    const std::vector<TypedParameter<SIRParam>>& parameters,
    std::shared_ptr<PenTrajectoryObserver> observer,
    RandGen& rng, std::map<ModelOptions,bool> opts, int block_cnt, int row_cnt);

#endif
