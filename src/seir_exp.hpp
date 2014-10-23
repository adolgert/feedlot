#ifndef _SIR_EXP_HPP_
#define _SIR_EXP_HPP_ 1
#include <map>
#include <random>
#include "boost/random/mersenne_twister.hpp"
#include "mt19937.hpp"
#include "rider_enums.hpp"
#include "parameter.hpp"
#include "trajectory.hpp"
#include "pen.hpp"
#include "model_options.hpp"

using RandGen=afidd::rng::mt19937;
//using RandGen=boost::mt19937;

int64_t SEIR_run(double end_time, const std::vector<int64_t>& seir_cnt,
    const std::vector<TypedParameter<SIRParam>>& parameters,
    std::map<ModelOptions,bool> opts,
    const PenContactGraph& pen_graph,
    std::shared_ptr<PenTrajectoryObserver> observer,
    RandGen& rng);

#endif
