#ifndef _MODEL_OPTIONS_HPP_
#define _MODEL_OPTIONS_HPP_ 1
#include <map>

/*! Boolean choices when building the model.
 */
enum class ModelOptions {
  ExponentialTransitions,
  DoubleGamma,
  Rider,
  AllToAllInfection,
  PenReplace
};

std::map<ModelOptions,bool> model_options();

#endif
