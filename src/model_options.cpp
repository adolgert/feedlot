#include "model_options.hpp"


std::map<ModelOptions,bool> model_options() {
	std::map<ModelOptions,bool> opts;
	opts[ModelOptions::ExponentialTransitions]=false;
	opts[ModelOptions::DoubleGamma]=false;
	opts[ModelOptions::Rider]=false;
	opts[ModelOptions::AllToAllInfection]=false;
  return opts;
}
