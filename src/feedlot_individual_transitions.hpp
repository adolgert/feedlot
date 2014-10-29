#ifndef _FEEDLOT_INDIVIDUAL_TRANSITIONS_HPP_
#define _FEEDLOT_INDIVIDUAL_TRANSITIONS_HPP_ 1


// Now make specific transitions.
template<typename BaseTransition>
class InfectiousExponential : public BaseTransition
{
public:
  virtual std::pair<bool, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>>
  Enabled(const typename BaseTransition::UserState& s,
    const typename BaseTransition::LocalMarking& lm,
    double te, double t0, typename BaseTransition::RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      double rate=I*s.params.at(SIRParam::Gamma);
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<
        afidd::smv::ExponentialDistribution<typename BaseTransition::RandGen>>(
        new afidd::smv::ExponentialDistribution<typename BaseTransition::RandGen>(rate, te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>(nullptr)};
    }
  }

  virtual void Fire(typename BaseTransition::UserState& s,
    typename BaseTransition::LocalMarking& lm, double t0,
      typename BaseTransition::RandGen& rng) override {
    lm.template Move<0, 0>(0, 2, 1); // Change the individual
    lm.template Move<0, 0>(1, 3, 1); // Change the summary count
  }
};

// Now make specific transitions.
template<typename BaseTransition>
class RecoverExponential : public BaseTransition
{
  virtual std::pair<bool, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>>
  Enabled(const typename BaseTransition::UserState& s,
    const typename BaseTransition::LocalMarking& lm,
    double te, double t0, typename BaseTransition::RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      double rate=I*s.params.at(SIRParam::Gamma);
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<afidd::smv::ExponentialDistribution<
        typename BaseTransition::RandGen>>(
        new afidd::smv::ExponentialDistribution<
        typename BaseTransition::RandGen>(rate, te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>(nullptr)};    }
  }

  virtual void Fire(typename BaseTransition::UserState& s,
    typename BaseTransition::LocalMarking& lm, double t0,
      typename BaseTransition::RandGen& rng) override {
    lm.template Move<0, 0>(0, 2, 1); // Move individual to summary count.
  }
};



// Now make specific transitions.
template<typename BaseTransition>
class Infectious : public BaseTransition
{
public:
  virtual std::pair<bool, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>>
  Enabled(const typename BaseTransition::UserState& s,
    const typename BaseTransition::LocalMarking& lm,
    double te, double t0, typename BaseTransition::RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<
        afidd::smv::WeibullDistribution<typename BaseTransition::RandGen>>(
        new afidd::smv::WeibullDistribution<typename BaseTransition::RandGen>(
          1/s.params.at(SIRParam::LatentBeta),
          s.params.at(SIRParam::LatentAlpha), te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>(nullptr)};
    }
  }

  virtual void Fire(typename BaseTransition::UserState& s,
    typename BaseTransition::LocalMarking& lm, double t0,
      typename BaseTransition::RandGen& rng) override {
    lm.template Move<0, 0>(0, 2, 1); // Change the individual
    lm.template Move<0, 0>(1, 3, 1); // Change the summary count
  }
};

// Now make specific transitions.
template<typename BaseTransition>
class Recover : public BaseTransition
{
  virtual std::pair<bool, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>>
  Enabled(const typename BaseTransition::UserState& s,
    const typename BaseTransition::LocalMarking& lm,
    double te, double t0, typename BaseTransition::RandGen& rng) override {
    int64_t I=lm.template Length<0>(0);
    if (I>0) {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<afidd::smv::GammaDistribution<
        typename BaseTransition::RandGen>>(
        new afidd::smv::GammaDistribution<
        typename BaseTransition::RandGen>(
          s.params.at(SIRParam::GammaAlpha),
          s.params.at(SIRParam::GammaBeta), te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>(nullptr)};    }
  }

  virtual void Fire(typename BaseTransition::UserState& s,
    typename BaseTransition::LocalMarking& lm, double t0,
      typename BaseTransition::RandGen& rng) override {
    lm.template Move<0, 0>(0, 2, 1); // Move individual to summary count.
    lm.template Move<0, 0>(1, 3, 1); // Change the summary count
  }
};


// Now make specific transitions.
template<typename BaseTransition>
class SubClinical : public BaseTransition
{
  virtual std::pair<bool, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>>
  Enabled(const typename BaseTransition::UserState& s,
    const typename BaseTransition::LocalMarking& lm,
    double te, double t0, typename BaseTransition::RandGen& rng) override {
    int64_t N=lm.template Length<0>(0);
    int64_t I=lm.template Length<0>(1);
    if (N>0 && I>0) {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate);
      return {true, std::unique_ptr<afidd::smv::GammaDistribution<
        typename BaseTransition::RandGen>>(
        new afidd::smv::GammaDistribution<
        typename BaseTransition::RandGen>(
          s.params.at(SIRParam::SubClinicalAlpha),
          s.params.at(SIRParam::SubClinicalBeta), te))};
    } else {
      //SMVLOG(BOOST_LOG_TRIVIAL(trace)<<"recover disable");
      return {false, std::unique_ptr<afidd::smv::TransitionDistribution< 
        typename BaseTransition::RandGen>>(nullptr)};
    }
  }

  virtual void Fire(typename BaseTransition::UserState& s,
    typename BaseTransition::LocalMarking& lm, double t0,
      typename BaseTransition::RandGen& rng) override {
    lm.template Move<0, 0>(0, 2, 1); // Move individual to summary count.
  }
};

#endif
