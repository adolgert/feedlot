#ifndef _PARAMETER_HPP_
#define _PARAMETER_HPP_ 1

struct Parameter {
  std::string name;
  double value;
  std::string description;
  Parameter(std::string n, double v, std::string desc)
  : name(n), value(v), description(desc) {}
  Parameter()=default;
  virtual ~Parameter() {}
};

template<typename EnumType>
struct TypedParameter : public Parameter {
  EnumType kind;
  TypedParameter(EnumType k, std::string n, double v, std::string desc)
  : Parameter(n, v, desc), kind(k) {}
  TypedParameter()=default;
  virtual ~TypedParameter() {}
};

#endif
