#include"functions.h"
#include<cmath>
using namespace Eigen;

double fe_function(double a, std::string b, double time){
  double result;

  if(b=="RAMP"){
    result = a*(time/t_end);
  }

  if(b=="SIN"){
    result = a*sin(time);
  }

  if(b=="COS"){
    result = a*cos(time);
  }

  if(b=="STEP"){
    result = a;
  }

  return result;
}
