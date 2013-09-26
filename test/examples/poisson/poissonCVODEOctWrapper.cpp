#include <octave/oct.h>

#include "poissonCVODE.h"

DEFUN_DLD (poissonCVODEOctWrapper, args, nargout,
	   "poissonCVODEOctWrapper(ivalues,rates,tv,events,options)") {

  int args_ok = 1;
  octave_value_list retval;

  int nargin = args.length ();
  if (nargin != 5) {
    args_ok = 0;
    print_usage();
    return retval;
  }

  Matrix arg_ivalues = args(0).matrix_value();
  Matrix arg_rates = args(1).matrix_value();
  Matrix arg_tv = args(2).matrix_value();
  Matrix arg_events = args(3).matrix_value();
  Matrix arg_options = args(4).matrix_value();

  if (arg_ivalues.rows() != 1 || arg_ivalues.cols() != 2) {
    args_ok = 0;
    printf("ERROR: poissonCVODEOctWrapper() -- bad ivalues arg\n"); 
    print_usage();
    return retval;
  }

  if (arg_rates.rows() != 1 || arg_rates.cols() != 1) {
    args_ok = 0;
    printf("ERROR: poissonCVODEOctWrapper() -- bad rates arg\n"); 
    print_usage();
    return retval;
  }

  if (arg_tv.rows() != 1 || arg_tv.cols() < 2) {
    args_ok = 0;
    printf("ERROR: poissonCVODEOctWrapper() -- bad tv arg\n"); 
    print_usage();
    return retval;
  }

  if (arg_events.rows() > 1) {
    args_ok = 0;
    printf("ERROR: poissonCVODEOctWrapper() -- bad events arg\n"); 
    print_usage();
    return retval;
  }

  if (arg_options.rows() != 1 || arg_options.cols() != 5) {
    args_ok = 0;
    printf("ERROR: poissonCVODEOctWrapper() -- bad options arg\n"); 
    print_usage();
    return retval;
  }

  if (args_ok) {
    double ivalues[1];
    ivalues[0]= arg_ivalues(0,0);
    ivalues[1]= arg_ivalues(0,1);

    double rates[0];
    rates[0] = arg_rates(0,0);

    unsigned int tv_cols = arg_tv.cols();
    std::vector<realtype> tv(tv_cols, 0.0);
    for (unsigned int i = 0; i < tv_cols; i++) { 
      tv[i] = arg_tv(0,i);
    }

    std::vector<realtype> ode_events;
    ode_events.push_back(arg_events(0,0));
    ode_events.push_back(arg_events(0,1));

    // options
    std::vector<realtype> ode_options(5,-1);
    ode_options[ODE_OPTION_RELTOL]   = (arg_options(0,0) < 0 ? 1e-3 : arg_options(0,0));
    ode_options[ODE_OPTION_ABSTOL]   = (arg_options(0,1) < 0 ? 1e-6 : arg_options(0,1));
    ode_options[ODE_OPTION_INITSTEP] = (arg_options(0,2) < 0 ? -1   : arg_options(0,2));
    ode_options[ODE_OPTION_MINSTEP]  = (arg_options(0,3) < 0 ? -1   : arg_options(0,3));
    ode_options[ODE_OPTION_MAXSTEP]  = (arg_options(0,4) < 0 ? -1   : arg_options(0,4));

    // allocate memory for solution (only once since static)
    static std::vector<realtype> aT;
    static std::vector<realtype> aY;

    // call solver
    int cvode_rval = cvode_sim_poisson(ivalues, rates, tv, aT, aY, ode_options, ode_events);

    unsigned int aT_size = aT.size();
    unsigned int NEQ = aY.size() / aT_size;
    Matrix mT(aT_size, 1);
    Matrix mY(aT_size, NEQ);

    for (int i=0; i < aT_size; i++) {
      mT(i,0) = aT[i];
      unsigned int offset = i * NEQ;
      for (int j=0; j < NEQ; j++) {
	mY(i,j) = aY[offset + j];
      }
    }

    // copy results to return values 
    retval(0) = mT;
    retval(1) = mY;

    return retval;
  }
}

