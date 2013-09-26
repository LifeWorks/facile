#include <octave/oct.h>

#include "steady_stateCVODE.h"

DEFUN_DLD (steady_stateCVODEOctWrapper, args, nargout,
	   "steady_stateCVODEOctWrapper(ivalues,rates,tv,events,options)\n") {

octave_stdout << "nargout=" << nargout << "\n";

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

  if (arg_ivalues.rows() != 1 || arg_ivalues.cols() != 3) {
    args_ok = 0;
    printf("ERROR: steady_stateCVODEOctWrapper() -- bad ivalues arg\n"); 
    print_usage();
    return retval;
  }

  if (arg_rates.rows() != 1 || arg_rates.cols() != 3) {
    args_ok = 0;
    printf("ERROR: steady_stateCVODEOctWrapper() -- bad rates arg\n"); 
    print_usage();
    return retval;
  }

  if (arg_tv.rows() != 1 || arg_tv.cols() < 2) {
    args_ok = 0;
    printf("ERROR: steady_stateCVODEOctWrapper() -- bad tv arg\n"); 
    print_usage();
    return retval;
  }

  if (arg_events.rows() > 1) {
    args_ok = 0;
    printf("ERROR: steady_stateCVODEOctWrapper() -- bad events arg\n"); 
    print_usage();
    return retval;
  }

  if (arg_options.rows() != 1 || arg_options.cols() != NUM_ODE_OPTIONS) {
    args_ok = 0;
    printf("ERROR: steady_stateCVODEOctWrapper() -- bad options arg\n"); 
    print_usage();
    return retval;
  }

  if (args_ok) {
    double ivalues[3];
    for(unsigned int i; i < 3; i++)
      ivalues[i]= arg_ivalues(0,i);

    double ode_rate_constants[3];
    for(unsigned int i; i < 3; i++)
      ode_rate_constants[i] = arg_rates(0,i);

    unsigned int tv_cols = arg_tv.cols();
    std::vector<realtype> tv(tv_cols, 0.0);
    for (unsigned int i = 0; i < tv_cols; i++) { 
      tv[i] = arg_tv(0,i);
    }

    unsigned int num_ode_events = arg_events.cols();
    std::vector<realtype>        ode_events(num_ode_events,0);
    for(unsigned int i = 0; i < num_ode_events; i++) {
      ode_events[i] = arg_events(i,0);
    }
    static std::vector<int>      event_flags;
    static std::vector<realtype> event_times;

    // options
    std::vector<realtype> ode_options(NUM_ODE_OPTIONS,-1);
    ode_options[ODE_OPTION_RELTOL]       = (arg_options(0,ODE_OPTION_RELTOL) < 0 ? 1e-3 : arg_options(0,ODE_OPTION_RELTOL));
    ode_options[ODE_OPTION_ABSTOL]       = (arg_options(0,ODE_OPTION_ABSTOL) < 0 ? 1e-6 : arg_options(0,ODE_OPTION_ABSTOL));
    ode_options[ODE_OPTION_INITSTEP]     = (arg_options(0,ODE_OPTION_INITSTEP) < 0 ? -1   : arg_options(0,ODE_OPTION_INITSTEP));
    ode_options[ODE_OPTION_MINSTEP]      = (arg_options(0,ODE_OPTION_MINSTEP) < 0 ? -1   : arg_options(0,ODE_OPTION_MINSTEP));
    ode_options[ODE_OPTION_MAXSTEP]      = (arg_options(0,ODE_OPTION_MAXSTEP) < 0 ? -1   : arg_options(0,ODE_OPTION_MAXSTEP));
    ode_options[ODE_OPTION_SS_TIMESCALE] = (arg_options(0,ODE_OPTION_SS_TIMESCALE) < 0 ? 1  : arg_options(0,ODE_OPTION_SS_TIMESCALE));
    ode_options[ODE_OPTION_SS_RELTOL]    = (arg_options(0,ODE_OPTION_SS_RELTOL) < 0 ? 1e-3  : arg_options(0,ODE_OPTION_SS_RELTOL));
    ode_options[ODE_OPTION_SS_ABSTOL]    = (arg_options(0,ODE_OPTION_SS_ABSTOL) < 0 ? 1e-6  : arg_options(0,ODE_OPTION_SS_ABSTOL));

    // allocate memory for solution (only once since static)
    static std::vector<realtype> aT;
    static std::vector<realtype> aY;

    // call solver
    int cvode_rval = cvode_sim_steady_state(ivalues, ode_rate_constants, tv, aT, aY,
                                 ode_options, ode_events, event_flags, event_times);

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

    Matrix m_event_flags(event_flags.size(), 1);
    Matrix m_event_times(event_times.size(), 1);

    for (int i=0; i < num_ode_events; i++) {
      m_event_flags(i,0) = event_flags[i];
      m_event_times(i,0) = event_times[i];
    }

    // copy results to return values
    retval(0) = mT;
    retval(1) = mY;
    retval(2) = m_event_flags;
    retval(3) = m_event_times;

    return retval;
  }
}
