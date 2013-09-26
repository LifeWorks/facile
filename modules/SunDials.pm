###############################################################################
# File:     SunDials.pm
#
# Copyright (C) 2005-2012 Ollivier, Siso-Nadal, Swain et al.
#
# This program comes with ABSOLUTELY NO WARRANTY.
# This is free software, and you are welcome to redistribute it
# under certain conditions. See file LICENSE.TXT for details.
#
# Synopsys: Export model to a C for integration using CVODE.
##############################################################################
# Detailed Description:
# ---------------------
#
###############################################################################

package SunDials;

use strict;
use diagnostics;          # equivalent to -w command-line switch
use warnings;

use vars qw(@ISA @EXPORT);

require Exporter;
@ISA = qw(Exporter);

# these subroutines are imported when you use Matlab.pm
@EXPORT = qw(
	     export_cvode_files
	    );

use Globals;
use Expression;
use Utils;

# CHECK LIBRARY TO SEE FOR JAC FCN TYPEDEF
#/usr/include/cvode/cvode_direct.h
#typedef int (*CVDlsDenseJacFn)(??long?? int N, realtype t
my @typedef_location = `find /usr/include /usr/local/include -name 'cvode_direct.h'`;
my $typedef_str = "";
my $int_N_str = "long int N"; # default typedef if no include file found
if (@typedef_location == 0) {
    #print "Warning: could not find cvode_direct.h -- assuming recent sundials lib\n";
}
if (@typedef_location == 1) {
    if (@typedef_location >= 2) {
	#print "Warning: found multiple versions of cvode_direct.h -- using first one at $typedef_location[0]\n";
    }
    chomp($typedef_location[0]);
    $typedef_str = `grep -P 'typedef.*CVDlsDenseJacFn.*int N' $typedef_location[0]`;
    chomp($typedef_str);
    if ($typedef_str =~ /typedef.*CVDlsDenseJacFn.*int N/) {
	#print "XXX: $typedef_location[0]\n";
	$int_N_str = ($typedef_str =~ /long\s+int\s+N/) ? "long int N" : "int N";
    }
}

# Prints python script for ODE integration using SciPy
sub export_cvode_files {
    my $self = shift;

    my %args = (
	output_file_prefix => undef,
	split_flag => 0,
	extern_flag => 0,
	jacobian_flag => 0,
	factor_flag => 0,
	@_,
       );

    my $output_file_prefix = $args{output_file_prefix};
    my $split_flag = $args{split_flag};
    my $extern_flag = $args{extern_flag};
    my $jacobian_flag = $args{jacobian_flag};
    my $factor_flag = $args{factor_flag};

    my $node_list_ref = $self->get_node_list();
    my $variable_list_ref = $self->get_variable_list();

    my $driver_file_contents;
    my $cvode_header_file_contents;
    my $cvode_file_contents;
    my %ode_func_contents = ();
    my $jac_func_contents;

    my $octwrap_file_contents;

    # get configuration vars affecting output
    my $config_ref = $self->get_config_ref();

    my $output_type = $config_ref->{output_type};

    my $compartment_volume = $config_ref->{compartment_volume};
    my $tf = $config_ref->{tf};
    my $tv = $config_ref->{tv};
    my $tk = $config_ref->{tk};
    my $solver = $config_ref->{"${output_type}_ode_solver"};
    my $solver_options = $config_ref->{"${output_type}_solver_options"}{$solver};
    my $ode_event_times = $config_ref->{ode_event_times};
    my $SS_timescale = $config_ref->{SS_timescale};
    my $SS_RelTol = $config_ref->{SS_RelTol};
    my $SS_AbsTol = $config_ref->{SS_AbsTol};
    my $plot_flag = $config_ref->{plot_flag};

    my @ode_events = defined $ode_event_times ? split(/\s*[, ]\s*/,$ode_event_times) : ();
    my $num_ode_events = @ode_events;

    # construct lists of nodes and variables to print out
    my @free_nodes = $node_list_ref->get_ordered_free_node_list();
    my @free_node_names = map ($_->get_name(), @free_nodes);
    my $num_free_nodes = @free_nodes;

    my @variables = $variable_list_ref->get_list();

    my @constant_rate_params = grep (($_->get_type() =~ /^(rate)|(other)$/ &&
				      $_->get_is_expression_flag() == 0), @variables);
    my @constant_rate_param_names = map($_->get_name(), @constant_rate_params);

    my @constant_rate_expressions = grep (($_->get_type() =~ /^(rate)|(other)$/ &&
					   $_->get_is_expression_flag() == 1 &&
					   $_->get_is_dynamic_flag() == 0), @variables);
    my @constant_rate_expression_names = map($_->get_name(), @constant_rate_expressions);

    my @dynamic_rate_expressions = grep ($_->get_type() =~ /^(rate)|(other)$/ &&
				 $_->get_is_dynamic_flag() == 1, @variables);

    my @ode_rate_constants = (@constant_rate_params,@constant_rate_expressions);
    my @ode_rate_constants_names = map($_->get_name(), @ode_rate_constants);
    my $num_ode_rate_constants = @ode_rate_constants;

    my @moiety_totals = grep ($_->get_type() eq "moiety_total", @variables);
    my @constrained_node_expressions = grep ($_->get_type() eq "constrained_node_expression",
					     @variables);

    #********************************************************#
    # generate cvode header                                  #
    #********************************************************#	
    $cvode_header_file_contents .= <<CVODE;
//
// File generated by Facile version $VERSION
//

#include <cvode/cvode.h>             // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h>  // serial N_Vector types, fcts., macros
#include <cvode/cvode_dense.h>       // prototype for CVDense
#include <sundials/sundials_dense.h> // definitions DlsMat DENSE_ELEM
#include <sundials/sundials_types.h> // definition of type realtype

#include <vector>

#define ODE_OPTION_RELTOL        0
#define ODE_OPTION_ABSTOL        1
#define ODE_OPTION_INITSTEP      2
#define ODE_OPTION_MINSTEP       3
#define ODE_OPTION_MAXSTEP       4
#define ODE_OPTION_SS_TIMESCALE  5
#define ODE_OPTION_SS_RELTOL     6
#define ODE_OPTION_SS_ABSTOL     7
#define NUM_ODE_OPTIONS          8

static const unsigned int NEQ = $num_free_nodes;

int cvode_sim_$output_file_prefix(double *ivalues, double *ode_rate_constants,
                                  std::vector<realtype> &tv,
                                  std::vector<realtype> &aT,
                                  std::vector<realtype> &aY,
                                  std::vector<realtype> &ode_options,
                                  std::vector<realtype> &ode_events,
                                  std::vector<int>      &event_flags,
                                  std::vector<realtype> &event_times
                                  );

int cvode_sim_print_output(std::vector<realtype> &aT,
                           std::vector<realtype> &aY);

CVODE

    #********************************************************#
    # generate Octave wrapper file                           #
    #********************************************************#	
    $octwrap_file_contents .= <<WRAP;
#include <octave/oct.h>

#include "${output_file_prefix}CVODE.h"

DEFUN_DLD (${output_file_prefix}CVODEOctWrapper, args, nargout,
	   "${output_file_prefix}CVODEOctWrapper(ivalues,rates,tv,events,options)\\n") {

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

  if (arg_ivalues.rows() != 1 || arg_ivalues.cols() != $num_free_nodes) {
    args_ok = 0;
    printf("ERROR: ${output_file_prefix}CVODEOctWrapper() -- bad ivalues arg\\n"); 
    print_usage();
    return retval;
  }

  if (arg_rates.rows() != 1 || arg_rates.cols() != $num_ode_rate_constants) {
    args_ok = 0;
    printf("ERROR: ${output_file_prefix}CVODEOctWrapper() -- bad rates arg\\n"); 
    print_usage();
    return retval;
  }

  if (arg_tv.rows() != 1 || arg_tv.cols() < 2) {
    args_ok = 0;
    printf("ERROR: ${output_file_prefix}CVODEOctWrapper() -- bad tv arg\\n"); 
    print_usage();
    return retval;
  }

  if (arg_events.rows() > 1) {
    args_ok = 0;
    printf("ERROR: ${output_file_prefix}CVODEOctWrapper() -- bad events arg\\n"); 
    print_usage();
    return retval;
  }

  if (arg_options.rows() != 1 || arg_options.cols() != NUM_ODE_OPTIONS) {
    args_ok = 0;
    printf("ERROR: ${output_file_prefix}CVODEOctWrapper() -- bad options arg\\n"); 
    print_usage();
    return retval;
  }

  if (args_ok) {
    double ivalues[$num_free_nodes];
    for(unsigned int i; i < $num_free_nodes; i++)
      ivalues[i]= arg_ivalues(0,i);

    double ode_rate_constants[$num_ode_rate_constants];
    for(unsigned int i; i < $num_ode_rate_constants; i++)
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
    int cvode_rval = cvode_sim_${output_file_prefix}(ivalues, ode_rate_constants, tv, aT, aY,
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
WRAP

    #********************************************************#
    # generate C/C++ driver file                             #
    #********************************************************#	
    $driver_file_contents .= <<DRIVER;
//
// File generated by Facile version $VERSION
//

#include <stdio.h>
#include <math.h>

// Custom header file
#include "${output_file_prefix}CVODE.h"

int main() {

  double ivalues[NEQ];
  double ode_rate_constants[$num_ode_rate_constants];

DRIVER

    # initial values
    $driver_file_contents .= "  // initial values (free nodes only)\n";
    for (my $i=0; $i < @free_nodes; $i++) {
      my $node_ref = $free_nodes[$i];
      my $node_name = $node_ref->get_name();
      my $initial_value = $node_ref->get_initial_value_in_molarity($compartment_volume);
      $driver_file_contents .= "  ivalues[$i] = $initial_value; // $node_name\n";
    }
    $driver_file_contents .= "\n";

    # rate constants
    $driver_file_contents .= "\n  // rate constants and constant expressions\n";
    for (my $i=0; $i < @ode_rate_constants; $i++) {
      my $variable_ref  = $ode_rate_constants[$i];
      my $variable_name = $variable_ref->get_name();
      my $variable_value = $variable_ref->get_value();
      $driver_file_contents .= sprintf("  realtype %-16s = ode_rate_constants[$i] = $variable_value;\n",$variable_name);
    }
    $driver_file_contents .= "\n";

    $driver_file_contents .= "  // time interval and sample times\n";
    $driver_file_contents .= "  realtype t0 = 0;\n";
    $driver_file_contents .= "  realtype tf = $tf;\n\n";
    $driver_file_contents .= "  std::vector<realtype> tv;\n";
    my ($tv_t0,$tv_dt,$tv_tf);
    if ($tv =~ /\[\s*(\S+)\s*:\s*(\S+)\s*:\s*(\S+)\]/) {
      $tv_t0 = $1;
      $tv_dt = $2;
      $tv_tf = $3;
    } elsif ($tv =~ /\[\s*(\S+)\s*(\S+)\]/) {
      $tv_t0 = $1;
      $tv_tf = $2;
      $tv_dt = "($tv_tf - $tv_t0)";
    }
    $driver_file_contents .= <<CVODE;
  realtype TVec_t0 = $tv_t0;
  realtype TVec_tf = $tv_tf;
  realtype TVec_dt = $tv_dt;
  unsigned int num_samples = ceil(TVec_tf / TVec_dt) + 1;
  for(unsigned int i=0; i < num_samples; i++) {
    realtype next_value = TVec_t0 + realtype(i) * TVec_dt;
    if (next_value <= TVec_tf) tv.push_back(next_value);
  };

CVODE

    $driver_file_contents .= "  // ode events\n";
    $driver_file_contents .= "  std::vector<realtype> ode_events;\n";
    $driver_file_contents .= join "", map {"  ode_events.push_back($_);\n"} (map {$_ eq '~' ? 0.0 : $_} @ode_events);
    $driver_file_contents .= "  unsigned int num_ode_events = ode_events.size();\n";
    $driver_file_contents .= "  std::vector<int>      event_flags;\n";
    $driver_file_contents .= "  std::vector<realtype> event_times;\n";
    $driver_file_contents .= "\n";

    $driver_file_contents .= "  // options\n";
    $driver_file_contents .= "  std::vector<realtype> ode_options(NUM_ODE_OPTIONS,-1);\n";
    $solver_options->{reltol} = 1.0e-3 if !defined $solver_options->{reltol};
    $solver_options->{abstol} = 1.0e-6 if !defined $solver_options->{abstol};
    $solver_options->{CVodeSetInitStep} = -1 if !defined $solver_options->{CVodeSetInitStep};
    $solver_options->{CVodeSetMinStep} = -1 if !defined $solver_options->{CVodeSetMinStep};
    $solver_options->{CVodeSetMaxStep} = -1 if !defined $solver_options->{CVodeSetMaxStep};
    $driver_file_contents .= sprintf("  ode_options[ODE_OPTION_RELTOL] = %e;\n",$solver_options->{reltol});
    $driver_file_contents .= sprintf("  ode_options[ODE_OPTION_ABSTOL] = %e;\n",$solver_options->{abstol});
    $driver_file_contents .= sprintf("  ode_options[ODE_OPTION_INITSTEP] = %g;\n",$solver_options->{CVodeSetInitStep});
    $driver_file_contents .= sprintf("  ode_options[ODE_OPTION_MINSTEP] = %g;\n",$solver_options->{CVodeSetMinStep});
    $driver_file_contents .= sprintf("  ode_options[ODE_OPTION_MAXSTEP] = %g;\n",$solver_options->{CVodeSetMaxStep});
    $driver_file_contents .= sprintf("  ode_options[ODE_OPTION_SS_TIMESCALE] = %g;\n",$SS_timescale);
    $driver_file_contents .= sprintf("  ode_options[ODE_OPTION_SS_RELTOL] = %g;\n",$SS_RelTol);
    $driver_file_contents .= sprintf("  ode_options[ODE_OPTION_SS_ABSTOL] = %g;\n",$SS_AbsTol);
    $driver_file_contents .= "\n";

    $driver_file_contents .= "  // call solver\n";
    $driver_file_contents .= "  std::vector<realtype> aT;\n";
    $driver_file_contents .= "  std::vector<realtype> aY;\n\n";
    $driver_file_contents .= "  cvode_sim_$output_file_prefix(ivalues, ode_rate_constants, tv, aT, aY,\n";
    $driver_file_contents .= "            ode_options, ode_events, event_flags, event_times);\n";
    $driver_file_contents .= "\n";

    $driver_file_contents .= "  // print output to screen\n";
    $driver_file_contents .= "  cvode_sim_print_output(aT, aY);\n";
    $driver_file_contents .= "\n";

    $driver_file_contents .= "  return (0);\n";
    $driver_file_contents .= "}\n";

    #********************************************************#
    # generate cvode solver file                             #
    #********************************************************#	
    $cvode_file_contents .= <<CVODE;
//
// File generated by Facile version $VERSION
//
// This code is adapted from the file "cvRoberts_dns.c", a usage example
// that is part of the the of the Sundials suite:
// * -----------------------------------------------------------------
// * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
// *                Radu Serban @ LLNL
// * -----------------------------------------------------------------
//

#include <stdio.h>
#include <math.h>

// Custom header file
#include "${output_file_prefix}CVODE.h"

#define Ith(v,i)    NV_Ith_S(v,i)
#define IJth(A,i,j) DENSE_ELEM(A,i,j)

// These macros needed to properly handle discontinuities
#if defined(SUNDIALS_SINGLE_PRECISION)
#define DT_BEFORE(x) nextafterf((float)x,x/2.0)
#define DT_AFTER(x) nextafterf((float)x,x*2.0)
#elif defined(SUNDIALS_DOUBLE_PRECISION)
#define DT_BEFORE(x) nextafter((double)x,x/2.0)
#define DT_AFTER(x) nextafter((double)x,x*2.0)
#elif defined(SUNDIALS_EXTENDED_PRECISION)
#define DT_BEFORE(x) nextafterl((long double)x,x/2.0)
#define DT_AFTER(x) nextafterl((long double)x,x*2.0)
#endif

// Private function to check function return values

static int check_flag(void *flagvalue, char *funcname, int opt);

// Private function to setup events
static int setup_next_event(void *cvode_mem, realtype t, realtype tf,
                            unsigned int &ode_event_index,
                            unsigned int num_ode_events,
                            std::vector<realtype> &ode_events,
                            std::vector<realtype> &ode_events_minus_dt,
                            int &steady_state_flag, realtype &ss_timeout
                           );

// Private functions to output results

static void StoreOutput(std::vector<realtype> &T,
                        std::vector<realtype> &Y,
                        realtype t, N_Vector y);

static void PrintOutput(realtype t, N_Vector y);

// Private function to print final statistics

static void PrintFinalStats(void *cvode_mem);

// User-defined model functions by the Solver

static int f_dydt(realtype t, N_Vector y, N_Vector ydot, void * user_data);

static int f_jac($int_N_str, realtype t,
                 N_Vector y, N_Vector fy, DlsMat J, void * user_data,
                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

// Constants and functions that may appear in user-defined model functions
static const realtype pi = atan(1.0)*4.0;

// Octave-style square() function -- (duty cycle between 0 and 1)
static realtype square(realtype t, realtype duty) {
  realtype r;

  if(t < 0) {
    realtype shift = ceil(-t / (2.0 * pi));
    t = t + shift * 2 * pi;
  }

  realtype phase = fmod(t, 2*pi);

  if (phase < 2*pi * duty) {
    r = +1;
  } else {
    r = -1;
  }
  return r;
}

static realtype square(realtype t) {
  realtype r;
  realtype duty = 0.5;

  if(t < 0) {
    realtype shift = ceil(-t / (2.0 * pi));
    t = t + shift * 2 * pi;
  }

  realtype phase = fmod(t, 2*pi);

  if (phase < 2*pi * duty) {
    r = +1;
  } else {
    r = -1;
  }
  return r;
}


// Main Program

int cvode_sim_$output_file_prefix(double *ivalues, double *ode_rate_constants,
                                  std::vector<realtype> &tv,
                                  std::vector<realtype> &aT,
                                  std::vector<realtype> &aY,
                                  std::vector<realtype> &ode_options,
                                  std::vector<realtype> &ode_events,
                                  std::vector<int>      &event_flags,
                                  std::vector<realtype> &event_times
                                 ) {

  if (tv.size() < 2) {
    printf("ERROR: invalid time vector (tv)\\n");
    printf("==> tv.size() = : %d\\n",tv.size());
    return(1);
  }

  realtype t0 = tv[0];
  N_Vector y = NULL;

  void *cvode_mem = NULL;

  int flag;

  // Create serial vector of length NEQ for I.C.
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, (char*)"N_VNew_Serial", 0)) return(1);

CVODE

    # initial values
    $cvode_file_contents .= "  // initial values (free nodes only)\n";
    for (my $i=0; $i < @free_nodes; $i++) {
      my $node_ref = $free_nodes[$i];
      my $node_name = $node_ref->get_name();
      my $initial_value = $node_ref->get_initial_value_in_molarity($compartment_volume);
      $cvode_file_contents .= "  Ith(y,$i) = ivalues[$i]; // $node_name\n";
    }
    $cvode_file_contents .= "\n";

    # tolerances and steady-state options
    $cvode_file_contents .= <<CVODE;
  // set absolute and relative integration tolerances
  realtype reltol = RCONST(ode_options[ODE_OPTION_RELTOL]);
  N_Vector abstol = N_VNew_Serial(NEQ);
  if (check_flag((void *)abstol, (char*)"N_VNew_Serial", 0)) return(1);
  for (int i=0; i<NEQ; i++) {
    Ith(abstol,i) = RCONST(ode_options[ODE_OPTION_ABSTOL]);
  }

  // set steady-state timescale and tolerances
  realtype SS_timescale = ode_options[ODE_OPTION_SS_TIMESCALE];
  realtype SS_RelTol    = ode_options[ODE_OPTION_SS_RELTOL];
  realtype SS_AbsTol    = ode_options[ODE_OPTION_SS_ABSTOL];

CVODE

    # setup and call ODE function
    $cvode_file_contents .= <<CVODE;

  // Call CVodeCreate to create the solver memory and specify the
  // Backward Differentiation Formula and the use of a Newton iteration
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, (char*)"CVodeCreate", 0)) return(1);
  
  // Call CVodeInit to initialize the integrator memory and specify the
  // user's right hand side function in y'=f(t,y), the inital time t0, and
  // the initial dependent variable vector y.
  flag = CVodeInit(cvode_mem, f_dydt, t0, y);
  if (check_flag(&flag, (char*)"CVodeInit", 1)) return(1);

  // Call CVodeSVtolerances to specify the scalar relative tolerance
  // and vector absolute tolerances
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, (char*)"CVodeSVtolerances", 1)) return(1);

  // Call CVodeSetUserData to specify rate constants
  void * user_data[3];
  user_data[0] = ode_rate_constants;
  user_data[1] = &event_flags;
  user_data[2] = &event_times;
  flag = CVodeSetUserData(cvode_mem, user_data);
  if (check_flag(&flag, (char*)"CVodeSetUserData", 1)) return(1);

  // Call CVDense to specify the CVDENSE dense linear solver
  flag = CVDense(cvode_mem, NEQ);
  if (check_flag(&flag, (char*)"CVDense", 1)) return(1);

CVODE

    my $jac_enable = $jacobian_flag ? "" : "//";

       $cvode_file_contents .= <<CVODE;
  // Set the Jacobian routine
$jac_enable  flag = CVDlsSetDenseJacFn(cvode_mem, f_jac);
$jac_enable  if (check_flag(&flag, (char*)"CVDlsSetDenseJacFn", 1)) return(1);

  // Set misc ODE solver options
  if (ode_options[ODE_OPTION_INITSTEP] >= 0.0) {
    flag = CVodeSetInitStep(cvode_mem, ode_options[ODE_OPTION_INITSTEP]);
    if (check_flag(&flag, (char*)"CVodeSetInitStep", 1)) return(1);
  }
  if (ode_options[ODE_OPTION_MINSTEP] >= 0.0) {
    flag = CVodeSetMinStep(cvode_mem, ode_options[ODE_OPTION_MINSTEP]);
    if (check_flag(&flag, (char*)"CVodeSetMinStep", 1)) return(1);
  }
  if (ode_options[ODE_OPTION_MAXSTEP] >= 0.0) {
    flag = CVodeSetMaxStep(cvode_mem, ode_options[ODE_OPTION_MAXSTEP]);
    if (check_flag(&flag, (char*)"CVodeSetMaxStep", 1)) return(1);
  }

  // Initialize ODE events vector
  unsigned int num_ode_events = ode_events.size();
  event_flags.resize(num_ode_events);
  event_times.resize(num_ode_events);
  std::vector<realtype> ode_events_minus_dt(num_ode_events);
  for(unsigned int i; i < num_ode_events; i++) {
    event_flags[i] = 0;
    event_times[i] = ode_events[i];
    if (ode_events[i] > 0)
      // absolute-time event (discontinuity)
      ode_events_minus_dt[i] = DT_BEFORE(ode_events[i]);
    else if (ode_events[i] < 0)
      // relative-time event
      ode_events_minus_dt[i] = 0;
    else
      // steady-state event
      ode_events_minus_dt[i] = 0;
  }

  // Setup mode and storage vectors
  unsigned int tv_size = tv.size();
  int cvode_mode = tv_size > 2 ? CV_NORMAL : CV_ONE_STEP;
  unsigned int num_samples = (cvode_mode == CV_NORMAL) ? tv_size : 100;
  aT.resize(0); aT.reserve(num_samples);
  aY.resize(0); aY.reserve((num_samples)*NEQ);

  // Initial and final time
  realtype t = t0;
  realtype tf = tv[tv_size-1];

  // Setup first event if any
  unsigned int ode_event_index = 0;
  int steady_state_flag = 0;
  realtype ss_timeout = 0;
  if (setup_next_event(cvode_mem, t0, tf,
                       ode_event_index,num_ode_events,
                       ode_events,ode_events_minus_dt,
                       steady_state_flag,ss_timeout)) return 1;

  // Print output at initial time
  //PrintOutput(t, y);
  StoreOutput(aT, aY, t, y);

  // Main loop: call CVode, output results, manage events and test for errors.

  unsigned int i_sample = 1;
  realtype tout = (cvode_mode == CV_NORMAL) ? tv[i_sample] : tf / realtype(100.0);

  while(1) {
    // call solver
//printf("xxx0 %.40g %.40g %.40g\\n",t,tout,DT_AFTER(t));
    flag = CVode(cvode_mem, tout, y, &t, cvode_mode);
    if (check_flag(&flag, (char*)"CVode", 1)) break;

    // determine next step
    if (steady_state_flag == 0 && flag == CV_SUCCESS) {
//printf("yyy1 t=%.40f steady_state_flag=%d\\n",t, steady_state_flag);
      //PrintOutput(t, y);
      StoreOutput(aT, aY, t, y);
      if (cvode_mode == CV_NORMAL) {
        i_sample++;
        tout = tv[i_sample];
      } else {
        tout = 0.0;
      }
    } else if (steady_state_flag == 0 && flag == CV_TSTOP_RETURN) {
//printf("yyy2 t=%.40f steady_state_flag=%d\\n",t, steady_state_flag);
      // Check that we stopped at expected time
      if (t != ode_events_minus_dt[ode_event_index]) { // should always be equal
        printf("ERROR: unexpected integrator CV_TSTOP_RETURN\\n");
        break;
      }

      // Output value at (event-dt) only if solver-defined steps
      if (cvode_mode == CV_ONE_STEP) {
        //PrintOutput(t, y);
        StoreOutput(aT, aY, t, y);
      }

      // Advance time to the discontinuity, extend and output y
      t = ode_events[ode_event_index];
      printf("ODE EVENT at t=%.20f\\n",t);

      // Ouput value at event time if solver-defined steps
      // or if on (or near) user-defined sample point
      if ((cvode_mode == CV_ONE_STEP) ||
          (cvode_mode == CV_NORMAL && t <= DT_AFTER(tout) && t>=DT_BEFORE(tout))) {
        //PrintOutput(t, y);
        if (cvode_mode == CV_NORMAL) {
          StoreOutput(aT, aY, tout, y);
          // user-defined sampling, so move to next sample point
          i_sample++;
          tout = tv[i_sample];
        } else {
          StoreOutput(aT, aY, t, y);
          // CV_ONE_STEP mode, need estimate of next output time
          tout = t + (tf - t0)/100.0;
        }
      }

      // Call CVodeReInit to re-initialize the integrator
      flag = CVodeReInit(cvode_mem, t, y);
      if (check_flag(&flag, (char*)"CVodeReInit", 1)) return(1);
      if (cvode_mode == CV_ONE_STEP) tout = t + (tf - t0)/100.0; // estimate output time

      // Update event_flags
      event_flags[ode_event_index] = 1;

      // Update event_times now that event time is known
      event_times[ode_event_index] = t;

      // Setup next event
      ode_event_index++;
      if (setup_next_event(cvode_mem, t, tf,
                           ode_event_index,num_ode_events,
                           ode_events,ode_events_minus_dt,
                           steady_state_flag,ss_timeout)) return 1;

    } else if (steady_state_flag == 1) {

      if (flag == CV_SUCCESS) {
        // steady-state has *NOT* timed out
//printf("yyy3 t=%.40f steady_state_flag=%d\\n",t, steady_state_flag);

        //PrintOutput(t, y);
        StoreOutput(aT, aY, t, y);
        if (cvode_mode == CV_NORMAL) {
          i_sample++;
          tout = tv[i_sample];
        } else {
          tout = 0.0;
        }

        // check for steady-state condition
        int steady_state_reached = 1;
        unsigned int before_last_index = aT.size() - 2;
        realtype delta_t = t - aT[before_last_index];
        unsigned int y_offset = before_last_index * NEQ;
        for (unsigned int i = 0; i < NEQ; i++) {
          // SS stopping condition:
          //   dy/dt * SS_timescale < SS_RelTol * max(abs(y), SS_AbsTol)
          realtype y_last = Ith(y,i);
          realtype delta_y = fabs((y_last - aY[y_offset + i]) * SS_timescale / delta_t);
          realtype abs_y_last = fabs(y_last);
          if (abs_y_last < SS_AbsTol) abs_y_last = SS_AbsTol;
          realtype delta_y_threshold = SS_RelTol * (abs_y_last > SS_AbsTol ? abs_y_last : SS_AbsTol);
          if (delta_y >= delta_y_threshold) {
            steady_state_reached = 0;
          }
        }

        if (steady_state_reached == 1) {
          printf("STEADY-STATE REACHED AT t=%.20f\\n", t);

          // reset steady_state_flag
          steady_state_flag = 0;

          // Update event_flags
          event_flags[ode_event_index] = 1;

          // Update event_times now that event time is known
          event_times[ode_event_index] = t;

          // Stop if steady-state event was the last event
          if (ode_event_index == num_ode_events - 1) break;

          // Call CVodeReInit to re-initialize the integrator
          flag = CVodeReInit(cvode_mem, t, y);
          if (check_flag(&flag, (char*)"CVodeReInit", 1)) return(1);
          if (cvode_mode == CV_ONE_STEP) tout = t + (tf - t0)/100.0; // estimate output time

          // Setup next event
          ode_event_index++;
          if (setup_next_event(cvode_mem, t, tf,
                               ode_event_index,num_ode_events,
                               ode_events,ode_events_minus_dt,
                               steady_state_flag,ss_timeout)) return 1;

        }

      } else if (flag == CV_TSTOP_RETURN) {
        // steady-state *HAS* timed out
//printf("yyy4 t=%.40f steady_state_flag=%d\\n",t, steady_state_flag);

        // Check that we stopped at expected time
        if (t != ss_timeout) { // should always be equal
          printf("ERROR: unexpected integrator CV_TSTOP_RETURN\\n");
          break;
        }

        // Advance time to the actual timeout value
        t = ss_timeout = DT_AFTER(ss_timeout);

        printf("STEADY-STATE TIMED-OUT AT t=%.20f\\n", t);
//printf("yyy5 t=%.40f steady_state_flag=%d\\n",t, steady_state_flag);

        // Ouput value at event time if solver-defined steps or if on user-defined sample point
        if (cvode_mode == CV_ONE_STEP || (cvode_mode == CV_NORMAL && t == tout)) {
          //PrintOutput(t, y);
          StoreOutput(aT, aY, t, y);
          if (cvode_mode == CV_NORMAL) {
            // user-defined sampling, so move to next sample point
            i_sample++;
            tout = tv[i_sample];
          } else {
            // CV_ONE_STEP mode, need estimate of next output time
            tout = t + (tf - t0)/100.0;
          }
        }

        // reset steady_state_flag
        steady_state_flag = 0;

        // Update event_flags
        event_flags[ode_event_index] = 1;

        // Update event_times now that event time is known
        event_times[ode_event_index] = t;

        // Stop if we timed out on tf
        if (t >= tf) break;

        // We timed-out, so advance to event corresponding to ss_timeout if any
        while(ode_event_index < num_ode_events &&
              ode_events[ode_event_index] != ss_timeout) {
          ode_event_index++;
        }

        // Call CVodeReInit to re-initialize the integrator
        flag = CVodeReInit(cvode_mem, t, y);
        if (check_flag(&flag, (char*)"CVodeReInit", 1)) return(1);
        if (cvode_mode == CV_ONE_STEP) tout = t + (tf - t0)/100.0; // estimate output time

        // Setup next event
        ode_event_index++;
        if (setup_next_event(cvode_mem, t, tf,
                             ode_event_index,num_ode_events,
                             ode_events,ode_events_minus_dt,
                             steady_state_flag,ss_timeout)) return 1;

      }

    } else {
      printf("ERROR: unexpected condition\\n");
      break;
    }

    // Stopping condition
    if (t >= tf) break;
  }

  // Print some final statistics
  PrintFinalStats(cvode_mem);

  // Free y and abstol vectors
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(abstol);

  // Free integrator memory
  CVodeFree(&cvode_mem);

  return(0);
}

// Private helper functions
static int setup_next_event(void *cvode_mem, realtype t, realtype tf,
                            unsigned int &ode_event_index,
                            unsigned int num_ode_events,
                            std::vector<realtype> &ode_events,
                            std::vector<realtype> &ode_events_minus_dt,
                            int &steady_state_flag, realtype &ss_timeout
                           ) {
  int flag;
  ss_timeout = 0;

  if ((ode_event_index < num_ode_events) && (ode_events[ode_event_index] > 0)) {
    // positive, absolute event time
    realtype next_event = ode_events_minus_dt[ode_event_index];
    flag = CVodeSetStopTime(cvode_mem, next_event);
    if (check_flag(&flag, (char*)"CVodeSetStopTime", 1)) return(1);
  } else if ((ode_event_index < num_ode_events) && (ode_events[ode_event_index] < 0)) {
    // negative event time, so the next event time is relative to current time
    ode_events[ode_event_index] = t - ode_events[ode_event_index]; // update event time now that it is known
    ode_events_minus_dt[ode_event_index] = DT_BEFORE(ode_events[ode_event_index]);
    realtype next_event = ode_events_minus_dt[ode_event_index];
    flag = CVodeSetStopTime(cvode_mem, next_event);
    if (check_flag(&flag, (char*)"CVodeSetStopTime", 1)) return(1);
  } else if ((ode_event_index < num_ode_events) && (ode_events[ode_event_index] == 0)) {
    // next event is steady-state
    steady_state_flag = 1;

    // compute time-out
    for (unsigned int i = ode_event_index + 1; i < num_ode_events; i++) {
      if (ode_events[i] > 0) {
        ss_timeout = ode_events[i];
        break;
      }
    }
    if (ss_timeout == 0) ss_timeout = tf;
    ss_timeout = DT_BEFORE(ss_timeout);
    // don't let integrator go past steady-state timeout
    flag = CVodeSetStopTime(cvode_mem, ss_timeout);
    if (check_flag(&flag, (char*)"CVodeSetStopTime", 1)) return(1);
  } else {
    // no events left
    flag = CVodeSetStopTime(cvode_mem, 2.0*tf);
    if (check_flag(&flag, (char*)"CVodeSetStopTime", 1)) return(1);
  }

  return 0;
}



static void StoreOutput(std::vector<realtype> &aT,
                        std::vector<realtype> &aY,
                        realtype t, N_Vector y) {
  // store t
  aT.push_back(t);

  // store y(t)
  realtype *y_ptr = &(Ith(y,0));
  for(int i=0; i < NEQ; i++) {
    aY.push_back(y_ptr[i]);
  }
}

CVODE

    # various auxiliary functions
    $cvode_file_contents .= <<'CVODE';
static void PrintOutput(realtype t, N_Vector y) {
  printf("At t = %0.20e      y[] = ",t);
  for (unsigned int i=0; i < NEQ; i++) {
    printf("%.20e  ", Ith(y,i));
  }
  printf("\n");
  return;
}

int cvode_sim_print_output(std::vector<realtype> &aT,
                           std::vector<realtype> &aY) {
  unsigned int num_samples = aT.size();
  for(unsigned int i = 0; i < num_samples; i++) {
    printf("At t = %0.20e      y[] = ",aT[i]);
    for (unsigned int j=0; j < NEQ; j++) {
      printf("%.20e  ", aY[i*NEQ + j]);
    }
    printf("\n");
  }
}

// Get and print some final statistics

static void PrintFinalStats(void *cvode_mem) {
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, (char*)"CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, (char*)"CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, (char*)"CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, (char*)"CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, (char*)"CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, (char*)"CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, (char*)"CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, (char*)"CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, (char*)"CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

//
// Check function return value...
//   opt == 0 means SUNDIALS function allocates memory so check if
//            returned NULL pointer
//   opt == 1 means SUNDIALS function returns a flag so check if
//            flag >= 0
//   opt == 2 means function allocates memory so check if returned
//            NULL pointer 
//

static int check_flag(void *flagvalue, char *funcname, int opt) {
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}

CVODE

    #********************************************************#
    # generate dydt function                                 #
    #********************************************************#
    # ODE function

    $ode_func_contents{header} .= <<'CVODE';
static int f_dydt(realtype t, N_Vector y, N_Vector ydot, void * user_data) {
//    printf("f_dydt call at t=%.20f\n",t);

    void** user_data_array = (void**) user_data;

    realtype *ode_rate_constants       = (realtype *)                user_data_array[0];
    std::vector<int> &event_flags      = *((std::vector<int> *)      user_data_array[1]);
    std::vector<realtype> &event_times = *((std::vector<realtype> *) user_data_array[2]);

CVODE

#    $ode_func_contents{header} .= "    global ode_tot_cputime\n";
#    $ode_func_contents{header} .= "    global ode_num_calls\n";
#    $ode_func_contents{header} .= "    ode_start_time = time.clock()\n";
#    $ode_func_contents{header} .= "    ode_num_calls = ode_num_calls + 1\n";
    $ode_func_contents{all} .= $ode_func_contents{header};

#    # Clock tick
#    # (incomplete implementation)
#    my $tick_line = ($tk == -1) ? "" : "print 'ode: sim time is t =',t\n\n";
#    $ode_func_contents{tick} .= "    $tick_line";
#    $ode_func_contents{all} .= $ode_func_contents{tick};
	
    # map state vector to free nodes
    $ode_func_contents{node_map} .= "    // state vector to node mapping\n" if (@constant_rate_params);
    for (my $j = 0; $j < @free_nodes; $j++) {
	my $node_ref = $free_nodes[$j];
	my $node_name = $node_ref->get_name();
	$ode_func_contents{node_map} .= "    realtype $node_name = Ith(y,$j);\n";
    }	
    $ode_func_contents{node_map} .= "\n";
    $ode_func_contents{all} .= $ode_func_contents{node_map};

    # ordinary rate constants (e.g. f1=1) and constant rate expressions (e.g. f2=2*f1)
    $ode_func_contents{consts} .= "    // rate constants and constant expressions\n" if (@constant_rate_params);
    for (my $i = 0; $i < @ode_rate_constants; $i++) {
	my $variable_ref = $ode_rate_constants[$i];
	my $variable_name = $variable_ref->get_name();
	my $variable_value = $variable_ref->get_value();
	$ode_func_contents{consts} .= "    realtype $variable_name = ode_rate_constants[$i];\n";
    }
    $ode_func_contents{consts} .= "\n";
    $ode_func_contents{all} .= $ode_func_contents{consts};

    # moiety totals, e.g. E_moiety = 123 (mol/L)
    $ode_func_contents{moiety} .= "    // moiety totals\n" if (@moiety_totals);
    foreach my $variable_ref (@moiety_totals) {
	my $variable_index = $variable_ref->get_index();
	my $variable_name = $variable_ref->get_name();
	my $variable_value = $variable_ref->get_value();
	$ode_func_contents{moiety} .= "    realtype $variable_name = $variable_value;\n";
    }
    $ode_func_contents{moiety} .= "\n";
    $ode_func_contents{all} .= $ode_func_contents{moiety};

    # contrained node expressions, e.g. E = C - E_moiety
    $ode_func_contents{cnodes} .= "    // dependent species\n" if (@constrained_node_expressions);
    foreach my $variable_ref (@constrained_node_expressions) {
	my $variable_index = $variable_ref->get_index();
	my $variable_name = $variable_ref->get_name();
	my $variable_value = $variable_ref->get_value();
	$ode_func_contents{cnodes} .= "realtype $variable_name = $variable_value;\n";
    }
    $ode_func_contents{cnodes} .= "\n";
    $ode_func_contents{all} .= $ode_func_contents{cnodes};

    # rate expressions
    $ode_func_contents{dynrates} .= "    // dynamic rate expressions\n" if (@dynamic_rate_expressions);
    foreach my $variable_ref (@dynamic_rate_expressions) {
	my $variable_index = $variable_ref->get_index();
	my $variable_name = $variable_ref->get_name();
	my $variable_value = $variable_ref->get_value();
	$ode_func_contents{dynrates} .= "    realtype $variable_name = $variable_value;\n";
    }
    $ode_func_contents{dynrates} .= "\n";
    $ode_func_contents{all} .= $ode_func_contents{dynrates};

    # print out differential equations for the free nodes
    $ode_func_contents{dydt} .= "    // differential equations for independent species\n";

    my @ode_rhs = ();
    for (my $j = 0; $j < @free_nodes; $j++) {
	my $node_ref = $free_nodes[$j];
	my $node_name = $node_ref->get_name();

	my $ode_rhs = "";

	my @create_reactions = $node_ref->get_create_reactions();
	my @destroy_reactions = $node_ref->get_destroy_reactions();

	# print positive terms
	foreach my $reaction_ref (@create_reactions) {
	    my $velocity = $reaction_ref->get_velocity();
	    $ode_rhs .= "+ $velocity ";
	}
	# print negative terms
	foreach my $reaction_ref (@destroy_reactions) {
	    my $velocity = $reaction_ref->get_velocity();
	    $ode_rhs .= "- $velocity ";
	}

	push @ode_rhs, $ode_rhs;

	if ($factor_flag) {
	    $ode_rhs = Expression->new({value=>$ode_rhs})->factor_expression();
 	}

	if ($ode_rhs ne "") {
	    $ode_func_contents{dydt} .= "    Ith(ydot,$j) = $ode_rhs;\n";
	} else {
	    $ode_func_contents{dydt} .= "    Ith(ydot,$j) = 0;\n";
	}
    }
    $ode_func_contents{all} .= $ode_func_contents{dydt};

#    $ode_func_contents{footer} .= "    ode_end_time = time.clock()\n";
#    $ode_func_contents{footer} .= "    ode_tot_cputime = ode_tot_cputime + (ode_end_time - ode_start_time);\n";
    $ode_func_contents{footer} .= "\n    return (0);\n}\n";
    $ode_func_contents{all} .= $ode_func_contents{footer};

    # substitute indexing operator for event_flags() and event_times()
    $ode_func_contents{all} =~ s/event_flags\((.*?)\)/event_flags\[$1 - 1\]/g;
    $ode_func_contents{all} =~ s/event_times\((.*?)\)/event_times\[$1 - 1\]/g;

    $cvode_file_contents .= $ode_func_contents{all};

    #********************************************************#
    # generate jacobian function                             #
    #********************************************************#
    $jac_func_contents .= <<CVODE;

static int f_jac($int_N_str, realtype t,
                 N_Vector y, N_Vector fy, DlsMat J, void * user_data,
                 N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {

    void** user_data_array = (void**) user_data;

    realtype *ode_rate_constants       = (realtype *)                user_data_array[0];
    std::vector<int> &event_flags      = *((std::vector<int> *)      user_data_array[1]);
    std::vector<realtype> &event_times = *((std::vector<realtype> *) user_data_array[2]);

CVODE

#    $jac_func_contents .= "    global jac_tot_cputime\n";
#    $jac_func_contents .= "    global jac_num_calls\n";
#    $jac_func_contents .= "    jac_start_time = time.clock()\n";
#    $jac_func_contents .= "    jac_num_calls = jac_num_calls + 1\n";

#    # Clock tick
#    if ($tk != -1) {
#	my $tick_line = ($tk == -1) ? "" : "print 'jac: sim time is t =',t\n\n";
#	$jac_func_contents .= "    $tick_line";
#    }

    # other sections are the same as the ode file
    $jac_func_contents .= $ode_func_contents{node_map};
    $jac_func_contents .= $ode_func_contents{consts};
    $jac_func_contents .= $ode_func_contents{moiety};
    $jac_func_contents .= $ode_func_contents{cnodes};
    $jac_func_contents .= $ode_func_contents{dynrates};

    # now generate equations for jacobian
    $jac_func_contents .= <<CVODE;
    // jacobian equations for independent species
    for (unsigned int i=0; i < $num_free_nodes; i++) {
        for (unsigned j=0; j < $num_free_nodes; j++) {
            IJth(J,i,j) = 0.0;
        }
    }
CVODE
    for (my $j = 0; $j < @free_nodes; $j++) {
	my $ode_rhs_ex_ref = Expression->new({value=>$ode_rhs[$j]});
	for (my $k = 0; $k < @free_nodes; $k++) {
	    my $dvar = $free_nodes[$k]->get_name();
	    my $jac_rhs = $ode_rhs_ex_ref->differentiate_expression($dvar);
	    if ($jac_rhs ne "0") {
		if ($factor_flag) {
		    $jac_rhs = Expression->new({value=>$jac_rhs})->factor_expression();
		}
		$jac_func_contents .= "    IJth(J,$j,$k) = $jac_rhs;\n";
	    }
	}
    }


#    # timings
#    $jac_func_contents .= "    jac_end_time = time.clock()\n";
#    $jac_func_contents .= "    jac_tot_cputime = jac_tot_cputime + (jac_end_time - jac_start_time);\n";

    $jac_func_contents .= "\n    return (0);\n}\n";

    # substitute indexing operator for event_flags() and event_times()
    $jac_func_contents =~ s/event_flags\((.*?)\)/event_flags\[$1 - 1\]/g;
    $jac_func_contents =~ s/event_times\((.*?)\)/event_times\[$1 - 1\]/g;

    $cvode_file_contents .= $jac_func_contents;

    # POST-PROCESSING


    # substitute exponentiation operator "^" -> "pow()";
    $driver_file_contents = hat_to_pow(\$driver_file_contents);
    $cvode_file_contents = hat_to_pow(\$cvode_file_contents);

    return ($driver_file_contents,
	    $cvode_header_file_contents,
	    $cvode_file_contents,
	    $octwrap_file_contents);
}

1;  # don't remove -- req'd for module to return true

