######################################################################################
# File:     Expression.pm
# Synopsys: Create, simplify, differentiate and print symbolic expressions.
#
# Copyright (C) 2005-2011 Ollivier JF
#
# This program comes with ABSOLUTELY NO WARRANTY.
# This is free software, and you are welcome to redistribute it
# under certain conditions. See file LICENSE.TXT for details.
#
######################################################################################
# Detailed Description:
# ---------------------
#
#
######################################################################################

use strict;
use diagnostics;		# equivalent to -w command-line switch
use warnings;

package Expression;
use Class::Std;
use base qw();
{
    use Carp;
    use Utils;

    #######################################################################################
    # CLASS ATTRIBUTES
    #######################################################################################

    #######################################################################################
    # ATTRIBUTES
    #######################################################################################
    my %terms_ref_of :ATTR(get => 'terms_ref');
    my %vars_ref_of :ATTR(get => 'vars_ref');

    my %value_of :ATTR(get => 'value');

    #######################################################################################
    # FUNCTIONS
    #######################################################################################

    #######################################################################################
    # CLASS METHODS
    #######################################################################################

    #######################################################################################
    # INSTANCE METHODS
    #######################################################################################
    #--------------------------------------------------------------------------------------
    # Function: START
    # Synopsys: 
    #--------------------------------------------------------------------------------------
    sub BUILD {
        my ($self, $obj_ID, $arg_ref) = @_;

	# check initializers
	# ...

	if (defined $arg_ref->{value}) {
	    $value_of{$obj_ID} = $arg_ref->{value};
	    $self->build_tables();
	} else {
	    if (defined $arg_ref->{terms_ref}) {
		$terms_ref_of{$obj_ID} = $arg_ref->{terms_ref};
	    }
	    if (defined $arg_ref->{vars_ref}) {
		$vars_ref_of{$obj_ID} = $arg_ref->{vars_ref};
	    }
	}
    }

    #--------------------------------------------------------------------------------------
    # Function: build_tables
    # Synopsys: Build table of vars and terms from string form of expression.
    #--------------------------------------------------------------------------------------
    sub build_tables {
	my $self = shift;

	my $expression = $value_of{ident $self};;

	$expression =~ s/\s//g;	# remove whitespace
	$expression = "+$expression" if ($expression !~ /^[+-]/); # insert leading +
	$expression =~ s/(?<!^)([+-])/ $1/g; # insert a space before each term for easier split

	my @terms = split(/ /, $expression);
	@terms = map {[split(/(?<=[+-])|\*/,$_)]} @terms;

	my $terms_ref = \@terms;
	my $vars_ref = {};

	for (my $i=0; $i < @$terms_ref; $i++) {
	    my $names_ref = $terms_ref->[$i];
	    for (my $j=1; $j < @$names_ref; $j++) {
		my $var_name = $names_ref->[$j];
		push @{$vars_ref->{$var_name}}, $i;
	    }
	}

	$vars_ref_of{ident $self} = $vars_ref;
	$terms_ref_of{ident $self} = $terms_ref;
    }

    #--------------------------------------------------------------------------------------
    # Function: sprint
    # Synopsys: Routine to print out terms to string
    #--------------------------------------------------------------------------------------
    sub sprint {
	my $self = shift;

	my $terms_ref = $terms_ref_of{ident $self};

	my $rval = "";
	foreach my $term_vars_ref (@$terms_ref) {
	    if (@$term_vars_ref > 1) {
		# print out sign
		$rval .= " $term_vars_ref->[0] ";
		# print out var products
		my @term_vars = @{$term_vars_ref}[1..$#{$term_vars_ref}];
		$rval .= join("*",@term_vars);
	    } elsif (@$term_vars_ref == 1) {
		# term had just one var, which got factored out, leaving only the sign
		$rval .= " $term_vars_ref->[0] 1 ";
	    }
	}

	return $rval;
    }

    #--------------------------------------------------------------------------------------
    # Function: factor_expression
    # Synopsys: Returns factored expression in string form.
    #--------------------------------------------------------------------------------------
    sub factor_expression {
	my $self = shift;

	my $vars_ref = $vars_ref_of{ident $self};
	my $terms_ref = $terms_ref_of{ident $self};

	my $rval = "";

	# add up constant terms
	my @constant_terms = grep {@{$_} == 1} @$terms_ref;
	my $constant = 0;
	map {$constant += ($_->[0] eq '+') ? +1 : -1} @constant_terms;
	if (@constant_terms) {
	    if (@$terms_ref == @constant_terms) {
		# there are only constant terms
		$rval .= $constant; # result of addition can be zero
	    } else {
		# there are other terms, append only if non-zero
		$rval .= $constant if $constant != 0; # result of addition can be zero
	    }
	}

	my @var_list = keys %$vars_ref;

	# sort vars according to first term they appear in
	# (this help preserve proper output order for vars that appear in equal number of terms)
	@var_list = sort {$vars_ref->{$a}[0] <=> $vars_ref->{$b}[0]} @var_list;

	while (1) {
	    # break out if no vars to process
	    last if !@var_list;

	    # sort vars according to how many terms they appear in
	    @var_list = sort {@{$vars_ref->{$b}} <=> @{$vars_ref->{$a}}} @var_list;

	    my $mcvar = $var_list[0];
	    my @term_no_list = @{$vars_ref->{$mcvar}};
	    last if (@term_no_list == 0); # done?

	    # remove one or more consecutive duplicates in term_no_list
	    # (can happen when same var appears twice in term e.g. b*a*a*a)
	    @term_no_list = map {$term_no_list[$_]} (grep {!defined $term_no_list[$_ + 1] ||
							   $term_no_list[$_] != $term_no_list[$_ + 1]} (0..$#term_no_list));

	    if (@term_no_list > 1) { # var appears in more than one term?
		my %sub_vars = ();
		my @sub_terms = ();
		for (my $i=0; $i < @term_no_list; $i++) {
		    my $term_no = $term_no_list[$i];
		    push @{$sub_terms[$i]}, $terms_ref->[$term_no][0]; # sign
		    my @term_vars = @{$terms_ref->[$term_no]}[1..$#{$terms_ref->[$term_no]}];
		    my $already_found = 0;
		    foreach my $var_name (@term_vars) {
			if (($var_name ne $mcvar) || $already_found) {
			    push @{$sub_terms[$i]}, $var_name;
			    push @{$sub_vars{$var_name}}, $i;
			} else {
			    $already_found = 1;
			}
			# delete term number from var's term list
			@{$vars_ref->{$var_name}} = grep {$_ != $term_no} @{$vars_ref->{$var_name}};
		    }
		}
		# option 1: recurse
		my $sub_ref = Expression->new({vars_ref=>\%sub_vars, terms_ref=>\@sub_terms});
		my $temp = $sub_ref->factor_expression();
		# option 2: no-recurse
		#my $temp = $sub_ref->sprint();
		if (!is_numeric($temp) || $temp != 0) {
		    $rval .= " + $mcvar";
		    $rval .= "*($temp)";
		}
	    } else {
		my $term_no = $term_no_list[0];
		my @term_vars = @{$terms_ref->[$term_no]}[1..$#{$terms_ref->[$term_no]}];
		foreach my $var_name (@term_vars) {
		    # delete term numbers from var's term list
		    @{$vars_ref->{$var_name}} = grep {$_ != $term_no} @{$vars_ref->{$var_name}};
		}
		my $ref = 
		$rval .= Expression->new({terms_ref => [$terms_ref->[$term_no]]})->sprint();
	    }
	}

	$rval = 0 if $rval eq "";
	return $rval;
    }

    #--------------------------------------------------------------------------------------
    # Function: differentiate_expression
    # Synopsys: Returns differentiated expression in string form.
    #--------------------------------------------------------------------------------------
    sub differentiate_expression {
	my $self = shift;
	my $dvar = shift;

	my $vars_ref = $vars_ref_of{ident $self};
	my $terms_ref = $terms_ref_of{ident $self};

	my $rval = "";

	if (defined $vars_ref->{$dvar} && @{$vars_ref->{$dvar}} > 0) {
	    my @rvals = ();
	    for (my $i=0; $i < @{$vars_ref->{$dvar}}; $i++) {
		my $term_no = ${$vars_ref->{$dvar}}[$i];

		# avoid repeats
		next if ($i != 0) && ($term_no == ${$vars_ref->{$dvar}}[$i-1]);

		my $num_dvars = 0;
		my @vars = ();
		foreach my $var (@{$terms_ref->[$term_no]}) {
		    push @vars, $var if ($var ne $dvar || $num_dvars >= 1);
		    $num_dvars++ if $var eq $dvar;
		}

		if (@vars > 1) {
		    push @rvals, "$vars[0] ".($num_dvars != 1 ? "$num_dvars*" : "").join("*",@vars[1..$#vars]);
		} else {
		    push @rvals, "$vars[0] ".$num_dvars;
		}
	    }
	    $rval = join " ", @rvals;
	} else {
	    $rval = "0";
	}
	return $rval;
    }
}


sub run_testcases {
    print "Running Factor testcase...\n";

    my $ex_ref;

    $ex_ref = Expression->new({value => "a*b + a*b"});
    print $ex_ref->get_value()." = ".$ex_ref->factor_expression() ."\n";

    $ex_ref = Expression->new({value => "+a*b + a*c + a"});
    print $ex_ref->get_value()." = ".$ex_ref->factor_expression() ."\n";

    $ex_ref = Expression->new({value => "a*b + a*c + a*a*b"});
    print $ex_ref->get_value()." = ".$ex_ref->factor_expression() ."\n";

    $ex_ref = Expression->new({value => "a*b + a*c + a + 5"});
    print $ex_ref->get_value()." = ".$ex_ref->factor_expression() ."\n";

    $ex_ref = Expression->new({value => "a*b + a*c + a*a*a*a + a*a*b*b"});
    print $ex_ref->get_value()." = ".$ex_ref->factor_expression() ."\n";

    $ex_ref = Expression->new({value => "a*b + a*c + a*a*a*a + a*b*a*b"});
    print $ex_ref->get_value()." = ".$ex_ref->factor_expression() ."\n";

    $ex_ref = Expression->new({value => "5+3+6+5+2+5"});
    print $ex_ref->get_value()." = ".$ex_ref->factor_expression() ."\n";

    $ex_ref = Expression->new({value => "+ D - D"});
    print $ex_ref->get_value()." = ".$ex_ref->factor_expression() ."\n";

    $ex_ref = Expression->new({value => "+ k*D - k*D"});
    print $ex_ref->get_value()." = ".$ex_ref->factor_expression() ."\n";

    $ex_ref = Expression->new({value => "f + D - D"});
    print $ex_ref->get_value()." = ".$ex_ref->factor_expression() ."\n";

    $ex_ref = Expression->new({value => "f + k*D - k*D"});
    print $ex_ref->get_value()." = ".$ex_ref->factor_expression() ."\n";

    $ex_ref = Expression->new({value => "f+k*D+f*D- k*D - f*D"});
    print $ex_ref->get_value()." = ".$ex_ref->factor_expression() ."\n";

    $ex_ref = Expression->new({value => "f + k*D + k*D*D + 2*k*D*A + 2*D*D"});
    print "d/dD ".$ex_ref->get_value()." = ".$ex_ref->differentiate_expression("D") ."\n";

    $ex_ref = Expression->new({value => "f+k"});
    print "d/dD ".$ex_ref->get_value()." = ".$ex_ref->differentiate_expression("D") ."\n";

    $ex_ref = Expression->new({value => "5"});
    print "d/dD ".$ex_ref->get_value()." = ".$ex_ref->differentiate_expression("D") ."\n";

}


# Package BEGIN must return true value
return 1;

