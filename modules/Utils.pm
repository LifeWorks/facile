###############################################################################
## File:     Utils.pm
## Synopsys: Miscellaneous small routines.
###############################################################################
## Detailed Description:
## ---------------------
##
###############################################################################

package Utils;

use strict;
use diagnostics;          # equivalent to -w command-line switch
use warnings;

use vars qw(@ISA @EXPORT);

require Exporter;
@ISA = qw(Exporter);

@EXPORT = qw(
	     is_numeric
	     hat_to_pow
	    );

#--------------------------------------------------- #
# Useful routine to check if a string is a number
#--------------------------------------------------- #
sub is_numeric {
    my $value = shift;
    # the following regexp is from Perl Cookbook (R2.1) to determine if string is a number
    if ($value =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/) {
	return 1;
    } else {
	return 0;
    }
}

sub hat_to_pow {
  my $ref = shift;
  my @in_lines = split(/^/,$$ref);
  my @out_lines = ();

  foreach my $line (@in_lines) {
    while ($line =~ /\^/) {
      #print "$line\n";
      my $hat_ptr = index($line,'^');

      my $r_paren_count = 0;
      my $r_paren_flag = 0;
      my $r_ptr = $hat_ptr + 1;
      while (substr($line,$r_ptr,1) =~ /\s/) {$r_ptr++};
      my $r_char = substr($line,$r_ptr,1);
      while (($r_ptr < length($line)) && 
	     (($r_char =~ /[_a-zA-Z0-9(]/ && $r_paren_flag == 0) || $r_paren_count != 0)) {
	$r_paren_count++ if ($r_char eq '(');
	$r_paren_count-- if ($r_char eq ')');
	$r_paren_flag = 1 if $r_paren_count;
	$r_ptr++;
	$r_char = substr($line,$r_ptr,1);
      }
      my $exponent = substr($line,$hat_ptr+1,$r_ptr-$hat_ptr-1);

      my $l_paren_count = 0;
      my $l_paren_flag = 0;
      my $l_ptr = $hat_ptr - 1;
      while (substr($line,$l_ptr,1) =~ /\s/) {$l_ptr--};
      my $l_char = substr($line,$l_ptr,1);
      while (($l_ptr < length($line)) && 
	     (($l_char =~ /[_a-zA-Z0-9)]/ && $l_paren_flag == 0) || $l_paren_count != 0)) {
	$l_paren_count++ if ($l_char eq ')');
	$l_paren_count-- if ($l_char eq '(');
	$l_paren_flag = 1 if $l_paren_count;
	$l_ptr--;
	$l_char = substr($line,$l_ptr,1);
      }
      my $base = substr($line,$l_ptr+1,$hat_ptr-$l_ptr-1);

      substr($line,$l_ptr+1,$r_ptr-$l_ptr-1) = "pow($base,$exponent)";
      
      #print "$line\n\n";
    }
    push @out_lines, $line;
  }

  return join("",@out_lines);
}
