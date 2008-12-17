#!/usr/bin/perl
# rotstev.pl 
#
# Generates rotation matrices to rotate Stevens' Operator equivalents by
# phi about y-axis and theta about z-axis. Based on the method of Buckmaster
# (Physica Status Solidi A vol 13, pp 9, 1972) and Rudowicz (J. Phys: C: 
# Solid State Physics, vol 18, pp 1415, 1985).
#
# NB: phi is denoted by f and theta by t in this program
#
# by Duc Le 2006 - duc.le@ucl.ac.uk

use Math::Algebra::Symbols trig=>1;     # exports trig function to my namespace
use Getopt::Long;

# The following subroutines were shamelessly stolen from the Perl Cookbook
# by O'Reily, and used for the symbolic matrix multiplications.

sub mmult {
  my ($m1,$m2) = @_;
  my ($m1rows,$m1cols) = matdim($m1);
  my ($m2rows,$m2cols) = matdim($m2);

  unless ($m1cols == $m2rows) {       # raise exception
    die "IndexError: matrices don't match: $m1cols != $m2rows";
  }

  my $result = [];
  my ($i, $j, $k, $m1el, $m2el, $prod_el);

  for $i (range($m1rows)) {
    for $j (range($m2cols)) {
      for $k (range($m1cols)) {
        $m1el = $m1->[$i][$k];
        $m2el = $m2->[$k][$j];
        $prod_el += $m1el * $m2el;
      }
      $result->[$i][$j] = $prod_el;
      $prod_el = 0;
    }
  }

  return $result;
}

sub range { 0 .. ($_[0] - 1) }

sub veclen {
    my $ary_ref = $_[0];
    my $type = ref $ary_ref;
    if ($type ne "ARRAY") { die "$type is bad array ref for $ary_ref" }
    return scalar(@$ary_ref);
}

sub matdim {
    my $matrix = $_[0];
    my $rows = veclen($matrix);
    my $cols = veclen($matrix->[0]);
    return ($rows, $cols);
}

# End of stolen code!

sub usage() {
  print STDERR << "EOF";

    $0:
    Calculates the relationship  between a set of crystal field  parameters
    which has  been rotated by phi about  the y-axis and theta about the z-
    axis and its unrotated counter-part.

    The  calculations  are  done  by  means  of  matrix  multiplication and 
    symbolic algebra  using the Math::Algebra::Symbols package and based on
    the method of  Buckmaster (phys. stat. sol. a, vol 13,  pp 9, 1972) and
    Rudowicz (J. Phys: Solid State Phys., vol 18, pp 1415, 1985).   

    usage: $0 [-h] [--help] 
              [-r] [-radians]
	      [-p phi_angle] [--phi phi_angle]
	      [-t theta_angle] [--theta theta_angle]
	      [-o output_file] [--output output_file]

     -h             : this (help) message
     -r             : specifies that phi and theta are given in radians
     -p phi_angle   : angle to rotate about y-axis (=0 by default)
     -t theta_angel : angle to rotate about z-axis (=90 degrees by default)  
     -o out_file    : output latex file name

    if -r is omitted, the program assumes phi and theta are in degrees.
    
    if either -p or -t are omitted, the program will use the default values 
    
    if -o is omitted the  program will output the  transformation relations
    and rotation matrix to screen.  Otherwise it will generate a Latex file
    of the transformations and matrix.

    (C) 2006 duc.le\@ucl.ac.uk
EOF
  exit;
}

# see http://aplawrence.com/Unix/perlgetops.html for details of GetOptions

GetOptions("help"=>\$helpflag,
	   "radians"=>\$radflag,
	   "phi=f"=>\$phi,
	   "theta=f"=>\$theta,
	   "output=s"=>\$output,
           "z_internal"=>\$z_int);
	   
usage() if $helpflag;

my ($S, $f, $t, $i, $pi) = symbols(qw(six f t i pi));
my ($l, $m, $n);

if (!$phi) { $phi = 0; }
if (!$theta) { $theta = 90; }

if (!$radflag) {
  $phi_val = $phi * $pi / 180;
  $theta_val = $theta * $pi / 180;
}
else {
  $phi_val = $phi;
  $theta_val = $theta;
}

# Transformation matrix between Stevens and Buckmaster operator equivalents
my $A2 = [
	[1,	0,	0,	   0,	  -1],
	[0,	1/2,	0,	   1/2,	   0],
	[0,	0,	sqrt($S),  0,	   0],
	[0,	1/2,	0,	  -1/2,	   0],
	[1,	0,	0,	   0,	   1]
     ];

my $invA2 = [
	[       1/2,         0,         0,         0,       1/2],
	[         0,         1,         0,         1,         0],
	[         0,         0, sqrt(1/$S),        0,         0],
	[         0,         1,         0,        -1,         0],
	[      -1/2,         0,         0,         0,       1/2]
      ];

# --------------------------------------- O2 ------------------------------------- #

my $D20p2 = (sqrt($S)/4)*sin($t)**2;
my $D20p1 = (sqrt($S)/2)*sin($t)*cos($t);
my $D20p0 = (1/2)*(3*cos($t)**2 - 1);

my $D2m1m2 = exp(-$i*$f) * -(1/2)*sin($t)*(1+cos($t));
my $D2m1p2 = exp(-$i*$f) *  (1/2)*sin($t)*(1-cos($t));
my $D2m1m1 = exp(-$i*$f) *  (1/2)*(1+cos($t))*(2*cos($t)-1);
my $D2m1p1 = exp(-$i*$f) *  (1/2)*(1-cos($t))*(2*cos($t)+1);
my $D2m1p0 = exp(-$i*$f) *  sqrt($S)/2*sin($t)*cos($t);

my $D2p1m2 = exp( $i*$f) * -(1/2)*sin($t)*(1-cos($t));
my $D2p1p2 = exp( $i*$f) *  (1/2)*sin($t)*(1+cos($t));
my $D2p1m1 = exp( $i*$f) *  (1/2)*(1-cos($t))*(2*cos($t)+1);
my $D2p1p1 = exp( $i*$f) *  (1/2)*(1+cos($t))*(2*cos($t)-1);
my $D2p1p0 = exp( $i*$f) * -sqrt($S)/2*sin($t)*cos($t);

my $D2m2m2 = exp(-2*$i*$f) *  (1/4)*(1+cos($t))**2;
my $D2m2p2 = exp(-2*$i*$f) *  (1/4)*(1-cos($t))**2;
my $D2m2m1 = exp(-2*$i*$f) *  (1/2)*sin($t)*(1+cos($t));
my $D2m2p1 = exp(-2*$i*$f) *  (1/2)*sin($t)*(1-cos($t));
my $D2m2p0 = exp(-2*$i*$f) *  sqrt($S)/4*sin($t)**2;

my $D2p2m2 = exp( 2*$i*$f) *  (1/4)*(1-cos($t))**2;
my $D2p2p2 = exp( 2*$i*$f) *  (1/4)*(1+cos($t))**2;
my $D2p2m1 = exp( 2*$i*$f) * -(1/2)*sin($t)*(1-cos($t));
my $D2p2p1 = exp( 2*$i*$f) * -(1/2)*sin($t)*(1+cos($t));
my $D2p2p0 = exp( 2*$i*$f) *  sqrt($S)/4*sin($t)**2;

my $D2 = [
	   [$D2m2m2, $D2m2m1, $D2m2p0, $D2m2p1, $D2m2p2],
           [$D2m1m2, $D2m1m1, $D2m1p0, $D2m1p1, $D2m1p2],
           [$D20p2, -$D20p1,  $D20p0,  $D20p1,  $D20p2 ],
           [$D2p1m2, $D2p1m1, $D2p1p0, $D2p1p1, $D2p1p2],
           [$D2p2m2, $D2p2m1, $D2p2p0, $D2p2p1, $D2p2p2]
         ];

my ($O2m2, $O2m1, $O2p0, $O2p1, $O2p2) = symbols(qw(Otmt Otmo Otpz Otpo Otpt));

# ------------------------------------ Finished !! ------------------------------- #

# The operator in the rotated frame. Rx(y) is equivalent to {Oxy} in Rudowicz notation
# and Ox_y, Oxmy is equivalent to [Oxy], [OxyM] in Rudowicz

my $B2 = [ [$O2m2],[$O2m1],[$O2p0],[$O2p1],[$O2p2] ];

my $tmp1 = mmult($invA2,$B2);
my $tmp2 = mmult($D2, $tmp1);
my $R2 = mmult($A2, $tmp2);

my $S2 = [];
for $l (0 .. 4) { 
  $S2->[$l] = $R2->[$l][0]->sub(f=>$phi_val, t=>$theta_val);
  $S2->[$l] = $S2->[$l]->sub(six=>6); 
  $M2->[0][$l] = $S2->[$l]->sub(Otmt=>1,Otmo=>0,Otpz=>0,Otpo=>0,Otpt=>0);
  $M2->[1][$l] = $S2->[$l]->sub(Otmt=>0,Otmo=>1,Otpz=>0,Otpo=>0,Otpt=>0);
  $M2->[2][$l] = $S2->[$l]->sub(Otmt=>0,Otmo=>0,Otpz=>1,Otpo=>0,Otpt=>0);
  $M2->[3][$l] = $S2->[$l]->sub(Otmt=>0,Otmo=>0,Otpz=>0,Otpo=>1,Otpt=>0);
  $M2->[4][$l] = $S2->[$l]->sub(Otmt=>0,Otmo=>0,Otpz=>0,Otpo=>0,Otpt=>1);
}

if (!$output) {
  %T2 = (0,'O2m2',1,'O2m1',2,'O2_0',3,'O2_1',4,'O2_2');
  print "Crystal field parameters under rotation of phi=$phi about y-axis and theta=$theta about z-axis are given by:\n\n";
  print "Transformations relations for O2q:\n";
  for $l (0 .. 4) { 
    print "\[$T2{$l}\] = "; 
    $_ = $S2->[$l]; 
    s/\$Otmt/\{O2m2\}/g;
    s/\$Otmo/\{O2m1\}/g;
    s/\$Otpz/\{O2_0\}/g;
    s/\$Otpo/\{O2_1\}/g;
    s/\$Otpt/\{O2_2\}/g;
    print "$_\n"; 
  }
  print "\nRotation Matrix for O2q:\n";
  for $l (0 .. 4) { for $m (0 .. 4) {
    print $M2->[$l][$m]; print "\t"; } print "\n"; }
}
else {
  %T2 = (0,'O_2^{-2}',1,'O_2^{-1}',2,'O_2^0',3,'O_2^1',4,'O_2^2');

  open (outfile, ">$output") or die "$0: cannot open $output for output.";
  print outfile << "EOF";

\\documentclass [12pt,a4paper,notitlepage]{article}
\\usepackage{latexsym}
\\usepackage{pslatex}
\\setlength{\\unitlength}{1cm}

\\begin{document}

Crystal field parameters under rotation of phi=$phi\$^{o}\$ about y-axis and theta=$theta\$^{o}\$ about z-axis are given by:

Transformation relations for \$O_2^q\$ Steven's equivalent operators:

EOF
  
  for $l (0 .. 4) {
    $_ = $S2->[$l];
    s/sqrt\(([\/\.\d]*)\)/\\sqrt\{$1}/g;
    s/([\.\d]+)\/([\.\d]+)/\\frac\{$1\}\{$2\}/g;
    s/\$Otmt/\\left\[ O\_2\^\{-2\} \\right\]/g;
    s/\$Otmo/\\left\[ O\_2\^\{-1\} \\right\]/g;
    s/\$Otpz/\\left\[ O\_2\^0 \\right\]/g;
    s/\$Otpo/\\left\[ O\_2\^1 \\right\]/g;
    s/\$Otpt/\\left\[ O\_2\^2 \\right\]/g;
    s/\*//g;
    print outfile "\\begin\{equation\}\n";
    print outfile "\\left\\\{ $T2{$l} \\right\\\} = $_\n";
    print outfile "\\end\{equation\}\n";
  }

  print outfile "\nRotation Matrix for \$O_2^q\$:\n\\begin{equation}\n\\left(\n\\begin{array}{ccccc}\n";
  for $l (0 .. 4) { 
    for $m (0 .. 4) {
      $_ = $M2->[$l][$m];
      s/sqrt\(([\/\.\d]*)\)/\\sqrt\{$1}/g;
      s/\/(\\sqrt\{[\/\.\d]*\})/\\frac\{\}\{$1\}/g;
      s/([\.\d]+)\/([\.\d]+)/\\frac\{$1\}\{$2\}/g;       # + matches 1 or more times
      s/\*//g;
      if ($m == 4) { print outfile "$_ "; } else {
      print outfile "$_ & "; }
    } 
    print outfile "\\\\\n"; 
  }
  print outfile "\\end{array}\n\\right)\n\\end{equation}\n";

  if (!$z_int) {
    print outfile "\\end\{document\}\n";
  }

  close (out_file);
}
