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

my ($f, $t, $i, $pi) = symbols(qw(f t i pi));
my ($T, $S, $s, $F) = symbols(qw(two seven seventy five));
my ($l, $m, $n);

my $A4 = [
	[2,   0,        0,      0,        0,      0,         0,     0,	     -2],
	[0, 1/sqrt($T), 0,      0,        0,      0,         0,	 1/sqrt($T),  0],
	[0,   0,    sqrt($S),   0,        0,      0,    -sqrt($S),  0,	      0],
	[0,   0,        0,   sqrt($S/$T), 0,   sqrt($S/$T),  0,     0,	      0],
	[0,   0,        0,      0,    2*sqrt($s), 0,         0,     0,	      0],
	[0,   0,        0,   sqrt($S/$T), 0,  -sqrt($S/$T),  0,     0,	      0],
	[0,   0,    sqrt($S),   0,	  0,      0,     sqrt($S),  0,	      0],
	[0, 1/sqrt($T), 0,      0,	  0,      0,         0,	-1/sqrt($T),  0],
	[2,   0,        0,      0,	  0,      0,         0,     0,	      2]
     ];

my $invA4 = [
	[ 1/4,     0,       0,       0,        0,      0,         0,     0,       1/4],
	[   0, 1/sqrt($T),  0,       0,        0,      0,         0,  1/sqrt($T),  0],
	[   0,     0,  0.5/sqrt($S), 0,        0,      0,  0.5/sqrt($S), 0,        0],
	[   0,     0,       0,  1/sqrt($T*$S), 0,  1/sqrt($T*$S), 0,     0,        0],
	[   0,     0,       0,       0,  0.5/sqrt($s), 0,         0,     0,        0],
	[   0,     0,       0,  1/sqrt($T*$S), 0, -1/sqrt($T*$S), 0,     0,        0],
	[   0,     0, -0.5/sqrt($S), 0,        0,      0,  0.5/sqrt($S), 0,        0],
	[   0, 1/sqrt($T),  0,       0,        0,      0,         0, -1/sqrt($T),  0],
	[-1/4,     0,       0,       0,        0,      0,         0,     0,       1/4]
      ];

# --------------------------------------- O4 ------------------------------------- #

my $D4p0p4 =  sqrt($s)/16   * sin($t)**4;
my $D4p0p3 =  sqrt($S*$T)/4 * sin($t)**3 * cos($t);
my $D4p0p2 =  sqrt($F*$T)/8 * sin($t)**2 * (7*cos($t)**2 - 1);
my $D4p0p1 =  sqrt($F)/4    * sin($t)    * cos($t) * (7*cos($t)**2 - 3);
my $D4p0p0 =  (1/8)         * (35*cos($t)**4 - 30*cos($t)**2 + 3);

my $D4p1p0 =  exp( $i*$f) * sqrt($F)/4    * sin($t)     * cos($t)     * (7*cos($t)**2-3);
my $D4p1p1 =  exp( $i*$f) * (1/8)         * (1+cos($t)) * (28*cos($t)**3 - 21*cos($t)**2 - 6*cos($t) + 3);
my $D4p1m1 =  exp(-$i*$f) * (1/8)         * (1-cos($t)) * (28*cos($t)**3 + 21*cos($t)**2 - 6*cos($t) - 3);
my $D4p1p2 =  exp( $i*$f) * sqrt($T)/8    * sin($t)     * (1+cos($t)) * (14*cos($t)**2 - 7*cos($t) - 1);
my $D4p1m2 = -exp(-$i*$f) * sqrt($T)/8    * sin($t)     * (1-cos($t)) * (14*cos($t)**2 + 7*cos($t) - 1);
my $D4p1p3 =  exp( $i*$f) * sqrt($S)/8    * sin($t)**2  * (1+cos($t)) * ( 4*cos($t) - 1);
my $D4p1m3 =  exp(-$i*$f) * sqrt($S)/8    * sin($t)**2  * (1-cos($t)) * ( 4*cos($t) + 1);
my $D4p1p4 =  exp( $i*$f) * sqrt($T*$S)/8 * sin($t)**3  * (1+cos($t));
my $D4p1m4 = -exp(-$i*$f) * sqrt($T*$S)/8 * sin($t)**3  * (1-cos($t));

my $D4p2p0 =  exp( 2*$i*$f) * sqrt($T*$F)/8 * sin($t)**2 * (7*cos($t)**2 - 1);
my $D4p2p1 =  exp( 2*$i*$f) * sqrt($T)/8    * sin($t)    * (1+cos($t)) * (14*cos($t)**2 - 7*cos($t) - 1);
my $D4p2m1 =  exp(-2*$i*$f) * sqrt($T)/8    * sin($t)    * (1-cos($t)) * (14*cos($t)**2 + 7*cos($t) - 1);
my $D4p2p2 =  exp( 2*$i*$f) * (1/4)         * (1+cos($t))**2 * (7*cos($t)**2 - 7*cos($t) + 1);
my $D4p2m2 =  exp(-2*$i*$f) * (1/4)         * (1-cos($t))**2 * (7*cos($t)**2 + 7*cos($t) + 1);
my $D4p2p3 =  exp( 2*$i*$f) * sqrt($T*$S)/8 * sin($t)    * (1+cos($t))**2 * (2*cos($t) - 1);
my $D4p2m3 = -exp(-2*$i*$f) * sqrt($T*$S)/8 * sin($t)    * (1-cos($t))**2 * (2*cos($t) + 1);
my $D4p2p4 =  exp( 2*$i*$f) * sqrt($S)/8    * sin($t)**2 * (1+cos($t))**2;
my $D4p2m4 =  exp(-2*$i*$f) * sqrt($S)/8    * sin($t)**2 * (1-cos($t))**2;

my $D4p3p0 =  exp( 3*$i*$f) * sqrt($F*$S)/4 * sin($t)**3 * cos($t);
my $D4p3p1 =  exp( 3*$i*$f) * sqrt($S)/8    * sin($t)**2 * (1+cos($t)) * (4*cos($t) - 1);
my $D4p3m1 =  exp(-3*$i*$f) * sqrt($S)/8    * sin($t)**2 * (1-cos($t)) * (4*cos($t) + 1);
my $D4p3p2 =  exp( 3*$i*$f) * sqrt($T*$S)/8 * sin($t)   * (1+cos($t))**2 * (2*cos($t) - 1);
my $D4p3m2 =  exp(-3*$i*$f) * sqrt($T*$S)/8 * sin($t)   * (1-cos($t))**2 * (2*cos($t) + 1);
my $D4p3p3 =  exp( 3*$i*$f) * (1/8)         * (1+cos($t))**3 * (4*cos($t) - 3);
my $D4p3m3 =  exp(-3*$i*$f) * (1/8)         * (1-cos($t))**3 * (4*cos($t) + 3);
my $D4p3p4 =  exp( 3*$i*$f) * sqrt($T)/8    * sin($t)   * (1+cos($t))**3;
my $D4p3m4 = -exp(-3*$i*$f) * sqrt($T)/8    * sin($t)   * (1-cos($t))**3;

my $D4p4p0 =  exp( 4*$i*$f) * sqrt($s)/16   * sin($t)**4;
my $D4p4p1 =  exp( 4*$i*$f) * sqrt($T*$S)/8 * sin($t)**3 * (1+cos($t));
my $D4p4m1 =  exp(-4*$i*$f) * sqrt($T*$S)/8 * sin($t)**3 * (1-cos($t));
my $D4p4p2 =  exp( 4*$i*$f) * sqrt($S)/8    * sin($t)**2 * (1+cos($t))**2;
my $D4p4m2 =  exp(-4*$i*$f) * sqrt($S)/8    * sin($t)**2 * (1-cos($t))**2;
my $D4p4p3 =  exp( 4*$i*$f) * sqrt($T)/8    * sin($t)   * (1+cos($t))**3;
my $D4p4m3 =  exp(-4*$i*$f) * sqrt($T)/8    * sin($t)   * (1-cos($t))**3;
my $D4p4p4 =  exp( 4*$i*$f) * (1/16)        * (1+cos($t))**4;
my $D4p4m4 =  exp(-4*$i*$f) * (1/16)        * (1-cos($t))**4;

my $D4 = [
             [ $D4p4p4,  $D4p4p3,  $D4p4p2,  $D4p4p1,  $D4p4p0,  $D4p4m1,  $D4p4m2,  $D4p4m3,  $D4p4m4], 
             [-$D4p3p4,  $D4p3p3,  $D4p3p2,  $D4p3p1,  $D4p3p0,  $D4p3m1,  $D4p3m2,  $D4p3m3, -$D4p3m4], 
             [ $D4p2p4, -$D4p2p3,  $D4p2p2,  $D4p2p1,  $D4p2p0,  $D4p2m1,  $D4p2m2, -$D4p2m3,  $D4p2m4], 
             [-$D4p1p4,  $D4p1p3, -$D4p1p2,  $D4p1p1,  $D4p1p0,  $D4p1m1, -$D4p1m2,  $D4p1m3, -$D4p1m4], 
             [ $D4p0p4, -$D4p0p3,  $D4p0p2, -$D4p0p1,  $D4p0p0,  $D4p0p1,  $D4p0p2,  $D4p0p3,  $D4p0p4], 
             [ $D4p1m4,  $D4p1m3,  $D4p1m2,  $D4p1m1, -$D4p1p0,  $D4p1p1,  $D4p1p2,  $D4p1p3,  $D4p1p4], 
             [ $D4p2m4,  $D4p2m3,  $D4p2m2, -$D4p2m1,  $D4p2p0, -$D4p2p1,  $D4p2p2,  $D4p2p3,  $D4p2p4], 
             [ $D4p3m4,  $D4p3m3, -$D4p3m2,  $D4p3m1, -$D4p3p0,  $D4p3p1, -$D4p3p2,  $D4p3p3,  $D4p3p4], 
             [ $D4p4m4, -$D4p4m3,  $D4p4m2, -$D4p4m1,  $D4p4p0, -$D4p4p1,  $D4p4p2, -$D4p4p3,  $D4p4p4]
         ];

my ($O4m4, $O4m3, $O4m2, $O4m1, $O4p0) = symbols(qw(Ofmf OfmT Ofmt Ofmo Ofpz));
my ($O4p4, $O4p3, $O4p2, $O4p1) = symbols(qw(Ofpf OfpT Ofpt Ofpo));

# ------------------------------------ Finished !! ------------------------------- #

# The operator in the rotated frame. Rx(y) is equivalent to {Oxy} in Rudowicz notation
# and Ox_y, Oxmy is equivalent to [Oxy], [OxyM] in Rudowicz

my $B4 = [ [$O4m4],[$O4m3],[$O4m2],[$O4m1],[$O4p0],[$O4p1],[$O4p2],[$O4p3],[$O4p4] ];

my $tmp1 = mmult($invA4,$B4);
my $tmp2 = mmult($D4, $tmp1);
my $R4 = mmult($A4, $tmp2);

my $S4 = [];
my $M4 = [];
for $l (0 .. 8) { 
  $S4->[$l] = $R4->[$l][0]->sub(f=>0, t=>$pi/2);
  $S4->[$l] = $S4->[$l]->sub(seven=>7,two=>2,seventy=>70,five=>5); 
  $M4->[0][$l] = $S4->[$l]->sub(Ofmf=>1,OfmT=>0,Ofmt=>0,Ofmo=>0,Ofpz=>0,Ofpo=>0,Ofpt=>0,OfpT=>0,Ofpf=>0);
  $M4->[1][$l] = $S4->[$l]->sub(Ofmf=>0,OfmT=>1,Ofmt=>0,Ofmo=>0,Ofpz=>0,Ofpo=>0,Ofpt=>0,OfpT=>0,Ofpf=>0);
  $M4->[2][$l] = $S4->[$l]->sub(Ofmf=>0,OfmT=>0,Ofmt=>1,Ofmo=>0,Ofpz=>0,Ofpo=>0,Ofpt=>0,OfpT=>0,Ofpf=>0);
  $M4->[3][$l] = $S4->[$l]->sub(Ofmf=>0,OfmT=>0,Ofmt=>0,Ofmo=>1,Ofpz=>0,Ofpo=>0,Ofpt=>0,OfpT=>0,Ofpf=>0);
  $M4->[4][$l] = $S4->[$l]->sub(Ofmf=>0,OfmT=>0,Ofmt=>0,Ofmo=>0,Ofpz=>1,Ofpo=>0,Ofpt=>0,OfpT=>0,Ofpf=>0);
  $M4->[5][$l] = $S4->[$l]->sub(Ofmf=>0,OfmT=>0,Ofmt=>0,Ofmo=>0,Ofpz=>0,Ofpo=>1,Ofpt=>0,OfpT=>0,Ofpf=>0);
  $M4->[6][$l] = $S4->[$l]->sub(Ofmf=>0,OfmT=>0,Ofmt=>0,Ofmo=>0,Ofpz=>0,Ofpo=>0,Ofpt=>1,OfpT=>0,Ofpf=>0);
  $M4->[7][$l] = $S4->[$l]->sub(Ofmf=>0,OfmT=>0,Ofmt=>0,Ofmo=>0,Ofpz=>0,Ofpo=>0,Ofpt=>0,OfpT=>1,Ofpf=>0);
  $M4->[8][$l] = $S4->[$l]->sub(Ofmf=>0,OfmT=>0,Ofmt=>0,Ofmo=>0,Ofpz=>0,Ofpo=>0,Ofpt=>0,OfpT=>0,Ofpf=>1);
}

if (!$output) {
  %T4 = (0,'O4m4',1,'O4m3',2,'O4m2',3,'O4m1',4,'O4_0',5,'O4_1',6,'O4_2',7,'O4_3',8,'O4_4');
  print "Crystal field parameters under rotation of phi=$phi about y-axis and theta=$theta about z-axis are given by:\n\n";
  print "Transformations relations for O4q:\n";
  for $l (0 .. 8) { 
    print "\[$T4{$l}\] = "; 
    $_ = $S4->[$l]; 
    s/\$Ofmf/\{O4m4\}/g;
    s/\$OfmT/\{O4m3\}/g;
    s/\$Ofmt/\{O4m2\}/g;
    s/\$Ofmo/\{O4m1\}/g;
    s/\$Ofpz/\{O4_0\}/g;
    s/\$Ofpo/\{O4_1\}/g;
    s/\$Ofpt/\{O4_2\}/g;
    s/\$OfpT/\{O4_3\}/g;
    s/\$Ofpf/\{O4_4\}/g;
    print "$_\n"; 
  }
  print "\nRotation Matrix for O4q:\n";
  for $l (0 .. 8) { for $m (0 .. 8) {
    print $M4->[$l][$m]; print "\t"; } print "\n"; }
}
else {
  %T4 = (0,'O_4^{-4}',1,'O_4^{-4}',2,'O_4^{-2}',3,'O_4^{-1}',4,'O_4^0',5,'O_4^1',6,'O_4^2',7,'O_4^3',8,'O_4^4');

  if (!$z_int) {
    open (outfile, ">$output") or die "$0: cannot open $output for output.";
    print outfile << "EOF";

\\documentclass [12pt,a4paper,notitlepage]{article}
\\usepackage{latexsym}
\\usepackage{pslatex}
\\setlength{\\unitlength}{1cm}

\\begin{document}

Crystal field parameters under rotation of phi=$phi\$^{o}\$ about y-axis and theta=$theta\$^{o}\$ about z-axis are given by:

Transformation relations for \$O_4^q\$ Steven's equivalent operators:

EOF
  
  }
  else {
    open (outfile, ">>$output") or die;
    print outfile "\nTransformation relations for \$O_4^q\$ Steven's equivalent operators:\n";
  }

  for $l (0 .. 8) {
    $_ = $S4->[$l];
    s/sqrt\(([\/\.\d]*)\)/\\sqrt\{$1}/g;
    s/([\.\d]+)\/([\.\d]+)/\\frac\{$1\}\{$2\}/g;
    s/\$Ofmf/\\left\[ O\_4\^\{-4\} \\right\]/g;
    s/\$OfmT/\\left\[ O\_4\^\{-3\} \\right\]/g;
    s/\$Ofmt/\\left\[ O\_4\^\{-2\} \\right\]/g;
    s/\$Ofmo/\\left\[ O\_4\^\{-1\} \\right\]/g;
    s/\$Ofpz/\\left\[ O\_4\^0 \\right\]/g;
    s/\$Ofpo/\\left\[ O\_4\^1 \\right\]/g;
    s/\$Ofpt/\\left\[ O\_4\^2 \\right\]/g;
    s/\$OfpT/\\left\[ O\_4\^3 \\right\]/g;
    s/\$Ofpf/\\left\[ O\_4\^4 \\right\]/g;
    s/\*//g;
    print outfile "\\begin\{equation\}\n";
    print outfile "\\left\\\{ $T4{$l} \\right\\\} = $_\n";
    print outfile "\\end\{equation\}\n";
  }

  print outfile "\nRotation Matrix for \$O_4^q\$:\n\\begin{equation}\n\\left(\n\\begin{array}{ccccccccc}\n";
  for $l (0 .. 8) { 
    for $m (0 .. 8) {
      $_ = $M4->[$l][$m];
      s/sqrt\(([\/\.\d]*)\)/\\sqrt\{$1}/g;
      s/\/(\\sqrt\{[\/\.\d]*\})/\\frac\{\}\{$1\}/g;
      s/([\.\d]+)\/([\.\d]+)/\\frac\{$1\}\{$2\}/g;       # + matches 1 or more times
      s/\*//g;
      if ($m == 8) { print outfile "$_ "; } else {
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