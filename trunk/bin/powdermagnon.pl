#! /usr/bin/perl
use FileHandle;
use PDL;
use PDL::Slatec;

$pi=3.1415;

unless ($#ARGV>=2) 
{print " program to create mcdisp.ini from mcphas.j for the calculation of neutron powder spectra\n\n";
print " usage: powdermagnon 0.3 2 0.1 10\n\n";
print " meaning take mcphas.j, generate a reflection list and put it to mcdisp.ini \n";
print " 0.3 ....qmin   [1/A] minimal q vector\n";
print " 2   ....qmax   [1/A] maximal q vector\n";
print " 0.1 ....deltaq [1/A] stepwidth in q\n";
print " 10  ....number of steps in polar coordinate theta\n";
print "        (steps in fi fixed by dfi=dtheta*pi/(4*sin(theta))\n";

print "\n";
print "... then run mcdisp\n";
print "\n";
print ".... then restart powdermagnon by: powdermagnon -r results/mcdisp.qei 0.5 30 0.2 \n";
print " 0.5 ....Emin   [meV] minimal energy\n";
print " 0.5 ....Emax   [meV] maximal energy\n";
print " 0.5 ....deltaE [meV] energy stepwidth\n";
print " The powder average is generated - output is printed to the console (STDOUT)\n";
 exit 0;}

my ($qmin) = $ARGV[0];
my ($qmax) = $ARGV[1];

if($qmin=="-r")
{#read output file mcdisp.qei and create powder average 
 print "#Emin=".$Emin." meV   Emax=".$Emax." meV deltaE=".$deltaE." meV\n";
 print "#Ha[T] Hb[T] Hc[T] T[K] Q[A^-1] energy[meV] powderintensity [barn/sr/f.u.]   f.u.=crystallogrpaphic unit cell (r1xr2xr3)\n";
    my $h = new FileHandle;
    open($h,$qmax);

my ($Emin) = $ARGV[2];
my ($Emax) = $ARGV[3];
my ($deltaE) = $ARGV[4];
$n=int(($Emax-$Emin)/$deltaE);
my (@ints)=();$#ints=$n;

$qold=0;$counter=1;
$qhold=0;$qkold=0;$qlold=0;
while(<$h>)
 {next if /^\s*#/;
  $line=$_;
  $line=~s/D/E/g;@numbers=split(" ",$line);
  $q=$numbers[7];
  if($q!=$qold&&$qold!=0){
     for($i=0;$i<=$n;++$i){
     $E=$Emin+($i-0.5)*$deltaE;$ints[$i]/=$counter;
     print $numbers[0]." ".$numbers[1]." ".$numbers[2]." ".$numbers[3]." ";
     print $qold." ".$E." ".$ints[$i]."\n";}
     $counter=0;@ints=0;
               }
  $i=int(($numbers[8]-$Emin)/$deltaE);
  $ints[$i]+=$numbers[9];
  $qold=$q;
  $qh=$numbers[4];
  $qk=$numbers[5];
  $ql=$numbers[6];
  if($qh!=$qhold||$qk!=$qkold||$ql!=$qlold){++$counter;}
  $qhold=$qh;
  $qkold=$qk;
  $qlold=$ql;
 }
close $h;
exit;
}

my ($deltaq) = $ARGV[2];
my ($nn)=$ARGV[3];
$dtheta=$pi/$nn+0.000001;

my ($latt,$p) = getlattice("./mcphas.j");
my ($a,$b,$c,$nofatoms) = @{$latt};
 print "lattice: a=".$a." b=".$b." c=".$c."\n";
 print "creating file mcdisp.ini\n";
    my $h = new FileHandle;
    open($h,">./mcdisp.ini");
print $h "# file mcdisp.ini created by program powdermagnon\n";
print $h "#qmin=".$qmin." A^-1  qmax=".$qmax." A^-1 deltaq=".$deltaq." A^-1 no of thetasteps=".$nn."\n";
print $h "# lattice: a=".$a." b=".$b." c=".$c."\n";
print $h "#\n";
print $h "# minimum/maximum energy of dispersion used for \n";
print $h "#               - calculation of standard deviation from data (see manual)\n";
print $h "#               - energy range when using option -r (full calculation of chi)\n";
print $h "emin=-0.02\n";
print $h "emax=0.02\n";
print $h "kf=2\n";
print $h "#\n";
#loop q
for ($q=$qmin;$q<=$qmax;$q+=$deltaq){ 

#loop sphere
 $h1=0;
 $h2=0;
 $h3=$c*$q/2/$pi;

 print $h "$h1 $h2 $h3\n";

 for ($theta=$dtheta;$theta<$pi;$theta+=$dtheta){ 
  $dfi=$dtheta*$pi/4/sin($theta);
 for ($fi=0;$fi<2*$pi;$fi+=$dfi){ 
 
 $h1=$a*$q*sin($theta)*cos($fi)/2/$pi;
 $h2=$b*$q*sin($theta)*sin($fi)/2/$pi;
 $h3=$c*$q*cos($theta)/2/$pi;

 print $h "$h1 $h2 $h3\n";
 }}
 $h1=0;
 $h2=0;
 $h3=-$c*$q/2/$pi;

 print $h "$h1 $h2 $h3\n";

}
    close $h; 
 print "please edit mcdisp.ini, mcdisp.mf and run mcdisp\n";
 print "then restart powdermagnon by: powdermagnon -r results\mcdisp.qei Emin Emax deltaE\n";
   exit;
#-----------------------------------------------------------------------
# Get lattic data, reading it from file 
sub getlattice {
    my ($file) = @_;
    my $h = new FileHandle;
    my $n = 0;
#     my @xlist = ();
  # input data int piddle
  if(open($h,$file))
  {      while(<$h>)
     {#next if /^\s*#/;
      # detect a= b= c= ...
      if (/^#\s*\Qa=\E/){($a)=(/^#\s*\Qa=\E\s*([^\s]*)\s/);}
      if (/^.*\Qb=\E/){($b)=(/^.*\Qb=\E\s*([^\s]*)\s/);}
      if (/^.*\Qc=\E/){($c)=(/^.*\Qc=\E\s*([^\s]*)\s/);}

      if (/^.*\Qr1x=\E/){($r1x)=(/^.*\Qr1x=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr1y=\E/){($r1y)=(/^.*\Qr1y=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr1z=\E/){($r1z)=(/^.*\Qr1z=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr2x=\E/){($r2x)=(/^.*\Qr2x=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr2y=\E/){($r2y)=(/^.*\Qr2y=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr2z=\E/){($r2z)=(/^.*\Qr2z=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr3x=\E/){($r3x)=(/^.*\Qr3x=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr3y=\E/){($r3y)=(/^.*\Qr3y=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr3z=\E/){($r3z)=(/^.*\Qr3z=\E\s*([^\s]*)\s/);}
      if (/^.*\Qnofatoms=\E/){($nofatoms)=(/^.*\Qnofatoms=\E\s*([^\s]*)\s/);}

      if (/^.*\Qnofneighbours=\E/){++$n;
                                   ($nofneighbours[$n])=(/^.*\Qnofneighbours=\E\s*([^\s]*)\s/);
				   ($x[$n])=(/^.*\Qx=\E\s*([^\s]*)\s/);
				   ($y[$n])=(/^.*\Qy=\E\s*([^\s]*)\s/);
				   ($z[$n])=(/^.*\Qz=\E\s*([^\s]*)\s/);
				   ($gJ[$n])=(/^.*\QgJ=\E\s*([^\s]*)\s/);
				  }
     }
     close $h; 
     if ($n!=$nofatoms) {print STDOUT "Failed to read data file \"$file\": wrong number of atoms\n";
                         return undef;}
     return ([$a,$b,$c,$nofatoms],pdl [[$r1x,$r1y,$r1z],[$r2x,$r2y,$r2z],[$r3x,$r3y,$r3z]]);
    } else {
	print STDOUT "Warning: failed to read data file \"$file\"\n";
	return undef;
    }
}


