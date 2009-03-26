#! /usr/bin/perl
use FileHandle;
use PDL;
use PDL::Slatec;

unless ($#ARGV>=0) 
{print " program to create table with neighbors from mcphas.j \n\n";
print " usage: makenn 23.3 [options] \n\n";
print " meaning take mcphas.j, generate all neighbors within sphere of 23.3A and\n";
print " put them into makenn.j, in interaction columns put the classical dipole interaction\n";
print " the output values are sorted by ascending distance\n";
print " option -rkky A(meV) kf(1/A) calculates the rkky interaction\n";
print "              according to J(R)=A.cos(2.kf.R)/(2.kf.R)^3\n";
print "              scaling A<0, kf should be the Fermi wavevector (usually\n";
print "              between 0.3-2.5 A^-1 depending on the electrondensity^0.333)\n";
print " option -rkky3d A(meV) ka(1/A) kb(1/A) kc(1/A) calculates the rkky interaction\n";
print "              according to J(R)=A.cos(2.kfR)/(2.kfR)^3\n";
print "              scaling A<0, kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2)\n";
print " option -rkkz A(meV) kf(1/A) calculates the rkky interaction\n";
print "              according to J(R)=A [sin(2.kf.R)-2.kf.R.cos(2.kf.R)]/(2.kf.R)^4\n";
print "              scaling A>0, kf should be the Fermi wavevector\n";
print " option -rkkz3d A(meV) ka(1/A) kb(1/A) kc(1/A)  calculates the rkky interaction\n";
print "              according to J(R)=A [sin(2.kfR)-2.kfR.cos(2.kfR)]/(2.kfR)^4\n";
print "              scaling A>0, kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2)\n";
print " option -kaneyoshi A(meV) D(A) alpha  calculates the kaneyoshi\n";
print  "             parametrization for the Bethe-Slater\n";
print "              curve: J(R)= A [-(R/D)^2+(R/D)^4].exp[-alpha.(R/D)^2]\n";
print "              with D corresponding to the orbital radius\n";
print "              the exponential alpha is conveniently put to  about 1\n";
print " option -kaneyoshi3d A(meV) Da(A) Db(A) Dc(A) alpha  calculates the 3d-kaneyoshi\n";
print  "             parametrization for the Bethe-Slater\n";
print "              curve: J(R)= A [-(RD)^2+(RD)^4].exp[-alpha.(RD)^2]\n";
print "              with RD=sqrt(Ra^2/Da^2+Rb^2/Db^2+Rc^2/Dc^2)\n";
print "              the exponential alpha is conveniently put to  about 1\n";
print " option -d puts to the last column the distance of the neighbors (A)\n";
 exit 0;}

my ($rmax) = $ARGV[0];
$rkky=0;$calcdist=0;
shift @ARGV; 
$_=$ARGV[0];
if(/-rkky3d/)
  {$rkky=4;shift @ARGV; $scale=$ARGV[0];shift @ARGV;
   $ka=$ARGV[0];shift @ARGV;  $kb=$ARGV[0];shift @ARGV;   $kc=$ARGV[0];shift @ARGV;
   print "calculating RKKY interaction J(R)=A.cos(2.kfR)/(2.kfR)^3 for scale A=$scale meV and\n";
   print "kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2) with ka=$ka A^-1 kb=$kb A^-1 kc=$kc A^-1\n";}
elsif(/-kaneyoshi3d/)
  {$rkky=5;shift @ARGV; $scale=$ARGV[0];shift @ARGV;
   $Da=$ARGV[0];shift @ARGV;   $Db=$ARGV[0];shift @ARGV;
   $Dc=$ARGV[0];shift @ARGV;$aa=$ARGV[0];shift @ARGV;
   print "calculating kaneyoshi parametrization for the Bethe-Slater curve\n";
   print "J(R)= A [-(RD)^2+(RD)^4].exp[-alpha.(RD)^2] for scale A=$scale meV \n";
   print "with RD=sqrt(Ra^2/Da^2+Rb^2/Db^2+Rc^2/Dc^2) with Da=$Da A^-1 Db=$Db A^-1 Dc=$Dc A^-1  and alpha=$aa\n";}
elsif(/-rkkz3d/)
  {$rkky=6;shift @ARGV; $scale=$ARGV[0];shift @ARGV;
   $ka=$ARGV[0];shift @ARGV;  $kb=$ARGV[0];shift @ARGV;   $kc=$ARGV[0];shift @ARGV;
   print "calculating RKKY interaction J(R)=A [sin(2.kf.R)-2.kf.R.cos(2.kf.R)]/(2.kf.R)^4 for scale A=$scale meV\n";
   print "kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2) with ka=$ka A^-1 kb=$kb A^-1 kc=$kc A^-1\n";}
elsif(/-rkky/)
  {$rkky=1;shift @ARGV; $scale=$ARGV[0];shift @ARGV;$kf=$ARGV[0];shift @ARGV;
   print "calculating RKKY interaction J(R)=A.cos(2.kf.R)/(2.kf.R)^3 for scale A=$scale meV and kf=$kf A^-1\n";}
elsif(/-kaneyoshi/)
  {$rkky=2;shift @ARGV; $scale=$ARGV[0];shift @ARGV;$D=$ARGV[0];shift @ARGV;$aa=$ARGV[0];shift @ARGV;
   print "calculating kaneyoshi parametrization for the Bethe-Slater curve\n";
   print "J(R)= A [-(R/D)^2+(R/D)^4].exp[-alpha.(R/D)^2] for scale A=$scale meV D=$D A^-1  alpha=$aa\n";}
elsif(/-rkkz/)
  {$rkky=3;shift @ARGV; $scale=$ARGV[0];shift @ARGV;$kf=$ARGV[0];shift @ARGV;
   print "calculating RKKY interaction J(R)=A [sin(2.kf.R)-2.kf.R.cos(2.kf.R)]/(2.kf.R)^4 for scale A=$scale meV and kf=$kf A^-1\n";}
$_=$ARGV[0];
if(/-d/)
  {$calcdist=1;print "putting distance of neighbors (A) to last column of makenn.j\n";}

my ($latt,$p) = getlattice("./mcphas.j");
my ($a,$b,$c,$nofatoms) = @{$latt};
 print "rmax=".$rmax." A\n";
 print "a=".$a." b=".$b." c=".$c."\n";
 print "primitive lattice[abc]:".$p."\n";
 $t=$p->slice("0,"); $t*=$a;
 $t=$p->slice("1,"); $t*=$b;
 $t=$p->slice("2,"); $t*=$c;
 print "primitive lattice[A]:".$p."\n";

    $r=0;
# determine $nmin,$nmax by looking at a cube with side 3rmax
     $inv=matinv(transpose($p)); #invert primitive lattice
# print "inverted primitive lattice[A]:".$inv."\n";
     #loop all corner points
  for ($i1=-1;$i1<=1;$i1+=2){
  for ($i2=-1;$i2<=1;$i2+=2){
  for ($i3=-1;$i3<=1;$i3+=2){
    $n=inner($inv , pdl[$i1*$rmax*1.5,$i2*$rmax*1.5,$i3*$rmax*1.5]); 
    if (($n->at(0))<$n1min) {$n1min=int($n->at(0))-1;}
    if (($n->at(1))<$n2min) {$n2min=int($n->at(1))-1;}
    if (($n->at(2))<$n3min) {$n3min=int($n->at(2))-1;}
    if (($n->at(0))>$n1max) {$n1max=int($n->at(0))+1;}
    if (($n->at(1))>$n2max) {$n2max=int($n->at(1))+1;}
    if (($n->at(2))>$n3max) {$n3max=int($n->at(2))+1;}
#    print"corner $i1 $i2 $i3 coordinates in prim bases:$n\n";
 #   print "$n1min to $n1max, $n2min to $n2max, $n3min to $n3max\n";
  }}}
print "$n1min to $n1max, $n2min to $n2max, $n3min to $n3max\n";

     #initialize output file results/makenn.j
 my ($h,$l)=printlattice("./mcphas.j",$n,$rn,$xn,$yn,$zn,$Jaa,$Jbb,$Jcc,$Jab,$Jba,$Jac,$Jca,$Jbc,$Jcb,">./results/makenn.j");

 for ($nnn=1;$nnn<=$nofatoms;++$nnn)    
 {   my $gJ=$gJ[$nnn];
     my ($rn)=new PDL ();
     my ($an)=new PDL ();

     my ($xn)=new PDL();    
     my ($yn)=new PDL();    
     my ($zn)=new PDL();    

     my ($Jaa)=new PDL();    
     my ($Jab)=new PDL();    
     my ($Jac)=new PDL();    

     my ($Jba)=new PDL();    
     my ($Jbb)=new PDL();    
     my ($Jbc)=new PDL();    

     my ($Jca)=new PDL();    
     my ($Jcb)=new PDL();    
     my ($Jcc)=new PDL();    


  for ($n1=$n1min;$n1<=$n1max;++$n1){ 
  for ($n2=$n2min;$n2<=$n2max;++$n2){ 
  for ($n3=$n3min;$n3<=$n3max;++$n3){  
  for ($nz=1;$nz<=$nofatoms;++$nz){  
   $rvec=pdl [$a*($x[$nz]-$x[$nnn]),$b*($y[$nz]-$y[$nnn]),$c*($z[$nz]-$z[$nnn])];
   $rvec+=$n1*$p->slice(",(0)")+$n2*$p->slice(",(1)")+$n3*$p->slice(",(2)");
   $rr=inner($rvec, $rvec);
   $r=sqrt($rr->at());
   $xx=$rvec->at(0)/$a;
   $yy=$rvec->at(1)/$b;
   $zz=$rvec->at(2)/$c;
#   if ($r<=$rmax && $r>0 && $x>=0 && $y>=0 && $z>=0){#save neighbour jjj format
   if ($r<=$rmax && $r>0){#save neighbour j format
#       print $r."   ($xx,$yy,$zz,$n1,$n2,$n3,$nz\n";
    $an=$an->append( pdl ([$nz]));
    $rn=$rn->append( pdl ([$r]));
    $xn=$xn->append( pdl ([$xx]));
    $yn=$yn->append( pdl ([$yy]));
    $zn=$zn->append( pdl ([$zz]));
    my ($interaction) = getinteraction($gJ,$r,$xx*$a,$yy*$b,$zz*$c);
    my ($jaa,$jab,$jac,$jba,$jbb,$jbc,$jca,$jcb,$jcc) = @{$interaction};

    $Jaa=$Jaa->append( pdl ([$jaa]));
    $Jab=$Jab->append( pdl ([$jab]));
    $Jac=$Jac->append( pdl ([$jac]));
    $Jba=$Jba->append( pdl ([$jba]));
    $Jbb=$Jbb->append( pdl ([$jbb]));
    $Jbc=$Jbc->append( pdl ([$jbc]));
    $Jca=$Jca->append( pdl ([$jca]));
    $Jcb=$Jcb->append( pdl ([$jcb]));
    $Jcc=$Jcc->append( pdl ([$jcc]));                         }
    }}}}  
   $n= qsorti($rn);
   $nofneighbours[$nnn]=(($rn->dims)[0]-1);
   printneighbourlist($h,$l,$nofneighbours[$nnn],$n,$an,$rn,$xn,$yn,$zn,$Jaa,$Jbb,$Jcc,$Jab,$Jba,$Jac,$Jca,$Jbc,$Jcb);   
 }

 endprint($h,$l);   
 print "file results/makenn.j created\n";
   exit;
#-----------------------------------------------------------------------

sub getinteraction {
   my ($gJ,$r,$rx,$ry,$rz)=@_;
   # calculate classical dipole interaction
if ($rkky==1)
{$jaa = $scale*cos(2*$kf*$r)/8/$kf/$kf/$kf/$r/$r/$r;
 $jbb =$jaa;$jcc =$jaa;$jab = 0;$jbc =0;$jac =0;$jba=0;$jcb=0;$jca=0;
}
elsif ($rkky==4)
{$kfr=sqrt($ka*$ka*$rx*$rx+$kb*$kb*$ry*$ry+$kc*$kc*$rz*$rz);
 $jaa = $scale*cos(2*$kfr)/8/$kfr/$kfr/$kfr;
 $jbb =$jaa;$jcc =$jaa;$jab = 0;$jbc =0;$jac =0;$jba=0;$jcb=0;$jca=0;
}
elsif ($rkky==3)
{$jaa = $scale*(sin(2*$kf*$r)-2*$kf*$r*cos(2*$kf*$r))/16/$kf/$kf/$kf/$kf/$r/$r/$r/$r;
 $jbb =$jaa;$jcc =$jaa;$jab = 0;$jbc =0;$jac =0;$jba=0;$jcb=0;$jca=0;
}
elsif ($rkky==6)
{$kfr=sqrt($ka*$ka*$rx*$rx+$kb*$kb*$ry*$ry+$kc*$kc*$rz*$rz);
 $jaa = $scale*(sin(2*$kfr)-2*$kfr*cos(2*$kfr))/16/$kfr/$kfr/$kfr/$kfr;
 $jbb =$jaa;$jcc =$jaa;$jab = 0;$jbc =0;$jac =0;$jba=0;$jcb=0;$jca=0;
}
elsif($rkky==2)
{my ($xx)=$r*$r/$D/$D;
$jaa= $scale*(-$xx+$xx*$xx)*exp(-$aa*$xx);
$jbb =$jaa;$jcc =$jaa;$jab = 0;$jbc =0;$jac =0;$jba=0;  $jcb=0;$jca=0;
}
elsif($rkky==5)
{my ($xx)=$rx*$rx/$Da/$Da+$ry*$ry/$Db/$Db+$rz*$rz/$Dc/$Dc;
$jaa= $scale*(-$xx+$xx*$xx)*exp(-$aa*$xx);
$jbb =$jaa;$jcc =$jaa;$jab = 0;$jbc =0;$jac =0;$jba=0;  $jcb=0;$jca=0;
}
else
{
my $c = $gJ * $gJ * .927405 * .927405 / 16.02183;  #[meV A^3]

$jaa = $c * (3 * $rx * $rx -$r *$r) /$r /$r /$r /$r /$r;
$jbb = $c * (3 * $ry * $ry -$r *$r) /$r /$r /$r /$r /$r;
$jcc = $c * (3 * $rz * $rz -$r *$r) /$r /$r /$r /$r /$r;

$jab = $c * (3 * $rx * $ry) /$r /$r /$r /$r /$r;
$jbc = $c * (3 * $ry * $rz) /$r /$r /$r /$r /$r;
$jac = $c * (3 * $rx * $rz) /$r /$r /$r /$r /$r;
   $jba=$jab;   $jcb=$jbc;   $jca=$jac;
}     
 return ([$jaa,$jab,$jac,$jba,$jbb,$jbc,$jca,$jcb,$jcc]);  
}


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
      if ($a==0&&/^#\s*\Qa=\E/){($a)=(/^#\s*\Qa=\E\s*([^\s]*)\s/);}
      if ($b==0&&/^.*\Qb=\E/){($b)=(/^.*\Qb=\E\s*([^\s]*)\s/);}
      if ($c==0&&/^.*\Qc=\E/){($c)=(/^.*\Qc=\E\s*([^\s]*)\s/);}

      if (/^.*\Qr1x=\E/){($r1x)=(/^.*\Qr1x=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr1y=\E/){($r1y)=(/^.*\Qr1y=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr1z=\E/){($r1z)=(/^.*\Qr1z=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr2x=\E/){($r2x)=(/^.*\Qr2x=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr2y=\E/){($r2y)=(/^.*\Qr2y=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr2z=\E/){($r2z)=(/^.*\Qr2z=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr3x=\E/){($r3x)=(/^.*\Qr3x=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr3y=\E/){($r3y)=(/^.*\Qr3y=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr3z=\E/){($r3z)=(/^.*\Qr3z=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr1a=\E/){($r1x)=(/^.*\Qr1a=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr1b=\E/){($r1y)=(/^.*\Qr1b=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr1c=\E/){($r1z)=(/^.*\Qr1c=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr2a=\E/){($r2x)=(/^.*\Qr2a=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr2b=\E/){($r2y)=(/^.*\Qr2b=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr2c=\E/){($r2z)=(/^.*\Qr2c=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr3a=\E/){($r3x)=(/^.*\Qr3a=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr3b=\E/){($r3y)=(/^.*\Qr3b=\E\s*([^\s]*)\s/);}
      if (/^.*\Qr3c=\E/){($r3z)=(/^.*\Qr3c=\E\s*([^\s]*)\s/);}
      if (/^.*\Qnofatoms=\E/){($nofatoms)=(/^.*\Qnofatoms=\E\s*([^\s]*)\s/);}

      if (/^.*\Qnofneighbours=\E/){++$n;
                                   ($nofneighbours[$n])=(/^.*\Qnofneighbours=\E\s*([^\s]*)\s/);
				   if (/^.*\Qx=\E/){($x[$n])=(/^.*\Qx=\E\s*([^\s]*)\s/);}
				   if (/^.*\Qy=\E/){($y[$n])=(/^.*\Qy=\E\s*([^\s]*)\s/);}
				   if (/^.*\Qz=\E/){($z[$n])=(/^.*\Qz=\E\s*([^\s]*)\s/);}
				   if (/^.*\Qda=\E/){($x[$n])=(/^.*\Qda=\E\s*([^\s]*)\s/);}
				   if (/^.*\Qdb=\E/){($y[$n])=(/^.*\Qdb=\E\s*([^\s]*)\s/);}
				   if (/^.*\Qdc=\E/){($z[$n])=(/^.*\Qdc=\E\s*([^\s]*)\s/);}
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

sub printlattice {
  my ($filein,$n,$rn,$xn,$yn,$zn,$Jaa,$Jbb,$Jcc,$Jab,$Jba,$Jac,$Jca,$Jbc,$Jcb,$fileout)=@_;   
#   print $rn->index($n)."\n";
    my $h = new FileHandle;
    my $l = new FileHandle;
    open($l,$fileout);
    open($h,$filein);
     while(<$h>)
     {#next if /^\s*#/;
      $text=$_;
      last if /^#.*\Q**********\E/;
      print $l ($text);      
     }
 return ($h,$l);
}

sub printneighbourlist {
  my ($h,$l,$nofn,$n,$an,$rn,$xn,$yn,$zn,$Jaa,$Jbb,$Jcc,$Jab,$Jba,$Jac,$Jca,$Jbc,$Jcb)=@_;   
     if ($nofn=="-1"){$nofn="0";}
     print $l ("#*************************************************************************\n");
     my $stopheader=0;
     my $stopprint=0;
     while(<$h>)
     {#next if /^\s*#/;
      last if /^#.*\Q**********\E/;
      $text=$_;
     if (/^.*\Qnofneighbours=\E/){$text=~s!\Qnofneighbours=\E\s*\d+!nofneighbours=$nofn!;$stopheader=1;}
     if (/^.*\Qdiagonalexchange=\E/){
      if ($rkky>=1)
       {$text=~s!\Qdiagonalexchange=\E\s*\d+!diagonalexchange=1!;$stopheader=1;}
      else
       {$text=~s!\Qdiagonalexchange=\E\s*\d+!diagonalexchange=0!;$stopheader=1;}
      }

      if ($stopprint==0){print $l ($text);}
      $stopprint=$stopheader;
     }
     print $l ("#da[a]    db[b]     dc[c]       Jaa[meV]  Jbb[meV]  Jcc[meV]  Jab[meV]  Jba[meV]  Jac[meV]  Jca[meV]  Jbc[meV]  Jcb[meV]\n");
if ($rkky==1)
{print $l ("# it follows output of RKKY interaction according to J(R)=A.cos(2.kf.R)/(2.kf.R)^3 with A=$scale meV and kf=$kf A^-1 generated by makenn\n");}
elsif ($rkky==3)
{print $l ("# it follows output of RKKY interaction according to J(R)=A.[sin(2.kf.R)-2.kf.R.cos(2.kf.R)]/(2.kf.R)^4 with A=$scale meV and kf=$kf A^-1 generated by makenn\n");}
elsif ($rkky==2)
{print $l ("# kaneyoshi parametrization for the Bethe-Slater curve J(R)= A [-(R/D)^2+(R/D)^4].exp[-alpha.(R/D)^2] for scale A=$scale meV kf=$D A^-1  alpha=$aa\n");}
elsif ($rkky==4)
{print $l ("# it follows output of RKKY interaction J(R)=A.cos(2.kfR)/(2.kfR)^3 for scale A=$scale meV and\n");
 print $l ("# kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2) with ka=$ka A^-1 kb=$kb A^-1 kc=$kc A^-1\n");}
elsif($rkky==5)
{print $l ("# kaneyoshi parametrization for the Bethe-Slater curve J(R)= A [-(RD)^2+(RD)^4].exp[-alpha.(RD)^2] for scale A=$scale meV\n");
 print $l ("# with RD=sqrt(Ra^2/Da^2+Rb^2/Db^2+Rc^2/Dc^2) with Da=$Da A^-1 Db=$Db A^-1 Dc=$Dc A^-1  and alpha=$aa\n");}
elsif($rkky==6)
 {print $l ("# it follows output of RKKY interaction  J(R)=A [sin(2.kfR)-2.kfR.cos(2.kfR)]/(2.kfR)^4 for scale A=$scale meV\n");
  print $l ("# kfR=sqrt(ka^2.Ra^2+kb^2.Rb^2+kc^2.Rc^2) with ka=$ka A^-1 kb=$kb A^-1 kc=$kc A^-1\n");}
else
{print $l ("# it follows output of classical DD interaction generated by makenn\n");}

 for ($n1=1;$n1<(($rn->dims)[0]);++$n1)
# the position xyz is relative position (not absolute coordinate of neighbour)
 {print $l sprintf("%4.4g %4.4g %4.4g ",$xn->index($n)->at($n1),$yn->index($n)->at($n1),$zn->index($n)->at($n1));
# {print $l ($xn->index($n)->at($n1)." ".$yn->index($n)->at($n1)." ".$zn->index($n)->at($n1))." ";
 print $l ($Jaa->index($n)->at($n1)." ".$Jbb->index($n)->at($n1)." ".$Jcc->index($n)->at($n1))." ";
 if ($rkky==0)
  { print $l ($Jab->index($n)->at($n1)." ".$Jba->index($n)->at($n1)." ".$Jac->index($n)->at($n1))." ";
    print $l ($Jca->index($n)->at($n1)." ".$Jbc->index($n)->at($n1)." ".$Jcb->index($n)->at($n1))." ";
  }
 $ddd=$an->index($n)->at($n1);
 if ($calcdist==1) {print $l ($rn->index($n)->at($n1)." a".$ddd);}
  print $l "\n";
 }
}

sub endprint {
  my ($h,$l)=@_;   
     close $h; 
     close $l;
}
