#! /usr/bin/perl
use FileHandle;
use PDL;
use PDL::Slatec;

unless ($#ARGV>=0) 
{print " program to create table with neighbors from mcphas.j \n\n";
print " usage: makenn 23.3 \n\n";
print " meaning take mcphase.j, generate all neighbors within sphere of 23.3A and\n";
print " put them into makenn.j, in interaction columns put the classical dipole interaction\n";
print " the output values are sorted by ascending distance\n";
print " perl packages to use\n";
 exit 0;}

my ($rmax) = $ARGV[0];
my ($latt,$p) = getlattice("./mcphas.j");
my ($a,$b,$c,$gJ) = @{$latt};
     print "rmax=".$rmax." A\n";
     print "a=".$a." b=".$b." c=".$c."\n";
     print "primitive lattice[abc]:".$p."\n";
     $t=$p->slice("0,"); $t*=$a;
     $t=$p->slice("1,"); $t*=$b;
     $t=$p->slice("2,"); $t*=$c;
     print "primitive lattice[A]:".$p."\n";
     my ($rn)=new PDL ();

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

    $r=0;
# determine $nmin,$nmax by looking at a cube with side 2rmax
     $inv=matinv($p);
  for ($i1=-1;$i1<=1;$i1+=2){
  for ($i2=-1;$i2<=1;$i2+=2){
  for ($i3=-1;$i3<=1;$i3+=2){
    $n=inner($inv , pdl[$i1*$rmax,$i2*$rmax,$i3*$rmax]); 
    if (($n->at(0))<$n1min) {$n1min=int($n->at(0))-1;}
    if (($n->at(1))<$n2min) {$n2min=int($n->at(1))-1;}
    if (($n->at(2))<$n3min) {$n3min=int($n->at(2))-1;}
    if (($n->at(0))>$n1max) {$n1max=int($n->at(0))+1;}
    if (($n->at(1))>$n2max) {$n2max=int($n->at(1))+1;}
    if (($n->at(2))>$n3max) {$n3max=int($n->at(2))+1;}
  }}}
print "$n1min to $n1max, $n2min to $n2max, $n3min to $n3max\n";

  for ($n1=$n1min;$n1<=$n1max;++$n1){ 
  for ($n2=$n2min;$n2<=$n2max;++$n2){ 
  for ($n3=$n3min;$n3<=$n3max;++$n3){  
  
   $rvec=$n1*$p->slice(",(0)")+$n2*$p->slice(",(1)")+$n3*$p->slice(",(2)");
   $rr=inner($rvec, $rvec);
   $r=sqrt($rr->at());
   $x=$rvec->at(0)/$a;
   $y=$rvec->at(1)/$b;
   $z=$rvec->at(2)/$c;
#       print $r."   ($x,$y,$z,$n1,$n2,$n3\n";
#   if ($r<=$rmax && $r>0 && $x>=0 && $y>=0 && $z>=0){#save neighbour jjj format
   if ($r<=$rmax && $r>0){#save neighbour j format
    $rn=$rn->append( pdl ([$r]));
    $xn=$xn->append( pdl ([$x]));
    $yn=$yn->append( pdl ([$y]));
    $zn=$zn->append( pdl ([$z]));
    my ($interaction) = getinteraction($gJ,$r,$x*$a,$y*$b,$z*$c);
    my ($jaa,$jab,$jac,$jba,$jbb,$jbc,$jca,$jcb,$jcc) = @{$interaction};

    $Jaa=$Jaa->append( pdl ([$jaa]));
    $Jab=$Jab->append( pdl ([$jab]));
    $Jac=$Jac->append( pdl ([$jac]));
    $Jba=$Jba->append( pdl ([$jba]));
    $Jbb=$Jbb->append( pdl ([$jbb]));
    $Jbc=$Jbc->append( pdl ([$jbc]));
    $Jca=$Jca->append( pdl ([$jca]));
    $Jcb=$Jcb->append( pdl ([$jcb]));
    $Jcc=$Jcc->append( pdl ([$jcc]));

   }
   }}}  
   $n= qsorti($rn);
   $nofneighbours+=(($rn->dims)[0]);
#   printlattice("./mcphas.jjj",$n,$rn,$xn,$yn,$zn,$Jaa,$Jbb,$Jcc,$Jab,$Jba,$Jac,$Jca,$Jbc,$Jcb,">./makenn.jjj");
   printlattice("./mcphas.jjj",$n,$rn,$xn,$yn,$zn,$Jaa,$Jbb,$Jcc,$Jab,$Jba,$Jac,$Jca,$Jbc,$Jcb,">./makenn.j");
   exit;
#-----------------------------------------------------------------------

sub getinteraction {
   my ($gJ,$r,$rx,$ry,$rz)=@_;
   # calculate classical dipole interaction
my $c = $gJ * $gJ * .927405 * .927405 / 16.02183;  #[meV A^3]

$jaa = $c * (3 * $rx * $rx -$r *$r) /$r /$r /$r /$r /$r
$jbb = $c * (3 * $ry * $ry -$r *$r) /$r /$r /$r /$r /$r
$jcc = $c * (3 * $rz * $rz -$r *$r) /$r /$r /$r /$r /$r


$jab = $c * (3 * $rx * $ry) /$r /$r /$r /$r /$r
$jbc = $c * (3 * $ry * $rz) /$r /$r /$r /$r /$r
$jac = $c * (3 * $rx * $rz) /$r /$r /$r /$r /$r
   $jba=$jab;
   $jcb=$jbc;
   $jca=$jac;
      
 return ([$jaa,$jab,$jac,$jba,$jbb,$jbc,$jca,$jcb,$jcc]);  
}


# Get lattic data, reading it from file 
sub getlattice {
    my ($file) = @_;
    my $h = new FileHandle;
#     my @xlist = ();
  # input data int piddle
  if(open($h,$file))
  {  while(<$h>)
     {next if /^\s*#/;
      # detect a= b= c= ...
      if (/^\s*\Qa=\E/){($a)=(/^\s*\Qa=\E\s*([^\s]*)\s/);}
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
      if (/^.*\Qnofneighbours=\E/){($nofneighbours)=(/^.*\Qnofneighbours=\E\s*([^\s]*)\s/);}
      if (/^.*\QgJ=\E/){($gJ)=(/^.*\QgJ=\E\s*([^\s]*)\s/);}
     }
     close $h; 
     return ([$a,$b,$c,$gJ],pdl [[$r1x,$r1y,$r1z],[$r2x,$r2y,$r2z],[$r3x,$r3y,$r3z]]);
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
     {next if /^\s*#/;
      $text=$_;
      if (/^.*\Qnofneighbours=\E/){$text=~s!\Qnofneighbours=\E\s*\d+!nofneighbours=$nofneighbours!;}
      print $l ($text);      
     }
     close $h; 
     print $l ("# it follows output of classical DD interaction generated by makenn\n");
 for ($n1=0;$n1<(($rn->dims)[0]);++$n1)
 {print $l ($xn->index($n)->at($n1)." ".$yn->index($n)->at($n1)." ".$zn->index($n)->at($n1))." ";
 print $l ($Jaa->index($n)->at($n1)." ".$Jbb->index($n)->at($n1)." ".$Jcc->index($n)->at($n1))." ";
 print $l ($Jab->index($n)->at($n1)." ".$Jba->index($n)->at($n1)." ".$Jac->index($n)->at($n1))." ";
 print $l ($Jca->index($n)->at($n1)." ".$Jbc->index($n)->at($n1)." ".$Jcb->index($n)->at($n1))."\n";
  }
 close $l;
}