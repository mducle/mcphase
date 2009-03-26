#!/usr/bin/perl
use FileHandle;
use PDL;
use PDL::Slatec;

unless ($#ARGV>=2) 
{print " program to extend crystallographic unit cell n times in r1 (or r2,r3) direction\n\n";
print " usage: extendunitcell 3 1 4\n\n";
print " meaning take mcphas.j, mcphas.tst and mcdiff.in and generate an extended description of the unit cell 3xr1,1xr2,4xr3\n";
print " put result into results/extend.j, results/extend.tst and results/extend.in\n";
 exit 0;}

my ($n1) = $ARGV[0];shift @ARGV; 
my ($n2) = $ARGV[0];shift @ARGV; 
my ($n3) = $ARGV[0]; 

print "reading mcphas.j ....\n";
 my ($latt,$p) = getlattice("./mcphas.j");
 my ($a,$b,$c,$nofa) = @{$latt};
 print "a=".$a." b=".$b." c=".$c."\n";
 print "original primitive lattice[abc]:".$p."\n";


 #initialize output file extendj.j
 my ($h,$l)=printlattice("./mcphas.j",$n1,$n2,$n3,">./results/extend.j");
  printneighbourlist($h,$l,$n1,$n2,$n3,$p,$nofa);   
  endprint($h,$l);   
 # extend lattice
 $p->slice(",0")*=$n1;
 $p->slice(",1")*=$n2;
 $p->slice(",2")*=$n3;

 print "new primitive lattice[abc]:".$p."\n";

print "reading mcphas.tst ...\n";

 open ($h,"mcphas.tst");
 open ($l,">./results/extend.tst"); 
 while(<$h>)
 {if (/^\s*#/)
  {$text=$_;
   if (/^.*\Qnofcomponents=\E/){($nofcomponents)=(/^.*\Qnofcomponents=\E\s*([^\s]*)\s/);}
   $text=~s/\Qnofatoms=\E\s*[\d.]+/nofatoms=$nofatoms /;
   print $l ($text);
  }
  else
  { # read configuration
   my @list =();
   push(@list,new PDL(split " "));
   @dims=$list[0]->dims; $nr1=$dims[0];
    for($i=2;$i<=$nofcomponents;++$i)
     {$_=<$h>;push(@list,new PDL(split " "));}
   $nr2=1;$nr3=1;
   while(<$h>)
   {last if !/^\s*\d/;
    push(@list,new PDL(split " "));
    for($i=2;$i<=$nofcomponents;++$i)
     {$_=<$h>;push(@list,new PDL(split " "));}
    $nr2++;
   } 
   if (!/^\s*#/)
   {while(<$h>)
    {last if /^\s*#/;
     $nr3++;
     for ($j=1;$j<=$nr2;++$j)
     {    for($i=1;$i<=$nofcomponents;++$i)
          {push(@list,new PDL(split " "));$_=<$h>;}
     }
     last if /^\s*#/;
    }
   }
   $conf=new PDL(cat @list);
   # extend configuration
   # write configuration
   for($ii3=1;$ii3<=$n3;++$ii3){for($i3=1;$i3<=$nr3;++$i3){
   for($ii2=1;$ii2<=$n2;++$ii2){for($i2=1;$i2<=$nr2;++$i2){
    for($i=1;$i<=$nofcomponents;++$i){
     for($ii1=1;$ii1<=$n1;++$ii1){for($i1=1;$i1<=$nr1;++$i1){
     $j=$i1;$offset=$nofcomponents*(($i3-1)*$nr2+$i2-1);
     print $l $conf->at($j-1,$offset+$i-1)." ";
     }} print $l "\n";
   }}}if($i3<$nr3||$ii3<$n3) {print $l "\n";}}}
   print $l $_;
  }

 } 
 close $h,$l;


print "reading mcdiff.in ...\n";
 open ($h,"mcdiff.in");
 open ($l,">./results/extend.head"); 
 while(<$h>)
 {$text=$_;
  if (/^\s*#/)
  {
   if (/^.*\Qnat=\E/){($nat)=(/^.*\Qnat=\E\s*([^\s]*)\s/);
                      $nofatoms=$nat*$n1*$n2*$n3;$text=~s/\Qnat=\E\s*[\d.]+/nat=$nofatoms /;}
   print $l ($text);
  }
  else
  {@numbers=split(" ");
   $x=$numbers[2];
   $y=$numbers[3];
   $z=$numbers[4];
   $dr1=$numbers[5];
   $dr2=$numbers[6];
   $dr3=$numbers[7];

   for ($i1=0;$i1<=$n1-1;++$i1){
   for ($i2=0;$i2<=$n2-1;++$i2){
   for ($i3=0;$i3<=$n3-1;++$i3){
     $da=$x+$i1*$p->at(0,0)+$i2*$p->at(0,1)+$i3*$p->at(0,2);
     $db=$y+$i1*$p->at(1,0)+$i2*$p->at(1,1)+$i3*$p->at(1,2);
     $dc=$z+$i1*$p->at(2,0)+$i2*$p->at(2,1)+$i3*$p->at(2,2);
     $dr1n=($dr1+$i1)/$n1;
     $dr2n=($dr2+$i2)/$n2;
     $dr3n=($dr3+$i3)/$n3;
   print $l $numbers[0]." ".$numbers[1]." ".$da." ".$db." ".$dc." ".$dr1n." ".$dr2n." ".$dr3n." ".$numbers[8]."\n";
   }}} 
  }
  last if /.*SECTION\s*3/;
 }
 
 close $h,$l;
 $exit;
  

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
  my ($filein,$n1,$n2,$n3,$fileout)=@_;   
#   print $rn->index($n)."\n";
    my $h = new FileHandle;
    my $l = new FileHandle;
    open($l,$fileout);
    open($h,$filein);
     while(<$h>)
     {#next if /^\s*#/;
      $text=$_;
      if (/^.*\Qr1x=\E/){($r1x)=(/^.*\Qr1x=\E\s*([^\s]*)\s/);$r1x*=$n1;$text=~s/\Qr1x=\E\s*[\d.]+/r1x=$r1x /;}
      if (/^.*\Qr1y=\E/){($r1y)=(/^.*\Qr1y=\E\s*([^\s]*)\s/);$r1y*=$n1;$text=~s/\Qr1y=\E\s*[\d.]+/r1y=$r1y /;}
      if (/^.*\Qr1z=\E/){($r1z)=(/^.*\Qr1z=\E\s*([^\s]*)\s/);$r1z*=$n1;$text=~s/\Qr1z=\E\s*[\d.]+/r1z=$r1z /;}
      if (/^.*\Qr2x=\E/){($r2x)=(/^.*\Qr2x=\E\s*([^\s]*)\s/);$r2x*=$n2;$text=~s/\Qr2x=\E\s*[\d.]+/r2x=$r2x /;}
      if (/^.*\Qr2y=\E/){($r2y)=(/^.*\Qr2y=\E\s*([^\s]*)\s/);$r2y*=$n2;$text=~s/\Qr2y=\E\s*[\d.]+/r2y=$r2y /;}
      if (/^.*\Qr2z=\E/){($r2z)=(/^.*\Qr2z=\E\s*([^\s]*)\s/);$r2z*=$n2;$text=~s/\Qr2z=\E\s*[\d.]+/r2z=$r2z /;}
      if (/^.*\Qr3x=\E/){($r3x)=(/^.*\Qr3x=\E\s*([^\s]*)\s/);$r3x*=$n3;$text=~s/\Qr3x=\E\s*[\d.]+/r3x=$r3x /;}
      if (/^.*\Qr3y=\E/){($r3y)=(/^.*\Qr3y=\E\s*([^\s]*)\s/);$r3y*=$n3;$text=~s/\Qr3y=\E\s*[\d.]+/r3y=$r3y /;}
      if (/^.*\Qr3z=\E/){($r3z)=(/^.*\Qr3z=\E\s*([^\s]*)\s/);$r3z*=$n3;$text=~s/\Qr3z=\E\s*[\d.]+/r3z=$r3z /;}
      if (/^.*\Qr1a=\E/){($r1a)=(/^.*\Qr1a=\E\s*([^\s]*)\s/);$r1a*=$n1;$text=~s/\Qr1a=\E\s*[\d.]+/r1a=$r1a /;}
      if (/^.*\Qr1b=\E/){($r1b)=(/^.*\Qr1b=\E\s*([^\s]*)\s/);$r1b*=$n1;$text=~s/\Qr1b=\E\s*[\d.]+/r1b=$r1b /;}
      if (/^.*\Qr1c=\E/){($r1c)=(/^.*\Qr1c=\E\s*([^\s]*)\s/);$r1c*=$n1;$text=~s/\Qr1c=\E\s*[\d.]+/r1c=$r1c /;}
      if (/^.*\Qr2a=\E/){($r2a)=(/^.*\Qr2a=\E\s*([^\s]*)\s/);$r2a*=$n2;$text=~s/\Qr2a=\E\s*[\d.]+/r2a=$r2a /;}
      if (/^.*\Qr2b=\E/){($r2b)=(/^.*\Qr2b=\E\s*([^\s]*)\s/);$r2b*=$n2;$text=~s/\Qr2b=\E\s*[\d.]+/r2b=$r2b /;}
      if (/^.*\Qr2c=\E/){($r2c)=(/^.*\Qr2c=\E\s*([^\s]*)\s/);$r2c*=$n2;$text=~s/\Qr2c=\E\s*[\d.]+/r2c=$r2c /;}
      if (/^.*\Qr3a=\E/){($r3a)=(/^.*\Qr3a=\E\s*([^\s]*)\s/);$r3a*=$n3;$text=~s/\Qr3a=\E\s*[\d.]+/r3a=$r3a /;}
      if (/^.*\Qr3b=\E/){($r3b)=(/^.*\Qr3b=\E\s*([^\s]*)\s/);$r3b*=$n3;$text=~s/\Qr3b=\E\s*[\d.]+/r3b=$r3b /;}
      if (/^.*\Qr3c=\E/){($r3c)=(/^.*\Qr3c=\E\s*([^\s]*)\s/);$r3c*=$n3;$text=~s/\Qr3c=\E\s*[\d.]+/r3c=$r3c /;}
      if (/^.*\Qnofatoms=\E/){($nofatoms)=(/^.*\Qnofatoms=\E\s*([^\s]*)\s/);$nofatoms*=$n1*$n2*$n3;$text=~s/\Qnofatoms=\E\s*[\d.]+/nofatoms=$nofatoms /;}

      last if /^#.*\Q**********\E/;
      print $l ($text);      
     }
 return ($h,$l);
}

sub printneighbourlist {
  my ($h,$l,$n1,$n2,$n3,$p,$nofatoms)=@_;   

  for($n=1;$n<=$nofatoms;++$n)
  { 

   # get a neighbour
   $nn=0;
     while(<$h>)
     {last if /^#.*\Q**********\E/;
      ++$nn;
      $text[$nn]=$_;
      }
    
   for ($i1=0;$i1<=$n1-1;++$i1){
   for ($i2=0;$i2<=$n2-1;++$i2){
   for ($i3=0;$i3<=$n3-1;++$i3){
     $da=$x[$n]+$i1*$p->at(0,0)+$i2*$p->at(0,1)+$i3*$p->at(0,2);
     $db=$y[$n]+$i1*$p->at(1,0)+$i2*$p->at(1,1)+$i3*$p->at(1,2);
     $dc=$z[$n]+$i1*$p->at(2,0)+$i2*$p->at(2,1)+$i3*$p->at(2,2);
     # print out neighbor
     print $l ("#*************************************************************************\n");
     for($i=1;$i<=$nn;++$i)
       {
        $text[$i]=~s!\Qda=\E\s*[\d.]+!da=$da!;
        $text[$i]=~s!\Qdb=\E\s*[\d.]+!db=$db!;
        $text[$i]=~s!\Qdc=\E\s*[\d.]+!dc=$dc!;
        print $l ($text[$i]) 
       }
   }}}

  } 
}

sub endprint {
  my ($h,$l)=@_;   
     close $h; 
     close $l;
}
