#! /usr/bin/perl

# perl packages to use
use FileHandle;
use PDL;
use PDL::Graphics::PGPLOT;
use PGPLOT;

if ($ARGV[0]=~m|-h.*|)
{print "program view: creates display and postscript plot of data file
		      
    use as view colx coly [colsymbsize] [colxerrorbars] [colyerrorbars] filename
          colx,coly .....  x- and y-axis
          colsymbolsize... column containing symbol size
          colyerrorbars ... column containing y-errorbars
          colxerrorbars ... column containing y-errorbars
	  filename... filename of data file
    
\n";

exit(0);
}

#get filename from command line if given
my (@file) = @ARGV;
my (@file1) = @ARGV;
$colx=$ARGV[0];shift @ARGV;
$coly=$ARGV[0];shift @ARGV;
if(($in)=($ARGV[0]=~m|\A(\d)\Z|)){$colsymbsize=$ARGV[0];shift @ARGV;}
if(($in)=($ARGV[0]=~m|\A(\d)\Z|)){$colxerror=$ARGV[0];shift @ARGV;}
if(($in)=($ARGV[0]=~m|\A(\d)\Z|)){$colyerror=$ARGV[0];shift @ARGV;}

    overview_plot("/xserv", @file); # plot on screen
    overview_plot("./dd.ps/vcps", @file1); #then plot on psfile
    system ("awk '{gsub(\"1 setlinewidth\",\" \"); print \$0}' ./dd.ps > ./view.ps\n");
    print "Wrote postscript file 'view.ps'\n";


sub overview_plot {
    my ($devspec, @file) = @_;
#open page on pgplot server to be able to display/print
    my $dev = dev("$devspec");
    die "PGOPEN failed!" unless $dev > 0;
#plot data
    pgpap(6.5,0.8); # set paper size
    pgslw 6; # set linewidth
    pgsch(2.5); # set character height
    pgvstd; # set standard (default) viewport (must come after pgsch !!!)
    hold;   # wait with display until release command is given

    my ($data,$Limits) = get_detector_data_2D($ARGV[0]);
    my ($x0,$x1,$y0,$y1,$xcol,$ycol,$xtext,$ytext) = @{$Limits};
    pgswin ($x0,$x1,$y0,$y1); #set window
#    pgwnad -100,100,-100,100;
    pgpage;
    pgbox("BCNSTI", 0, 0, "BCNTI", 0, 0); #draw labeled frame arount viewport
    #pgtbox("ZYHBCNST",0,0,"ZYDBCNST",0,0);
#pgline (2,0.6667,0.6667,0,2);
    
    #label plot
    pglab($xtext,$ytext, '\frNdCu\d2\u T=2 K');
#    pgmtxt("T", 2.5, 0.5, 0.5, "Excitation Energy ".$files);

#      pgarro(0.4, 0.3, 0.4, 0.1);
#      pgptxt(0.4, 0.3, 0.0, 0.0, '(0.4 1 0)');

#      pgarro(0.3333, 0.25, 0.3333, 0.05);  # plot an arrow and write info
#      pgptxt(0.2, 0.3, 0.0, 0.0, '(1/3 1 0)');
#      pgptxt(0.9, -0.15, 0.0, 0.0, '(1 1 0)');
      
      #    shift @file;
#      pgptxt(0.98, 2.1, 0.0, 0.0, '\fr\gG'); # Z B Gamma
#      pgptxt(0.32, 2.1, 0.0, 0.0, '\fr\gG');
#      pgptxt(-0.03, 2.1, 0.0, 0.0, '\frZB');
#      pgptxt(0.65, 2.1, 0.0, 0.0, '\frZB');

#      pgptxt(0.7, 1.8, 0.0, 0.0, '\frT=2 K');  #some info
#      pgptxt(0.7, 0.1, 0.0, 0.0, '\fr\gm\d0\uH||b=1T');
    $col=0;$files="";
    while($file[0])
    {
    $colx=$file[0];shift @file;
    $coly=$file[0];shift @file;
    if(($in)=($file[0]=~m|\A(\d)\Z|)){$colsymbsize=$file[0];shift @file;}
    if(($in)=($file[0]=~m|\A(\d)\Z|)){$colxerror=$file[0];shift @file;}
    if(($in)=($file[0]=~m|\A(\d)\Z|)){$colyerror=$file[0];shift @file;}
    ++$col;if ($col==1){pgsci(3);}else{pgsci(1);}
    $files=$files." ".$file[0];
    print $colx." ".$coly." ".$file[0]."\n";
     $xcol=$colx-1;
     $ycol=$coly-1;
     $intcol=$colsymbsize-1;

    my ($data,$Limits) = get_detector_data_2D($file[0]);
    shift @file;
    $lw=($data->slice("$intcol,:"))**0.5;
    $maxint=max($lw);$minint=min($lw);
    if ($maxint==$minint) {$lw*=5;}# nop int column
    else {$lw/=$maxint/30;} # if there was a intcolumn
    $_=$devspec;
    $symb=$lw/30*3.9;
    $ch=$lw/20;
    $lwerr=$lw/5;
    if (/.*cps/&&$col>1){$lw*=0.2;$ch*=.6;$lwerr*=1;}
    
    
    #draw datapoints from slice of piddle "$data"
    if ($col==2){$symbol=new PDL(227,846,624,524);}
    if ($col==3){$symbol=new PDL(902,903,904,905);}
    for ($k=1;$k<(($data->dims)[1]);++$k)
     {pgslw $lw->at(0,$k)+1.0;
      pgsch $ch->at(0,$k)+1.0;
#      if ($col>1){
#      points $data->at($xcol,$k),$data->at($ycol,$k),$symbol->at(int($symb->at(0,$k)));}
#      else
      {points $data->at($xcol,$k),$data->at($ycol,$k),1;}
      $fwhmx=$data->at($colxerror-1,$k);
      $fwhmy=$data->at($colyerror-1,$k);
      pgslw $lwerr->at(0,$k)+1.5;
      if ($colxerror)
      {pgerrx 1,$data->at($xcol,$k)-$fwhmx/2,$data->at($xcol,$k)+$fwhmx/2,$data->at($ycol,$k),1.5;}
      if ($colyerror)
      {pgerry 1,$data->at($xcol,$k),$data->at($ycol,$k)-$fwhmy/2,$data->at($ycol,$k)+$fwhmy/2,1.5;}
      if ($k>1){
            
                 pgslw 20; #some vertical lines
                 pgsch 3.0; # set character height
                 pgsls 1; # set linestyle
                my $s1=new PDL($data->at($xcol,$k-1),$data->at($xcol,$k));
		my $s2=new PDL($data->at($ycol,$k-1),$data->at($ycol,$k));
		poly $s1,$s2;
                } 
     }
    }
#    pgslw 3; #some vertical lines
#    pgsch 1.0; # set character height
#    $s1=new PDL(0.333333,0.33333);$s2=new PDL(0.5,2);poly $s1,$s2;
#    pgsls 2; # set linestyle
#    $s1=new PDL(0.666667,0.666667);$s2=new PDL(0,2);poly $s1,$s2;

    pgebuf; #end the buffer and
    release; #display it

#close page on pgplot server
#	pgclos;
	return ();
}


# Get numerical data, reading it from file 
sub get_detector_data_2D {
    my ($info) = @_;
    # unless ($info->{'Numeric Data'})
    {($info->{'Numeric Data'}) = read_data_file_2D($info);
    # unless($in) {# let user select axes
    #              print "\n 1. h\n 2. k\n 3. l\n";
    #              print "select x axis (1-3) ? ";
    #              $in=<STDIN>;
    #             }
     $xcol=$colx-1;
     $ycol=$coly-1;
     $intcol=$colsymbsize-1;
     if ($in==1) {$xtext='\fr(h10)';}
     if ($in==2) {$xtext='\frk';}
     if ($in==3) {$xtext='\frl';}
     $xtext='\fr \gm\d0\uH\d||b\u(T)';
     $ytext='\fr M\d||b\u(\gm\dB\u/f.u.)';
     
     # calculate limits
     $x0=min(($info->{'Numeric Data'})->slice("$xcol,:"));         
     $x1=max(($info->{'Numeric Data'})->slice("$xcol,:"));         
     $y0=min(($info->{'Numeric Data'})->slice("$ycol,:"));         
     $y1=max(($info->{'Numeric Data'})->slice("$ycol,:"));         
     ($info->{'Limits'})=[$x0,$x1,$y0,$y1,$xcol,$ycol,$xtext,$ytext];
    
     } 
    return ($info->{'Numeric Data'},$info->{'Limits'});
}


# Read 2D numeric data, skipping comment lines.
sub read_data_file_2D {
    my ($file) = @_;
    my $h = new FileHandle;
    my @v=();
     my @xlist = ();
  # input data int piddle
  print $file."\n";
  if(open($h,$file))
  {  while(<$h>)
     {next if /^\s*#/;
      # detect > to see where int data starts [h k l en1 en2 en3 ... > int1 int2 int3 ...
       $x=new PDL(split " ",$_);
       push(@xlist,$x);
     }
     
     close $h; 
     return (cat @xlist);
    } else {
	print STDOUT "Warning: failed to read data file \"$file\"\n";
	return undef;
    }
}

