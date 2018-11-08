# this is an inputfile for gnuplot used to plot results of 
# 1. free energy
# 2. magnetization
# 3. magnetostriction

# use as: gnuplot dlfemu.gnu

# NOTE: this file has probably to be modified for
#  	specific plots you are interested in

set terminal X11

fromline=1
toline=160

#set xrange [0:3]
#set yrange [-0.3:-0.0]
set xlabel  "H||b[T]"
set ylabel "energy u[meV]"
set title "Calculated Magnetic Energy T=2.0K"
plot "./mcphas.fum" every ::fromline::toline using 2:9 with lines

pause -1
set terminal postscript
set output "./uT2.ps"
replot

#*********************************************************************

set terminal X11

#set xrange [0:3]
#set yrange [-1.0:0.0]
#set dgrid3d 50,50,1
#set contour
set xlabel  "H||b[T]"
set ylabel "free energy [meV]"
set title "Calculated Free energy at T=2.0K"
plot "./mcphas.fum" every ::fromline::toline using 2:8 with lines

pause -1
set terminal postscript
set output "./feT0_5.ps"
replot

#*********************************************************************
set terminal X11

#set xrange [0:13]
#set yrange [0:2.5]
set xlabel  "H||b[T]"
set ylabel "magnetic moment [mub/f.u.]"
set title "Calculated Magnetic Moment T=2.0K"
plot "./mcphas.fum" every ::fromline::toline using 2:10
#with lines

pause -1
set terminal postscript
set output "./momT0_5.ps"
replot



#*********************************************************************
set terminal X11

set xrange [0:13]
set yrange [-5:10]
set xlabel  "H||c[T]"
set ylabel "JbJb(2),JcJc(2)"
set title "Calculated Magnetostriction T=0.5K"
plot "./mcphas.j2" every ::fromline::toline using 2:7 \
,"./mcphas.j2" every ::fromline::toline using 2:8 
#with lines

pause -1
set terminal postscript
set output "./jj2T0_5.ps"
replot
