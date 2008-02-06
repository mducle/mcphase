\begin{verbatim}

For dos/windows operating systems the following programs can be used
in addition.

Note the required file format in the following example:

----------------------------------------------------------------
{ this is a test file the fileheader is enclosed in brackets
then there follow some columns (up to 30) "1-3" or "a-c" of data}
1 2.3 3.145
3 2.5 3.1451
4 2.7 3.1457
5 3   3.1458
----------------------------------------------------------------


In order to manipulate such files you have the following commands available:

--------------------------------------
1. FOR GRAPHICAL DISPLAY AND PRINTING:
--------------------------------------

VIEW *.* 13         - means view columns 1 vs.3, add option /c if you
                      want to continue without typing anything after
                      the graph has been viewed, you can view also more
                      than one file by typing view *.* ij *.* ij ...)
                      after running view there is created a batch viewprn.bat
                      which can be started to create automatically
                      a postscript file of the viewed data (it uses
                      the programs PRINTF,GRAPH and PLOT described below)
                      points can be deleted - then the program view
                      creates a file delpoint.del, which contains the
                      new data set.
		      
GRAPH {xaxis=axistext] {yaxis=axistext]
                      - makes the scales and the axis TITLES (axistext ...)
                        in the postscript file graph.ps ... remember only
                        when the option /s is added the postscript file is
                        finished using the 'showpage' command          
                        you MUST have defined dos variables minx maxx miny
                        maxy (range of graph) and locx locy heigth width
                        (location and size of the GRAPH (in mm))

PLOT *.* 13 {size=7]  {xbars=4] {ybars=5] {symbol=a(or 0 for empty circle
                                      . for full circel or - for line)] {/s]
                     - the example plots data stored in file *.* xaxis
                     is column1 yaxis is column 3 symbols are 'a' (or circle
                     or lines) with size column7 using column 4 for the
                     x-errorbars  and col5 for y-errorbars to the postscript                                                    
                     file plot.ps ... remember only when the option /s is
                     added the postscript file is finished using the
                     'showpage' command you must have defined dos variables
                     minx maxx miny ymaxy (range of plot) and locx locy
                     heigth width (location and size of the plot (in mm))
                     use program graph to create the scales for such a plot                                                      "
                     the filename of the text is just written, where the
                     cursor is located on the postskript page
		     
DEOVERLP *.* 13 {size=7]{xbars=4]{xlow=2 xup=4]{ybars=5]{ylow=2 yup=3]{/tx/ty/dx=6/dy=8/b]
                 ->DEOVERLaPs data stored in file *.*, xaxis is column1, yaxis is column 3"
                files DEOVERLP.n are created such that the data points do not overlap in plot"
                size: size of symbol(<1)"
                xbars: x-error bar column   ybars: y-errorbars column"
                xlow/xup: lower and upper measurement limits in x direction"
                ylow/yup: lower and upper measurement limits in y direction"
                option /dx=6 /dy=8 displays numbers on plot corresponding to col6 and 8"
                option /tx, /ty displays comment text for each data point (horizontal,vertical)"
                option /b allows overlap within data blocks separated by {multiline comments}"

PRINTPS *1.* *2.* ... - sends some postskript files to the printer

PRINTF *1.* 23 *2.* 32 ...
                      - uses programs GRAPH and PLOT to make a postskript
                        file of a standard figure containig the data of files
                        *1.* column 2 vs col.3 ...
                        MIND: before startin this program the dos variables
                        minx maxx miny maxy must be set according to the
                        axis ranges desired ...

------------------------
2.  COLUMN MANIPULATIONS:
------------------------

FFORM *.* 3 ##.##   - formats 3rd column of file with basic number format
DELCOL *.* 25       - deletes columns 2 to 5 from file *.*
NEWCOL *.* 3        - inserts new column nr.3 in file *.* (line number col)
SWAP *.* 13         - swap column 1 and 3
SWAPF *.* *.*       - swaps the first columns of two files


------------------------------
3. LINE/COMMENT MANIPULATIONS:
------------------------------
  
CUTHEAD *.* 13      - cuts the fileheader in file *.* to 13 lines
SORTF *.* 2         - sort data points according to ascending column 2

REDUCE *.* 2 dmin=0.1 -- reduces points in data file such that minimal
                      difference between neighbouring data points in column
                      2 is bigger than dmin=0.1
AVERAGE *.* 4 [options] .... takes 4 points and averages data
                      option
                               middle point is taken
                      /av      points are averaged
                      /med     median of points is calculated and kept

---------------------------
4. MATHEMATICAL OPERATIONS:
---------------------------

DIFF *.* 23         - differentiate column 3 with respect to column 2, adding
                      option "stp=2.4" sets a minimum averaging interval(col2)
                      otherwise the differentiation is done point by point.
                      the differentiation is done by sampling the xaxis
                      (2nd column) in steps of 2.4(or point by point) and
                      caculating the linear regression with respect to such
                      an xaxis interval for all other columns !!
                      to the yaxis (3rd column in our example) the slope of
                      its linear regression is written, to all other columns
                      (z1,z2..axis) the value of the regression at <x> is
                      written (corresponds to <z1>,<z2>..)
                      NOTE: all columns are changed by this operation !!!

POT *.* 2 24.1      - col2=col2^24.1
EXP *.* 1           - col1=exp(col1)

INVERT *.* 3        - invert columns 2
FACT *.* 3 2.42     - multiply col. 2 with the constant 2.42
SHIFTC *.* 3 3.47   - add 3.47 to column 3 option /+ shifts more files with
                      increasing the shift in each file
ZSHIFT *.* 23 23.1  - add a number to col 3 such that function col3(col2) is
                      zero at 23.1 (used for instance to shift dl/l values)
RSHIFT *1.* 23 *2.* 14 24.4 - shifts the a function relativ to another function
                    - y(x)=col3(col2)  in file *1.* is shifted in y by adding a
                      constant, the constant is chosen such that the minimum
                      difference between y=col3 (in file *1.*) and y=col4(col1)
                      (in file *2.*) is 24.4         

MULT *1.* 13 *2.* 21  - multiply function col3(col1) of file *1 with
                        function col1(col2) of file *2 (using linear
                        interpolation to find corresponding values of col1 in
                        file *1 and col2 in file *2,NOTE: for linear interpol.
                        of file *2.* the data of  col2 must be sorted
                        according to increasing or decreasing values)
                         -> the result is written into col3 of file *1.*
ADD  *1.* 13 *2.* 21  - add function col3(col1) of file *1 to function
                        col1(col2) of file *2 (using linear interpolation
                        to find corresponding values of col1 in file *1 and
                        col2 in file *2, NOTE: for linear interpolation
                        of file *2.* the data of col2 must be sorted
                        according to increasing or decreasing values)
                         -> the result is written into col3 of file *1.*
                        options:  /noex   ... do not extrapolate
ADDC *.* 13           - column 3 of file *.* is added to column 1 point by
                        point - the result is written in column 1
MULTC *.* 13          - column 3 of file *.* is multiplied with col 1 point by
                        point - the result is written to column 1
SPLINE *.* 1  spl=2   - column 1(2) are taken as x(y)-axis for calculating
                        the curvature columns of a cubic spline function
                        this column is added to the other columns as the last
                        column
INTERPOL *.* 1 stp=0.2 spl=3
                      - the column 1 is taken as xaxis for interpolating
                        the other columns with a stepwidth of 0.2. the
                        option spl=3 means that the 3rd column will be
                        interpolated by a cubic spline function the curvature
                        of which is taken from the last column (which is del)
                        (option newton=5 means do newton interpolation using 5
                        nearest neighbouring points)
SPLINT *.* 1 stp=0.2 spl=3
                      - just combines spline and interpol commands by
                        spline interpolating y(x) as col3(col1) with
                        stepwidth .2




CONVOLUT *.* 23 stp=0.014 mode=gauss fwhm=0.2
                      - (means take file *.* second column as xaxis and
                        calculate CONVOLUTation for y axis=3rd column  with
                        stepwidth 0.014) --->the result is written in file *.cvt 
                        use the mode= switch to specify the type
                        of convolution function  used. 
                        convolution is done according to the formula:
                        conv(x')=sum_{x} col3(x)*f(x'-x)

                        possible modes are: gauss, fullprof2t, fullprofq
                        the function parameters of f(x) : 
                              for gauss: fwhm=, 
                              for fullprof2t: u= v= w=, 
                              for fullprofq:  u= v= w=,lamdba=[A]
                        [u,v,w used to calculate fwhm from 2theta: 
                              for D9-ILL u=16.08575, v=-6.2195, w=0.83903]
                        
                        example: this procedure can be used to calculate a diffraction
                        pattern from the output of program ELN

LINREG *.* 23         - 23 means calc.LINREG of 3rd columns with respect
                        to second column)---> the result is written in
                        file *.*, the LINear REGression is calculated by
                        least squares method, the linear regression is
                        written to the last column of file *.*    
SMOOTHE *.* 23 stp=4.46
                      - means SMOOTHE 3rd column with respect to second
                        column) ---> the result is written in file *.*,                                 "
                        smoothing is done in col2-intervals of 2.34 (+1
                        datapoint extra) the SMOOTHing  is done by sampling
                        the xaxis(2nd column) in steps of 2.34 and
                        calculating the linear regression at the average <x>
                        of such an xaxis interval, if a data point exceeds
                        the standard deviation of this linear regression then
                        its value is changed to the standard deviation but
                        only within a small interval around <x>
FOURIER *.* 12 xmax=3.47 deltax=0.1
                      - means fouriertransform functions coln(x)=(other columns)
                        x=col1 -> in the range x={-3.47;3.47] y={-10,10]
                        datapoint sampling is done in intervals deltax and
                        linear interpolation;the complex fouriertransform is
                        written into files FOURIER.re and FOURIER.im
                        format: {]..{]..xcolumn{1/{x]]..{]...


-----------------------------
5. NEUTRON POWDER DIFFRACTION:
-----------------------------

ELN *.*               - calculate neutron reflection and intensity list
                        with input file *.*  - such a file can easily be generated
                        from the output of the program spins (only some
                        structural data of nonmagnetic atoms and wavelength
                        and maximal angle have to be inserted)


here follows an example input file (see examples/ndcu2b_new/eln/spins.out):

--------------------------------------------------
# this file is the input file for program eln 
# eln is a program for the calculation of elastic
# neutron diffraction pattern
# 
# it contains 4 sections corresponding to different
# groups of parameters
#
# -all lines have to start with a # sign with the 
# exception of the lines containing atomic positional
# parameters
# -the other parameters have to be defined in the
# corresponding section by statements such as 
# parameter=value
# -the sequence of the parameters within a section is
# arbitrary
# -commas must not be used
# 
#
# %SECTION 1%  OVERALL PARAMETERS
#
#lambda=0.9 wavelength (A)
#thetamax=19 maximum bragg angle (deg)
#ovalltemp=0.3 overall temperature factor (A^2)
#lorentz=1   type of lorentzfactor to be used
#            1.....neutron powder flat sample
#            2.....neutron powder cylindrical sample
#            3.....neutron single crystal
#            4.....neutron TOF powder cyl. sample - d-pattern log scaled
#            5.....neutron TOF powder cyl. sample - d-pattern normal scaled
#
#
#
# %SECTION 2% LIST OF NONMAGNETIC ATOMS IN CRYSTALLOGRAPHIC UNIT CELL
#
#
#nat=4                  nr of nonmagnetic atoms in primitive crystalographic unit cell
#atom number x(a) y(b) z(c) dr1(r1) dr2(r2) dr3(r3)   lines with nonmagnetic atoms
# note: x y and z are not used by the program !
{Cu} 0 0.3021 0.7022  0 0.3021 0.7022
{Cu} 0 0.6979 0.7022  0 0.6979 0.7022
{Cu} 0 0.8021 0.3718  0 0.8021 0.3718
{Cu} 0 0.1979 0.3718  0 0.1979 0.3718
# ... note: if an occupancy less than 1.0 is needed - it can be inserted such as {0.6Cu}
#
#
#
# %SECTION 3% DESCRIPTION OF THE LATTICE
#
#
# ... note:what follows here may directly be taken from the output of program spins 
#          (file spins.out)
#
# lattice constants:                 lattice constants and primitive lattice vectors
# a=4.3843 b=7.0326 c=7.4194 alpha=  90 beta=  90 gamma=  90
# r1x= 0.5 r2x=   0 r3x=   0
# r1y= 0.5 r2y=   1 r3y=   0   primitive lattice vectors (a)(b)(c)
# r1z= 0.5 r2z=   0 r3z=   1
#
#
#
# %SECTION 4% DESCRIPTION OF MAGNETIC UNIT CELL AND LIST OF MAGNETIC ATOMS
#
#
# here follows the description of the magnetic unit cell with respect
# to the primitive crystallographic unit cell
# the crystallographic unit cell has to be taken nr1 nr2 and nr3 times
# along r1 r2 and r3 respectively to get magnetic unit cell
# nat denotes the number of magnetic atoms in magnetic unit cell
#
# nat lines follow to describe the magnetic moment configuration:
# 'atom number' means the type of atom (element symbol) followed
# by its coordinates with respect to abc (xyz) and to the primitive
# lattice (dr1 dr2 dr3) - note: x y and z are not used by the program !
# Col 8-10 contain the expectation values of the angular momentum components
# Ja Jb and Jc for this atom. 
#
#nr1=10 nr2=1 nr3=1 nat=20 atoms in primitive magnetic unit cell:
#atom number x(a) y(b) z(c) dr1(r1) dr2(r2) dr3(r3)  <Ja> <Jb> <Jc> ...
# note: x y and z are not used by the program !
{Nd3+} 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000  0.0000 2.8900 0.0000
{Nd3+} 0.0000 0.5000 0.0766 0.0000 0.5000 0.0766  0.0000 2.8900 0.0000
{Nd3+} 0.5000 0.5000 0.5000 1.0000 -0.0000 -0.0000  0.0000 2.9500 0.0000
{Nd3+} 0.5000 1.0000 0.5766 1.0000 0.5000 0.0766  0.0000 2.9500 0.0000
{Nd3+} 1.0000 1.0000 1.0000 2.0000 -0.0000 -0.0000  0.0000 -2.9500 0.0000
{Nd3+} 1.0000 1.5000 1.0766 2.0000 0.5000 0.0766  0.0000 -2.9500 0.0000
{Nd3+} 1.5000 1.5000 1.5000 3.0000 -0.0000 -0.0000  0.0000 -2.8900 0.0000
{Nd3+} 1.5000 2.0000 1.5766 3.0000 0.5000 0.0766  0.0000 -2.8900 0.0000
{Nd3+} 2.0000 2.0000 2.0000 4.0000 -0.0000 -0.0000  0.0000 2.9100 0.0000
{Nd3+} 2.0000 2.5000 2.0766 4.0000 0.5000 0.0766  0.0000 2.9100 0.0000
{Nd3+} 2.5000 2.5000 2.5000 5.0000 -0.0000 0.0000  0.0000 -2.8900 0.0000
{Nd3+} 2.5000 3.0000 2.5766 5.0000 0.5000 0.0766  0.0000 -2.8900 0.0000
{Nd3+} 3.0000 3.0000 3.0000 6.0000 -0.0000 -0.0000  0.0000 -2.9500 0.0000
{Nd3+} 3.0000 3.5000 3.0766 6.0000 0.5000 0.0766  0.0000 -2.9500 0.0000
{Nd3+} 3.5000 3.5000 3.5000 7.0000 -0.0000 0.0000  0.0000 2.9500 0.0000
{Nd3+} 3.5000 4.0000 3.5766 7.0000 0.5000 0.0766  0.0000 2.9500 0.0000
{Nd3+} 4.0000 4.0000 4.0000 8.0000 -0.0000 -0.0000  0.0000 2.8900 0.0000
{Nd3+} 4.0000 4.5000 4.0766 8.0000 0.5000 0.0766  0.0000 2.8900 0.0000
{Nd3+} 4.5000 4.5000 4.5000 9.0000 -0.0000 -0.0000  0.0000 -2.9100 0.0000
{Nd3+} 4.5000 5.0000 4.5766 9.0000 0.5000 0.0766  0.0000 -2.9100 0.0000
------------------------------------------
\end{verbatim}
After issuing the command {\prg eln spins.out} 
this input file gives the following output file:
\begin{verbatim}
--------------------------------------------------
#{.\eln.out 20:35:3007-18-2004   unit cell:
#r1x= 21.9215  A r2x= 0  A r3x= 0 
#r1y= 35.163  A r2y= 7.0326  A r3y= 0 
#r1z= 37.097  A r2z= 0  A r3z= 7.4194 
#wavelength= .9 A   number of atoms: 60 
#overall temperature factor: exp(-2* .3 A*(sin(theta)/lambda)^2)
#Lorentz Factor: 1 / sin^2(2theta)   neutron powder flat sample
#
#    x(r1)    y(r2)  z(r3)  mx(mb) my(mb) mz(mb) occupancy  sl     gJ   ff
#   0.000   0.000   0.000   0.000   2.102   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.000   0.500   0.077   0.000   2.102   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.100   0.000   0.000   0.000   2.145   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.100   0.500   0.077   0.000   2.145   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.200   0.000   0.000   0.000  -2.145   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.200   0.500   0.077   0.000  -2.145   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.300   0.000   0.000   0.000  -2.102   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.300   0.500   0.077   0.000  -2.102   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.400   0.000   0.000   0.000   2.116   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.400   0.500   0.077   0.000   2.116   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.500   0.000   0.000   0.000  -2.102   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.500   0.500   0.077   0.000  -2.102   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.600   0.000   0.000   0.000  -2.145   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.600   0.500   0.077   0.000  -2.145   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.700   0.000   0.000   0.000   2.145   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.700   0.500   0.077   0.000   2.145   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.800   0.000   0.000   0.000   2.102   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.800   0.500   0.077   0.000   2.102   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.900   0.000   0.000   0.000  -2.116   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.900   0.500   0.077   0.000  -2.116   0.000   1.000   0.769   0.727   0.054  25.029   0.310  12.102   0.658   4.722  -0.022   0.675  18.342   1.627   7.260   0.964   2.602   0.015#
#   0.000   0.302   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.000   0.698   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.000   0.802   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.000   0.198   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.100   0.302   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.100   0.698   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.100   0.802   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.100   0.198   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.200   0.302   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.200   0.698   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.200   0.802   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.200   0.198   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.300   0.302   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.300   0.698   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.300   0.802   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.300   0.198   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.400   0.302   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.400   0.698   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.400   0.802   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.400   0.198   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.500   0.302   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.500   0.698   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.500   0.802   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.500   0.198   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.600   0.302   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.600   0.698   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.600   0.802   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.600   0.198   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.700   0.302   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.700   0.698   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.700   0.802   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.700   0.198   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.800   0.302   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.800   0.698   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.800   0.802   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.800   0.198   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.900   0.302   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.900   0.698   0.702   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.900   0.802   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
#   0.900   0.198   0.372   0.000   0.000   0.000   1.000   0.772   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000#
# h       k       l      d[A]    |Q|[A^-1]    2theta    Inuc(2th) Imag(2th) Itot(2t) |sf|   lpg }
  0.000   0.000   0.000 100.00000   0.06283    0.000     0.000     0.000     0.000     0.000     0.000
  0.200   0.000   0.000  21.92150   0.28662    2.352     0.000     1.302     1.302     0.000   593.525
 -0.200   0.000   0.000  21.92150   0.28662    2.352     0.000     1.302     1.302     0.000   593.525
 -0.000   0.000   1.000   7.41940   0.84684    6.954     0.000     0.083     0.083     0.000    68.211
  0.000   0.000  -1.000   7.41940   0.84684    6.954     0.000     0.083     0.083     0.000    68.211
  0.600   0.000   0.000   7.30717   0.85985    7.061     0.000     0.987     0.987     0.000    66.170
 -0.600   0.000   0.000   7.30717   0.85985    7.061     0.000     0.987     0.987     0.000    66.170
  0.400   0.000   1.000   6.14412   1.02261    8.400     0.000     0.652     0.652     0.000    46.857
 -0.400   0.000  -1.000   6.14412   1.02261    8.400     0.000     0.652     0.652     0.000    46.857
 -0.400   0.000   1.000   6.14412   1.02261    8.400     0.000     0.652     0.652     0.000    46.857
  0.400   0.000  -1.000   6.14412   1.02261    8.400     0.000     0.652     0.652     0.000    46.857
  0.000   1.000   1.000   5.10407   1.23099   10.116     0.214     0.000     0.214     4.885    32.414
  0.000  -1.000  -1.000   5.10407   1.23099   10.116     0.214     0.000     0.214     4.885    32.414
  0.000   1.000  -1.000   5.10407   1.23099   10.116     0.214     0.000     0.214     4.885    32.414
 -0.000  -1.000   1.000   5.10407   1.23099   10.116     0.214     0.000     0.214     4.885    32.414
  0.800   0.000   1.000   4.40819   1.42531   11.718     0.000     0.047     0.047     0.000    24.243
 -0.800   0.000  -1.000   4.40819   1.42531   11.718     0.000     0.047     0.047     0.000    24.243
 -0.800   0.000   1.000   4.40819   1.42531   11.718     0.000     0.047     0.047     0.000    24.243
  0.800   0.000  -1.000   4.40819   1.42531   11.718     0.000     0.047     0.047     0.000    24.243
  1.000   0.000   0.000   4.38430   1.43308   11.782     0.000     0.030     0.030     0.000    23.984
 -1.000   0.000   0.000   4.38430   1.43308   11.782     0.000     0.030     0.030     0.000    23.984
  0.600   1.000   1.000   4.18436   1.50155   12.347     0.000     0.011     0.011     0.000    21.869
 -0.600  -1.000  -1.000   4.18436   1.50155   12.347     0.000     0.011     0.011     0.000    21.869
 -0.600   1.000   1.000   4.18436   1.50155   12.347     0.000     0.011     0.011     0.000    21.869
  0.600  -1.000  -1.000   4.18436   1.50155   12.347     0.000     0.011     0.011     0.000    21.869
  0.600  -1.000   1.000   4.18436   1.50155   12.347     0.000     0.011     0.011     0.000    21.869
 -0.600   1.000  -1.000   4.18436   1.50155   12.347     0.000     0.011     0.011     0.000    21.869
  0.600   1.000  -1.000   4.18436   1.50155   12.347     0.000     0.011     0.011     0.000    21.869
 -0.600  -1.000   1.000   4.18436   1.50155   12.347     0.000     0.011     0.011     0.000    21.869
  0.200   0.000   2.000   3.65770   1.71776   14.134     0.000     0.026     0.026     0.000    16.771
 -0.200   0.000  -2.000   3.65770   1.71776   14.134     0.000     0.026     0.026     0.000    16.771
 -0.200   0.000   2.000   3.65770   1.71776   14.134     0.000     0.026     0.026     0.000    16.771
  0.200   0.000  -2.000   3.65770   1.71776   14.134     0.000     0.026     0.026     0.000    16.771
  0.000   2.000   0.000   3.51630   1.78683   14.705     0.353    -0.000     0.353     9.108    15.519
 -0.000  -2.000   0.000   3.51630   1.78683   14.705     0.353    -0.000     0.353     9.108    15.519
  0.600   0.000   2.000   3.30783   1.89944   15.638     0.000     0.146     0.146     0.000    13.763
 -0.600   0.000  -2.000   3.30783   1.89944   15.638     0.000     0.146     0.146     0.000    13.763
 -0.600   0.000   2.000   3.30783   1.89944   15.638     0.000     0.146     0.146     0.000    13.763
  0.600   0.000  -2.000   3.30783   1.89944   15.638     0.000     0.146     0.146     0.000    13.763
  1.200   0.000   1.000   3.27772   1.91689   15.782     0.000     0.025     0.025     0.000    13.518
 -1.200   0.000  -1.000   3.27772   1.91689   15.782     0.000     0.025     0.025     0.000    13.518
  1.200   0.000  -1.000   3.27772   1.91689   15.782     0.000     0.025     0.025     0.000    13.518
 -1.200   0.000   1.000   3.27772   1.91689   15.782     0.000     0.025     0.025     0.000    13.518
 -0.600   2.000   0.000   3.16853   1.98295   16.330     0.000     0.032     0.032     0.000    12.650

etc ....
-----------------------------------------------------------------------
\end{verbatim}
.. such a reflection list can be easily transformed into a diffraction pattern
by convolution with the resolution function of a neutron diffractometer.
Use programs {\prg convolut} (windows only) or {\prg convolute} (windows and linux) for 
this task. Here is an example command

{prg convolut ELN.OUT 68 STP=0.1 MODE=GAUSS FWHM=0.25}

which gives the magnetic pattern shown in fig.~\ref{elnpattern}



