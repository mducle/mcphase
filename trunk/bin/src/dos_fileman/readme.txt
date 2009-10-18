{\footnotesize
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
                      new data set.
		      
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

[deactivated and substituted by perl script of the same name] FFORM *.* 3 ##.##   - formats 3rd column of file with basic number format
[deactivated and substituted by perl script of the same name] DELCOL *.* 25       - deletes columns 2 to 5 from file *.*
[deactivated and substituted by perl script of the same name] NEWCOL *.* 3        - inserts new column nr.3 in file *.* (line number col)
SWAP *.* 13         - swap column 1 and 3
[deactivated and substituted by perl script of the same name] SWAPF *.* *.*       - swaps the first columns of two files


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

[deactivated and substituted by perl script of the same name] POT *.* 2 24.1      - col2=col2^24.1
[deactivated and substituted by perl script of the same name] EXP *.* 1           - col1=exp(col1)

INVERT *.* 3        - invert columns 2
[deactivated and substituted by perl script of the same name] FACT *.* 3 2.42     - multiply col. 2 with the constant 2.42
[deactivated and substituted by perl script of the same name] SHIFTC *.* 3 3.47   - add 3.47 to column 3 option /+ shifts more files with
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
[deactivated and substituted by perl script of the same name] ADDC *.* 13           - column 3 of file *.* is added to column 1 point by
                        point - the result is written in column 1
[deactivated and substituted by perl script of the same name] MULTC *.* 13          - column 3 of file *.* is multiplied with col 1 point by
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
                        pattern from the output of program mcdiff

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
\end{verbatim}
}
