import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import javachart.chart.*;
import java.io.*;
import java.lang.*;
import java.util.StringTokenizer;
import com.sun.image.codec.jpeg.JPEGCodec;
import com.sun.image.codec.jpeg.JPEGImageEncoder;


public class display extends Panel implements Runnable {
 public LineChart chart = new LineChart("display");
 Button bRot=new Button("save display.jpg");                       //erstellt einen Button
 Thread myThread = null;
 static double vals[]={0.,1.};
 static String[] file;
 static int[] colx;
 static int[] coly;
 static String [] legend; 
 static String xText = "";
 static String yText = "";
 static FileInputStream ff;
 
 public void start(){ myThread = new Thread (this); myThread.start();}

 public void stop(){myThread = null;}

 public void run(){ while(myThread!=null){
    try{Thread.sleep(500);
       }catch(Exception ignored){}
       // here do something

 File fileIni;
 try{
 for (int i=0;i<file.length;++i)
 {String s="";
  legend[i]="false";
      Dataset ds = chart.getDataset(file[i]+s.valueOf(i));
      ds.getData().removeAllElements();

 fileIni = new File(file[i]);
 ff = new FileInputStream(fileIni);
    //ffnen der Datei
    DataInputStream inStream = new DataInputStream(ff);
    String strLine;
    String sx;
    String sy;
    int clx = colx[i];
    int cly = coly[i];   

//    display app = new display();
//    app.setSize(640, 640);

//    Dataset ds = chart.getDataset("xy");
    ds.getGc().setMarkerStyle(Gc.MK_DIAMOND);
    ds.getGc().setMarkerSize(7);
    ds.getGc().setLineWidth(1);
//    ds.getGc().setFillColor(Color.blue);

     //Auslesen der Datei
    while (inStream.available() > 0)
    {
      strLine = inStream.readLine();
      if ((strLine.length() == 0)
        ||(TrimString(strLine).substring(0, 1).equalsIgnoreCase("#")))
      {
      for(int i1=0;i1<=strLine.length();++i1)
       {if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylegend=true")){legend[i]="true";chart.setLegendVisible(true);}}
        if(i1<=strLine.length()-19){if(strLine.substring(i1,i1+19).equalsIgnoreCase("displaylegend=false")){legend[i]="false";chart.setLegendVisible(false);}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displayxtext=")){xText=strLine.substring(i1+13,strLine.length());}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displayytext=")){yText=strLine.substring(i1+13,strLine.length());}}
        if(i1<=strLine.length()-17){if(strLine.substring(i1,i1+17).equalsIgnoreCase("displaylines=true")){chart.setLineVisible(true);}}
        if(i1<=strLine.length()-18){if(strLine.substring(i1,i1+18).equalsIgnoreCase("displaylines=false")){chart.setLineVisible(false);}}
        if(i1<=strLine.length()-13){if(strLine.substring(i1,i1+13).equalsIgnoreCase("displaytitle=")){chart.getBackground().setTitleString(strLine.substring(i1+13,strLine.length()));}}
        }

        continue;
       }

     //      System.out.println(strLine);
     // replace tabs by spaces
      strLine=strLine.replaceAll("[\t\n\u000B\u0009\f]"," ");

    //      System.out.println(strLine);

      // select colx and coly
      sx=TrimString(strLine);
      sy=TrimString(strLine);
      int cx =clx-1;
      int cy =cly-1;      

      while (cx>0)
      {--cx;
       int iPos = sx.indexOf(" ");
       if (iPos < 0)
       {
         continue;
       }
       sx=sx.substring(iPos);
       sx=TrimString(sx); 
      }

      while (cy>0)
      {--cy;
       int iPos = sy.indexOf(" ");
       if (iPos < 0)
       {
         continue;
       }
       sy=sy.substring(iPos);
       sy=TrimString(sy); 
      }
       cx=sx.indexOf(" ");
       cy=sy.indexOf(" ");
       if (cx>0) {sx=sx.substring(0,cx);}
       if (cy>0) {sy=sy.substring(0,cy);}

       sx=sx.replace('D','E');
       sy=sy.replace('D','E');

      Double p = new Double(0.0);
//      System.out.println(sx+" "+sy);
//      p.valueOf(strLine);
//    double[] myDatax = {p.parseDouble(sx)};
//    double[] myDatay = {p.parseDouble(sy)};
   try{
      Datum d = new Datum(p.parseDouble(sx),p.parseDouble(sy),null);
      ds.addDatum(d);}
      catch(NumberFormatException e){;}
    }
    ff.close();
    //double[] myDatay = {stringToDouble(strLine,0),stringToDouble(strLine,0)};
    }
// double[] myDatax = {1, 3, 2, 3,33};
// double[] myDatay = {123, 432, 223, 345,33};

   chart.getXAxis().setTitleString(xText); 
   chart.getYAxis().setTitleString(yText); 

 
// app.setVisible(true);
//  repaint();

 }
 catch(EOFException e)
    {
      System.out.println("EOF: " + e.getLocalizedMessage());
      //EntSession.CWatch("Unplanned 'End Of File' in DSN-Konfigurationsdatei!");
    }

    catch (FileNotFoundException e)
    {
      System.out.println("File not found: " + e.getLocalizedMessage());
      //EntSession.CWatch("Konfigurationsdatei cti_listener.ini nicht gefunden!");
    }

    //Sonstiger Dateifehler
    catch (IOException e)
    {
      System.out.println("Dateifehler: " + e.getLocalizedMessage());
      //EntSession.CWatch("Fehler beim Zugriff auf Datei cti_listener.ini!");
    }

 

       //  
        repaint();



 }}


 protected void initChart(){ 
 // Give it a title
 //     chart.getBackground().setTitleFont(new Font("Serif", Font.PLAIN, 24));
 //     chart.getBackground().setTitleString("Comparing Apples and Oranges");
chart.getXAxis().setTitleString("hallo"); 
chart.getYAxis().setTitleString(yText); 
 String s="";
 for (int i=0;i<file.length;++i)
   {
   chart.addDataset(file[i]+s.valueOf(i),vals,vals);
   }  

    bRot.addActionListener(new ActionListener(){
    public void actionPerformed(ActionEvent ed){
    try{
         FileOutputStream fos=new FileOutputStream("display.jpg");
         BufferedImage image= new BufferedImage(chart.getWidth(),chart.getHeight(), BufferedImage.TYPE_INT_RGB); 
         Graphics g=image.getGraphics();
         paint(g);
         JPEGImageEncoder encoder= JPEGCodec.createJPEGEncoder(fos); 
         encoder.encode(image);
         fos.close();
    }    catch (FileNotFoundException e)
    {
         System.out.println("File not found: " + e.getLocalizedMessage());
         //EntSession.CWatch("Konfigurationsdatei cti_listener.ini nicht gefunden!");
    }
         //Sonstiger Dateifehler
         catch (IOException e)
    {
         System.out.println("Dateifehler: " + e.getLocalizedMessage());
         //EntSession.CWatch("Fehler beim Zugriff auf Datei cti_listener.ini!");
    }

      }
                                                 });
 }



 public void update(Graphics g){paint(g);}

 public void paint(Graphics g){try{chart.paint(this,g);}catch(ArrayIndexOutOfBoundsException e){;}}

 public static void main(String[] args){ 
  String ss;
  if (args.length<2)
  {System.out.println("- too few arguments...\n");
    System.out.println("  program display - show and watch data file by viewing a xy graphic on screen\n\n");
    System.out.println("use as:  display xcol ycol filename [xcol1 ycol1 filename1 ...]\n\n");
    System.out.println("         xcol,ycol ... column to be taken as x- and y-axis\n");
    System.out.println("	 filename ..... filename of datafile\n\n");
    System.out.println("	 Data files may contain lines to tune the display output, such as\n");
    System.out.println("	 # displayytext=intensity\n");
    System.out.println("	 # displayxtext=meV \n");
    System.out.println("	 # displaylines=true \n");
    System.out.println("	 # displaytitle=My new Graph\n");
    System.out.println("	 # displaylegend=false \n\n");
  System.exit(0);
  }

  file = new String[args.length/3];
  legend = new String[args.length/3];
  colx = new int[args.length/3];
  coly = new int[args.length/3];
   Double p = new Double(0.0);
//      System.out.println(sx+" "+sy);
//      p.valueOf(strLine);
//    double[] myDatax = {};
 int j=0;
 String title="";
 for(int i=0; i<args.length-1;	i+=3)
 {file[j]=args[i+2];
  legend[j]="false";
  Integer pp;
  ss=args[i];
  colx[j]=p.valueOf(ss).intValue();
  ss=args[i+1];
  coly[j]=p.valueOf(ss).intValue();
  ++j; 
 title=title+args[i]+" "+args[i+1]+" "+args[i+2]+" ";
}

 display myPanel = new display();
         myPanel.initChart();
 

 
 Frame myFrame = new Frame(title);
       myFrame.addWindowListener(new WindowAdapter() {
	    public void windowClosing(WindowEvent e) {System.exit(0);}
	});
       myFrame.add(myPanel.bRot);                                          //fügt dem JFrame den Button hinzu
       myFrame.pack();
       myFrame.add(myPanel);
       myFrame.setSize(300,300);
       myFrame.setVisible(true);
       myPanel.start();	

//	frame = new JFrame("FileChooserDemo");
//	frame.addWindowListener(new WindowAdapter() {
//	    public void windowClosing(WindowEvent e) {System.exit(0);}
//	});
//	frame.getContentPane().add("Center", panel);
//	frame.pack();
//	frame.setVisible(true);

//	panel.updateState();
 }

 static private String FirstWord(String strSource)
 {String fw;
  fw=TrimString(strSource);
       int iPos = fw.indexOf(" ");
       if (iPos >= 0)
       {
       fw=strSource.substring(0,iPos);
       fw=TrimString(fw); 
       }
 return(fw); 
 }

 static private String DropWord(String strSource)
 {String fw;
  fw=TrimString(strSource);
       int iPos = fw.indexOf(" ");
       if (iPos >= 0)
       {
       fw=strSource.substring(iPos);
       fw=TrimString(fw); 
       }
 return(fw); 
 }

 static private String TrimString(String strSource)
 {
    while ((strSource.startsWith(" "))
      && (strSource.length() > 0))
      {
        strSource = strSource.substring(1, strSource.length());
      }

    while ((strSource.endsWith(" "))
        && (strSource.length() > 0))
      {
        strSource = strSource.substring(0, strSource.length() - 1);
      }

    return(strSource);
 }




}