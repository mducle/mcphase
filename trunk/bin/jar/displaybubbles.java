import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import javachart.chart.*;
import java.io.*;
import java.lang.*;
import com.sun.image.codec.jpeg.JPEGCodec;
import com.sun.image.codec.jpeg.JPEGImageEncoder;

public class displaybubbles extends Panel implements Runnable {
 public BubbleChart chart = new BubbleChart("displaybubbles");
 Button bRot=new Button("save displaybubbles.jpg");                       //erstellt einen Button
 Thread myThread = null;
 static   double vals[]={0.,1.};
 static String[] file;
 static int[] colx;
 static int[] coly;
 static int[] colint;
 
 public void start(){ myThread = new Thread (this); myThread.start();}

 public void stop(){myThread = null;}

 public void run(){ while(myThread!=null){
    try{Thread.sleep(5000);
       }catch(Exception ignored){}
       // here do something

 File fileIni;
 try{
 for (int i=0;i<file.length;++i)
 {String s="";
      Dataset ds = chart.getDataset(file[i]+s.valueOf(i));
      ds.getData().removeAllElements();

 fileIni = new File(file[i]);

    //ï¿½ffnen der Datei
    DataInputStream inStream = new DataInputStream(new FileInputStream(fileIni));
    String strLine;
    String sx;
    String sy;
    String sint;
    int clx = colx[i];
    int cly = coly[i];   
    int clint = colint[i];

//    display app = new display();
//    app.setSize(640, 640);

//    Dataset ds = chart.getDataset("xy");
//    ds.getGc().setMarkerStyle(Gc.MK_DIAMOND);
//    ds.getGc().setMarkerSize(7);
    
    
//    ds.getGc().setLineWidth(0);
//    ds.getGc().setFillColor(Color.blue);

     //Auslesen der Datei
    while (inStream.available() > 0)
    {
      strLine = inStream.readLine();
      if ((strLine.length() == 0)
        ||(TrimString(strLine).substring(0, 1).equalsIgnoreCase("#")))
      {
        continue;
      }
      
      // select colx and coly
      sx=TrimString(strLine);
      sy=TrimString(strLine);
      sint=TrimString(strLine);
      int cx =clx-1;
      int cy =cly-1;      
      int cint =clint-1;      

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

      while (cint>0)
      {--cint;
       int iPos = sint.indexOf(" ");
       if (iPos < 0)
       {
         continue;
       }
       sint=sint.substring(iPos);
       sint=TrimString(sint); 
      }
       cx=sx.indexOf(" ");
       cy=sy.indexOf(" ");
       cint=sint.indexOf(" ");
       if (cx>0) {sx=sx.substring(0,cx);}
       if (cy>0) {sy=sy.substring(0,cy);}
       if (cint>0) {sint=sint.substring(0,cint);}

      Double p = new Double(0.0);
      
      System.out.println(sx+" "+sy+" "+sint);
//      p.valueOf(strLine);
//    double[] myDatax = {p.parseDouble(sx)};
//    double[] myDatay = {p.parseDouble(sy)};
   try{
      Datum d = new Datum(p.parseDouble(sx), p.parseDouble(sy), p.parseDouble(sint),0,null);
      //ds.getGc().setMarkerSize((int)p.parseDouble(sint));
      ds.addDatum(d);}
      catch(NumberFormatException e){System.exit(1);}
    }   
  
    //double[] myDatay = {stringToDouble(strLine,0),stringToDouble(strLine,0)};
    }
// double[] myDatax = {1, 3, 2, 3,33};
// double[] myDatay = {123, 432, 223, 345,33};

 
// app.setVisible(true);
//  repaint();

 }
 catch(EOFException e)
    {
      System.out.println("EOF: " + e.getLocalizedMessage());
      //EntSession.CWatch("Unplanmï¿½ï¿½iges 'End Of File' in DSN-Konfigurationsdatei!");
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
 String s="";
       chart.setZScale(0.05);
 for (int i=0;i<file.length;++i)
      { chart.addDataset(file[i]+s.valueOf(i),vals,vals);
   }  
    bRot.addActionListener(new ActionListener(){
    public void actionPerformed(ActionEvent ed){
    try{
         FileOutputStream fos=new FileOutputStream("displaybubbles.jpg");
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
  if (args.length<3)
  {System.out.println("- too few arguments...\n");
    System.out.println("  program displaybubbles - show and watch data file by viewing a xy graphic on screen\n\n");
    System.out.println("use as:  displaybubbles xcol ycol intcol filename [xcol1 ycol1 intcol filename1 ...]\n\n");
    System.out.println("         xcol,ycol ... column to be taken as x-, y- and intensity-axis\n");
    System.out.println("	 filename ..... filename of datafile\n\n");
  System.exit(0);
  }

  file = new String[args.length/3];
  colx = new int[args.length/3];
  coly = new int[args.length/3];
  colint = new int[args.length/3];
   Double p = new Double(0.0);
//      System.out.println(sx+" "+sy);
//      p.valueOf(strLine);
//    double[] myDatax = {};
 int j=0;
 String title="";
 for(int i=0; i<args.length-1;	i+=4)
 {file[j]=args[i+3];
  Integer pp;
  ss=args[i];
  colx[j]=p.valueOf(ss).intValue();
  ss=args[i+1];
  coly[j]=p.valueOf(ss).intValue();
  ss=args[i+2];
  colint[j]=p.valueOf(ss).intValue();
  ++j; 
 title=title+args[i]+" "+args[i+1]+" "+args[i+2]+" "+args[i+3]+" ";
}

 Frame myFrame = new Frame(title);
 displaybubbles myPanel = new displaybubbles();
	myFrame.addWindowListener(new WindowAdapter() {
	    public void windowClosing(WindowEvent e) {System.exit(0);}
	});
 myPanel.initChart();
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