import java.awt.*;
import java.awt.event.*;
import javachart.chart.*;
import java.io.*;
import java.lang.*;

public class displaymag extends Panel implements Runnable {
 public LineChart chart = new LineChart("sample");
 Thread myThread = null;
 static   double vals[]={0.,1.};
 static String[] file;
 static int[] colx;
 static int[] coly;
 static int ctr=0;
 
 public void start(){ myThread = new Thread (this); myThread.start();}

 public void stop(){myThread = null;}

 public void run(){ while(myThread!=null&&ctr<=10){
    try{Thread.sleep(500);
       }catch(Exception ignored){}
       // here do something

 File fileIni;
 try{
 for (int i=0;i<file.length;++i)
 {String s="abcdefghikl";
      Dataset ds = chart.getDataset(s.substring(i,i+1));
      ds.getData().removeAllElements();

 fileIni = new File(file[i]);

    //ffnen der Datei
    DataInputStream inStream = new DataInputStream(new FileInputStream(fileIni));
    String strLine;
    String sx;
    String sy;
    int clx = colx[i];
    int cly = coly[i];   

//    displaymag app = new displaymag();
//    app.setSize(640, 640);

//    Dataset ds = chart.getDataset("xy");
    ds.getGc().setMarkerStyle(Gc.MK_DIAMOND);
    ds.getGc().setMarkerSize(7);
//    ds.getGc().setLineWidth(0);
//    ds.getGc().setFillColor(Color.blue);

     //Auslesen der Datei
    while (inStream.available() > 0)
    { strLine = inStream.readLine();
        ctr=0;
      if ((strLine.length() == 0)
        ||(strLine.substring(0, 1).equalsIgnoreCase("#")))
      {
        continue;
      }
      
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

      Double p = new Double(0.0);
//      System.out.println(sx+" "+sy);
//      p.valueOf(strLine);
//    double[] myDatax = {p.parseDouble(sx)};
//    double[] myDatay = {p.parseDouble(sy)};
 try{
      Datum d = new Datum(p.parseDouble(sx),p.parseDouble(sy),null);
      ds.addDatum(d);
      }
      catch(NumberFormatException e){;}
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
      //EntSession.CWatch("Unplanned 'End Of File' in DSN-Konfigurationsdatei!");
    }

    catch (FileNotFoundException e)
    {++ctr;
      System.out.println("File not found: " + e.getLocalizedMessage());
      //EntSession.CWatch("Konfigurationsdatei cti_listener.ini nicht gefunden!");
    }

    //Sonstiger Dateifehler
    catch (IOException e)
    {++ctr;
      System.out.println("Dateifehler: " + e.getLocalizedMessage());
      //EntSession.CWatch("Fehler beim Zugriff auf Datei cti_listener.ini!");
    }

 

       //  
        repaint();
 }}

 protected void initChart(){ 
    chart.setLineVisible(false);
    chart.setLegendVisible(true);
    String s="abcdefghijkl";
 for (int i=0;i<file.length;++i)
   {//char ii=i;   
    chart.addDataset(s.substring(i,i+1),vals,vals);
   }  
 }


 
 public void update(Graphics g){paint(g);}

 public void paint(Graphics g){try{chart.paint(this,g);}catch(ArrayIndexOutOfBoundsException e){;}}

 public static void main(String[] args){ 
  String ss;
  file = new String[args.length/3];
  colx = new int[args.length/3];
  coly = new int[args.length/3];
   Double p = new Double(0.0);
//      System.out.println(sx+" "+sy);
//      p.valueOf(strLine);
//    double[] myDatax = {};
 int j=0;
 String title="Magnetisation";
 for(int i=0; i<args.length-1;	i+=3)
 {file[j]=args[i+2];
  Integer pp;
  ss=args[i];
  colx[j]=p.valueOf(ss).intValue();
  ss=args[i+1];
  coly[j]=p.valueOf(ss).intValue();
  ++j; 
// title=title+args[i]+" "+args[i+1]+" "+args[i+2]+" ";
}

 Frame myFrame = new Frame(title);
 displaymag myPanel = new displaymag();
	myFrame.addWindowListener(new WindowAdapter() {
	    public void windowClosing(WindowEvent e) {System.exit(0);}
	});
 myPanel.initChart();
 myFrame.add(myPanel);
 myFrame.setSize(400,400);
 myFrame.setVisible(true);
 myPanel.start();	
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
