/* 
 * Viewer.java
 */

import java.awt.*;
import java.awt.event.*;

public class Viewer extends Frame {
	private Image image;

	public Viewer(String fileName) {
		Toolkit toolkit = Toolkit.getDefaultToolkit();
		image = toolkit.getImage(fileName);
		MediaTracker mediaTracker = new MediaTracker(this);
		mediaTracker.addImage(image, 0);
		try
		{
			mediaTracker.waitForID(0);
		}
		catch (InterruptedException ie)
		{
			System.err.println(ie);
			System.exit(1);
		}
		addWindowListener(new WindowAdapter() {
      		public void windowClosing(WindowEvent e) {
        		System.exit(0);
      		}
		});
		setSize(image.getWidth(null), image.getHeight(null));
		setTitle(fileName);
		show();
	}

	public void paint(Graphics graphics) {
		graphics.drawImage(image, 0, 0, null);
	}

	public static void main(String[] args) {
		new Viewer(args[0]);
	}
}