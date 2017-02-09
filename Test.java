import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import objects.TriangleMesh;

public class Test {
	public static void main(String args[]) throws FileNotFoundException, IOException, InterruptedException {
		STLReader reader = new STLReader();
		BufferedInputStream is = new BufferedInputStream(new FileInputStream("test.stl"));
		TriangleMesh mesh = reader.importStream(is);
		is.close();
		if(mesh != null)
		{
			STLWriter writer = new STLWriter();
			OutputStream os = new FileOutputStream("out.stl");
			writer.exportStream(mesh, os);
			os.close();
		}
		else
			System.out.println("MESH RETURNED NULL");
	}
}
