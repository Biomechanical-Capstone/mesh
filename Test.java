import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import objects.TriangleMesh;

public class Test {
	public static void main(String args[]) throws FileNotFoundException, IOException {
		STLReader reader = new STLReader();
		BufferedInputStream is = new BufferedInputStream(new FileInputStream("test.stl"));
		TriangleMesh mesh = reader.importStream(is);
		is.close();
	}
}
