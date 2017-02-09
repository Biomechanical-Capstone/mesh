package mesh.io;
import java.io.DataInput;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.HashMap;

import mesh.math.BoundingBox;
import mesh.math.Vec3;
import mesh.objects.TriangleMesh;

public class STLReader {
	
	/**
     *  Import a new Scene object from a (BINARY STL) stream
     */
	public TriangleMesh importStream(InputStream is) throws IOException
	{
		int x, vertno=0, face[], faceArray[][] = new int[0][0], len;
		Vec3 vert, norm, calcNorm, vertArray[], v1, v2, v3;
		Integer index;
		ArrayList vlist = new ArrayList(1024*128);
		ArrayList flist = new ArrayList(1024*64);
		HashMap vmap = new HashMap(1024*128);
		byte[] buff = new byte[80];
		
		norm = new Vec3();
		vert = new Vec3();
		vertArray = new Vec3[0];
		faceArray = new int[0][0];
		face = new int[3];
		BoundingBox bounds = new BoundingBox(vert, vert);
		
		// binary STL is always little-endian
		DataInput in = new LittleEndianDataInputStream(is);
		
		try {
		
		    int count = 0;
		
			try {
			    in.readFully(buff);
			} catch (EOFException e) {
			    System.out.println("STL: " + e);
			}
		
			vlist.clear();
			flist.clear();
			vmap.clear();
		
			// get the number of faces
			len = in.readInt();
			//HAD TO SPECIFY NUMBER OF FACES FOR MY MESH
			//len = 3;
			System.out.println("STL; faces=" + len);
		
			count++;
				
			// read every face
			for (int faceno = 0; faceno < len; faceno++) {
		
			    //System.out.println("STL: face#" + faceno);
		
			    face = new int[3];
		
			    // read the normal
			    readVec(in, norm);
		
			    // read the 3 vertices
			    for (vertno = 0; vertno < 3; vertno++) {
				readVec(in, vert);
				index = (Integer) vmap.get(vert.toString());
		
				if (index != null) {
				    face[vertno] = index.intValue();
				}
				else {
				    if (vert.x < bounds.minx) bounds.minx = vert.x;
				    if (vert.x > bounds.maxx) bounds.maxx = vert.x;
				    if (vert.y < bounds.miny) bounds.miny = vert.y;
				    if (vert.y > bounds.maxy) bounds.maxy = vert.y;
				    if (vert.z < bounds.minz) bounds.minz = vert.z;
				    if (vert.z > bounds.maxz) bounds.maxz = vert.z;
		
				    x = vlist.size();
				    face[vertno] = x;
				    vmap.put(vert.toString(), new Integer(x));
				    vlist.add(vert);
		
				    vert = new Vec3();	// new Vec3
				}
			    }
		
			    // read padding
			    in.skipBytes(2);
		
			    // calculate the normal, and compare
			    v1 = (Vec3) vlist.get(face[0]);
			    v2 = (Vec3) vlist.get(face[1]);
			    v3 = (Vec3) vlist.get(face[2]);
		
			    calcNorm = v2.minus(v1).cross(v3.minus(v1));
			    double length = calcNorm.length();
			    if (length > 0.0) calcNorm.scale(1.0/length);
		
			    if (!calcNorm.equals(norm))
		
				System.out.println("\nWarning: facet normal mismatch." +
					      "read: " + norm + 
					      "; calculated: " + calcNorm);
				    
			    flist.add(face);
			}
		
			//System.out.println("STL: building mesh");
		
			// build the mesh
			faceArray = (int[][]) flist.toArray(faceArray);
			vertArray = (Vec3[]) vlist.toArray(vertArray);
		
			TriangleMesh mesh = new TriangleMesh(vertArray, faceArray);
		
			validate(mesh);
		
			//System.out.println("STL: new object added");
			//}
		
		    if (count == 0)
			System.out.println("\nNo object created");
		    
		    return mesh;
		}
		catch (Exception e) {
		
		    return null;
		}
		
	}
	
	/**
	 *  read a 3D vector from a (binary) STL stream.
	 */
	protected static Vec3 readVec(DataInput in, Vec3 vec) throws IOException
	{
		if (vec == null) vec = new Vec3();
		
		vec.x = in.readFloat();
		vec.y = in.readFloat();
		vec.z = in.readFloat();
		
		return vec;
	}
	
	/**
	 *  validate that a trianglemesh is consistent with STL
	 */
	public static boolean validate(TriangleMesh mesh)
	{
		boolean valid = true;
		
		
		    // number of faces must be even
		    int fc = mesh.getFaces().length;
		    if (((fc/2) * 2) != fc) {
			valid = false;
			System.out.println("\nvalidate: number of faces (" + fc +
				  ") is not even");
		    }
		
		    // number of edges must be a multiple of 3
		    int ec = mesh.getEdges().length;
		    if (((ec/3) * 3) != ec) {
			valid = false;
			System.out.println("\nvalidate: number of edges (" + ec +
				  ") is not a multiple of 3");
		    }
		
		    // number of edges must be 2/3 number of faces
		    if (2*ec != 3*fc) {
			valid = false;
			System.out.println("\nvalidate: number of faces (" + fc +
				  ") is not 2/3 the number of edges (" + ec + ")");
		    }
		
		    // calculate the number of holes we've found...
		    int vc = mesh.getVertices().length;
		    int holes = -((fc - ec + vc - 2) / 2);
		    if (holes > 0)
			System.out.println("\nvalidate: calculated holes (using Euler's) = " +
				  holes);
		
		    System.out.println("validate: fc=" +fc+ "; ec=" +ec+ "; vc=" +vc);
		
		return valid;
	}
}
