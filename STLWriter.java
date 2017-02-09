import java.io.CharArrayWriter;
import java.io.DataInput;
import java.io.DataOutput;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;


import math.BoundingBox;
import math.Mat4;
import math.Vec3;
import objects.MeshVertex;
//import objects.ObjectInfo;
import objects.TriangleMesh;
import streams.LittleEndianDataOutputStream;

public class STLWriter {
	public static final String VERSION = "3.0";
	protected NumberFormat nf = NumberFormat.getInstance(Locale.US);
	public int decimal = 12;
	protected int action, faces=0, rendermode;
	public double surfError = 0.05;
	protected Mat4 move;
	/**
     *  export the Scene (in BINARY STL) to the specified stream.
     */
	  public void exportStream(TriangleMesh mesh, OutputStream os)
			throws IOException, InterruptedException
	    {
			System.out.println("export to stream");

			// find the bounds, and set the transform
			//findBounds(list);

			// binary STL is always little-endian
			DataOutput out = new LittleEndianDataOutputStream(os);

			nf.setMaximumFractionDigits(decimal);
			nf.setGroupingUsed(false);

			//ObjectInfo info = null;
			//TriangleMesh mesh = null;
			MeshVertex vert[] = null;
			TriangleMesh.Face face[] = null;
			Vec3 v;
			//Mat4 trans = null;
			double length;
			
			String hdr = ".";
			for (int i = 0; i < 79; i++) {
				hdr += ".";
			}

				out.writeBytes(hdr);
				out.writeInt(mesh.getFaceCount());

			    vert = mesh.getVertices();
			    face = mesh.getFaces();

			    //trans = info.coords.fromLocal().times(move);

			    // print all faces to file
			    for (int i = 0; i < face.length; i++) {
				v = vert[face[i].v2].r.minus(vert[face[i].v1].r).cross(vert[face[i].v3].r.minus(vert[face[i].v1].r));

				length = v.length();
				if (length > 0.0) v.scale(1.0/length);

				//System.out.println("STL; norm before trans=" + v);
				// = info.coords.fromLocal().timesDirection(v);

				writeVec(out, v, null);
				writeVec(out, vert[face[i].v1].r, null);
				writeVec(out, vert[face[i].v2].r, null);
				writeVec(out, vert[face[i].v3].r, null);

				// two byte padding (ho-hum...)
				out.writeBytes("  ");
			    }
//		}

			System.out.println("stream complete");
	    }
	  
	   /**
	     *  find the bounds of the list of meshes
	     */
//	    protected void findBounds(TriangleMesh mesh)
//	    {
//		BoundingBox bb = null;
//		double dx=0, dy=0, dz=0;
//		ObjectInfo info;
//		//TriangleMesh mesh = null;
//
//		faces = 0;
//
//		//int max = list.size();
//		//for (int x = 0; x < max; x++) {
//		//    info = (ObjectInfo) list.get(x);
//
//		    bb = info.getBounds().transformAndOutset(info.coords.fromLocal());
//
//		    //System.out.println("export: bounds=" + bb);
//
//		    if (bb.minx < dx) dx = bb.minx;
//		    if (bb.miny < dy) dy = bb.miny;
//		    if (bb.minz < dz) dz = bb.minz;
//
//		//    mesh = info.object.convertToTriangleMesh(surfError);
//		    if (mesh != null) faces += mesh.getFaces().length;
//		//}
//
//		// calculate the transform to ensure no negative numbers...
//		dx = (dx < 0 ? -dx : 0.0);
//		dy = (dy < 0 ? -dy : 0.0);
//		dz = (dz < 0 ? -dz : 0.0);
//		move = Mat4.translation(dx, dy, dz);
//	    }
	    
	    /**
	     *  write a 3D vector to a (binary) STL stream
	     */
	    protected void writeVec(DataOutput out, Vec3 p, Mat4 trans)
		throws IOException
	    {
		if (trans != null) p = trans.times(p);

		out.writeFloat((float) p.x);
		out.writeFloat((float) p.y);
		out.writeFloat((float) p.z);
	    }
}
