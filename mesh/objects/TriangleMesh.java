package mesh.objects;
/* Copyright (C) 1999-2015 by Peter Eastman

   This program is free software; you can redistribute it and/or modify it under the
   terms of the GNU General Public License as published by the Free Software
   Foundation; either version 2 of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
   PARTICULAR PURPOSE.  See the GNU General Public License for more details. */

import java.awt.*;
import java.util.*;

import mesh.math.BoundingBox;
import mesh.math.Vec3;

/** The TriangleMesh class represents an aritrary surface defined by a mesh of triangular
    faces.  Depending on the selected smoothing method, the surface may simply consist of
    the triangular faces, or it may be a smooth subdivision surface which either interpolates
    or approximates the vertices of the control mesh. */

public class TriangleMesh
{
	
  /** A vertex specifies a position vector, the number of edges which share the vertex, and
      the "first" edge.  If the vertex is in the interior of the mesh, any edge can be the
      first one.  If it is on the boundary, then the first edge must be one of the two boundary
      edges.  A vertex also has a "smoothness" parameter associated with it. */

  public class Vertex extends MeshVertex
  {
    public int edges, firstEdge;
    public float smoothness;

    public Vertex(Vec3 p)
    {
      super(p);
      edges = 0;
      firstEdge = -1;
      smoothness = 1.0f;
    }

    public Vertex(Vertex v)
    {
      super(v);
      edges = v.edges;
      firstEdge = v.firstEdge;
      smoothness = v.smoothness;
    }

    /** Make this vertex identical to another one. */

    public void copy(Vertex v)
    {
      r.set(v.r);
      edges = v.edges;
      firstEdge = v.firstEdge;
      smoothness = v.smoothness;
      ikJoint = v.ikJoint;
      ikWeight = v.ikWeight;
    }

    /** Multiple the fields of this vertex by a constant. */

    public void scale(double d)
    {
      r.scale(d);
      smoothness *= d;
      ikWeight *= d;
    }

    /** Set the various fields to zero. */

    public void clear()
    {
      r.set(0.0, 0.0, 0.0);
      smoothness = 0.0f;
      ikWeight = 0.0;
    }

    /** Construct a list of all edges which share the vertex. */

    public int[] getEdges()
    {
      int e[] = new int [edges], i;
      if (edges == 0)
        return e;
      Face f = face[edge[firstEdge].f1];
      e[0] = firstEdge;
      for (i = 1; i < edges; i++)
        {
          if (vertex[f.v1] == this)
            {
              if (f.e1 == e[i-1])
                e[i] = f.e3;
              else
                e[i] = f.e1;
            }
          else if (vertex[f.v2] == this)
            {
              if (f.e1 == e[i-1])
                e[i] = f.e2;
              else
                e[i] = f.e1;
            }
          else
            {
              if (f.e2 == e[i-1])
                e[i] = f.e3;
              else
                e[i] = f.e2;
            }
          if (face[edge[e[i]].f1] == f && edge[e[i]].f2 > -1)
            f = face[edge[e[i]].f2];
          else
            f = face[edge[e[i]].f1];
        }
      return e;
    }

    /** This method tells whether the list of edges returned by getEdges() are ordered
        clockwise or counter-clockwise. */

    public boolean clockwise()
    {
      Face f = face[edge[firstEdge].f1];

      if (f.e1 == firstEdge)
        return (vertex[f.v2] == this);
      else if (f.e2 == firstEdge)
        return (vertex[f.v3] == this);
      else
        return (vertex[f.v1] == this);
    }
  }

  /** An edge is defined by the two vertices which it connects, and the two faces it is
      adjacent to.  For a boundary edge, f2 must equal -1.  An edge also has a "smoothness"
      parameter. */

  public class Edge
  {
    public int v1, v2, f1, f2;
    public float smoothness;

    public Edge(int vertex1, int vertex2, int face1)
    {
      v1 = vertex1;
      v2 = vertex2;
      f1 = face1;
      f2 = -1;
      smoothness = 1.0f;
    }
  }

  /** A face is defined by its three vertices and three edges.  The vertices must be arranged
      in counter-clockwise order, when viewed from the outside.  Edges 1, 2, and 3 connect
      vertices 1 and 2, 2 and 3, and 3 and 1 respectively. */

  public class Face
  {
    public int v1, v2, v3, e1, e2, e3;

    public Face(int vertex1, int vertex2, int vertex3, int edge1, int edge2, int edge3)
    {
      v1 = vertex1;
      v2 = vertex2;
      v3 = vertex3;
      e1 = edge1;
      e2 = edge2;
      e3 = edge3;
    }

    /** Given another face, return the index of the edge it shares with this one, or -1 if
        they do not share an edge. */

    public int getSharedFace(Face f)
    {
      if (f.e1 == e1 || f.e2 == e1 || f.e3 == e1)
        return e1;
      if (f.e1 == e2 || f.e2 == e2 || f.e3 == e2)
        return e2;
      if (f.e1 == e3 || f.e2 == e3 || f.e3 == e3)
        return e3;
      return -1;
    }
  }

  /* Beginning of TriangleMesh's variables and methods. */

  private Vertex vertex[];
  private Edge edge[];
  private Face face[];
  private boolean closed;
  private BoundingBox bounds;

  private static double LOOP_BETA[], BUTTERFLY_COEFF[][];

  /* The following constants are used during subdivision for recording parameter types. */

  private static final int PER_VERTEX = 1;
  private static final int PER_FACE = 2;
  private static final int PER_FACE_VERTEX = 3;

  /* Precalculate coefficients for Loop and Butterfly subdivision. */

  static {
    double beta;
    int i, j;

    LOOP_BETA = new double [32];
    for (i = 3; i < LOOP_BETA.length; i++)
      {
        beta = 0.375+0.25*Math.cos(2.0*Math.PI/i);
        LOOP_BETA[i] = (0.625-beta*beta)/i;
      }
    BUTTERFLY_COEFF = new double [32][];
    for (i = 5; i < BUTTERFLY_COEFF.length; i++)
      {
        BUTTERFLY_COEFF[i] = new double [i+1];
        BUTTERFLY_COEFF[i][i] = 1.0;
        beta = 2.0*Math.PI/i;
        for (j = 0; j < i; j++)
          {
            BUTTERFLY_COEFF[i][j] = (0.25+Math.cos(beta*j)+0.5*Math.cos(2.0*beta*j))/i;
            BUTTERFLY_COEFF[i][i] -= BUTTERFLY_COEFF[i][j];
          }
      }
    BUTTERFLY_COEFF[3] = new double [] {5.0/12.0, -1.0/12.0, -1.0/12.0, 0.75};
    BUTTERFLY_COEFF[4] = new double [] {.375, 0.0, -0.125, 0.0, 0.75};
    BUTTERFLY_COEFF[6] = new double [] {1.0, 0.125, -0.125, 0.0, -0.125, 0.125, 0.0};
  }

  /** The constructor takes three arguments.  v[] is an array containing the vertices.
      faces[][] is an N by 3 array containing the indices of the vertices which define each
      face.  The vertices for each face must be listed in order, such that they go
      counter-clockwise when viewed from the outside of the mesh.  All faces must have a
      consistent vertex order, such that the object has a well defined outer surface.
      This is true even if the mesh does not form a closed surface.  It is an error to
      call the constructor with a faces[][] array which does not meet this condition, and
      the results are undefined.  norm[] contains the normal vector at each vertex.  If
      any element of norm[] is null, flat shading will be used around that vertex. */

  public TriangleMesh(Vec3 v[], int faces[][])
  {
    Vertex vt[] = new Vertex [v.length];
    for (int i = 0; i < v.length; i++)
      vt[i] = new Vertex(v[i]);
    setShape(vt, faces);
  }

  public TriangleMesh(Vertex v[], int faces[][])
  {
    setShape(v, faces);
  }

  protected TriangleMesh()
  {
  }

  /** Construct the list of edges. */

  void findEdges(int faces[][])
  {
    int i, numEdges1 = 0, numEdges2 = 0, numCopied = 0;
    int faceEdges[][] = new int [faces.length][3], copiedEdges[] = new int [faces.length*3];
    Edge edges1[] = new Edge [faces.length*3], edges2[] = new Edge [faces.length*3];

    // If the mesh is closed, then each edge should be traversed twice, once in each
    // direction.  If the mesh is open, some edges will be traversed only once, which
    // could be in either direction.

    closed = true;
    for (i = 0; i < faces.length; i++)
      {
        if (faces[i][0] > faces[i][1])
          edges1[faceEdges[i][0] = numEdges1++] = new Edge(faces[i][0], faces[i][1], i);
        else
          {
            edges2[faceEdges[i][0] = numEdges2++] = new Edge(faces[i][0], faces[i][1], i);
            faceEdges[i][0] += edges1.length;
          }
        if (faces[i][1] > faces[i][2])
          edges1[faceEdges[i][1] = numEdges1++] = new Edge(faces[i][1], faces[i][2], i);
        else
          {
            edges2[faceEdges[i][1] = numEdges2++] = new Edge(faces[i][1], faces[i][2], i);
            faceEdges[i][1] += edges1.length;
          }
        if (faces[i][2] > faces[i][0])
          edges1[faceEdges[i][2] = numEdges1++] = new Edge(faces[i][2], faces[i][0], i);
        else
          {
            edges2[faceEdges[i][2] = numEdges2++] = new Edge(faces[i][2], faces[i][0], i);
            faceEdges[i][2] += edges1.length;
          }
      }
    if (numEdges1 != numEdges2)
      closed = false;

    // We now have two lists of edges: one for each direction of traversal.  Determine which
    // which ones are duplicates, and add any unique edges from edges2 into edges1.

    Hashtable<Point, Integer> edgeTable = new Hashtable<Point, Integer>();
    for (i = 0; i < numEdges1; i++)
      edgeTable.put(new Point(edges1[i].v1, edges1[i].v2), i);
    for (i = 0; i < numEdges2; i++)
      {
        Integer index = edgeTable.get(new Point(edges2[i].v2, edges2[i].v1));
        if (index == null)
          {
            copiedEdges[i] = numEdges1+numCopied++;
            edges1[copiedEdges[i]] = edges2[i];
          }
        else
          {
            copiedEdges[i] = index;
            edges1[index].f2 = edges2[i].f1;
          }
      }
    if (numCopied > 0)
      closed = false;

    // Record the edges for each face.

    for (i = 0; i < faces.length; i++)
      {
        if (faceEdges[i][0] >= edges1.length)
          faceEdges[i][0] = copiedEdges[faceEdges[i][0]-edges1.length];
        if (faceEdges[i][1] >= edges1.length)
          faceEdges[i][1] = copiedEdges[faceEdges[i][1]-edges1.length];
        if (faceEdges[i][2] >= edges1.length)
          faceEdges[i][2] = copiedEdges[faceEdges[i][2]-edges1.length];
      }

    // Construct the edges and faces.

    edge = new Edge [numEdges1+numCopied];
    for (i = 0; i < numEdges1+numCopied; i++)
      edge[i] = edges1[i];
    face = new Face [faces.length];
    for (i = 0; i < faces.length; i++)
      face[i] = new Face(faces[i][0], faces[i][1], faces[i][2], faceEdges[i][0], faceEdges[i][1], faceEdges[i][2]);
  }

  /** Calculate the (approximate) bounding box for the mesh. */

  private void findBounds()
  {
    double minx, miny, minz, maxx, maxy, maxz;
    Vec3 vert[] = null;
    int i;
    if (vert.length == 0)
      minx = maxx = miny = maxy = minz = maxz = 0.0;
    else
    {
      minx = maxx = vert[0].x;
      miny = maxy = vert[0].y;
      minz = maxz = vert[0].z;
      for (i = 1; i < vert.length; i++)
      {
        if (vert[i].x < minx) minx = vert[i].x;
        if (vert[i].x > maxx) maxx = vert[i].x;
        if (vert[i].y < miny) miny = vert[i].y;
        if (vert[i].y > maxy) maxy = vert[i].y;
        if (vert[i].z < minz) minz = vert[i].z;
        if (vert[i].z > maxz) maxz = vert[i].z;
      }
    }
    bounds = new BoundingBox(minx, maxx, miny, maxy, minz, maxz);
  }

  /** Get the bounding box for the mesh.  This is always the bounding box for the unsmoothed
      control mesh.  If the smoothing method is set to approximating, the final surface may not
      actually touch the sides of this box.  If the smoothing method is set to interpolating,
      the final surface may actually extend outside this box. */

  public BoundingBox getBounds()
  {
    if (bounds == null)
      findBounds();
    return bounds;
  }

  /** These methods return the lists of vertices, edges, and faces for the mesh. */

  
  public MeshVertex[] getVertices()
  {
    return vertex;
  }

  public Vertex getVertex(int i)
  {
    return vertex[i];
  }

  public Edge[] getEdges()
  {
    return edge;
  }

  public Face[] getFaces()
  {
    return face;
  }


  /** Get a list of the positions of all vertices which define the mesh. */

  public Vec3 [] getVertexPositions()
  {
    Vec3 v[] = new Vec3 [vertex.length];
    for (int i = 0; i < v.length; i++)
      v[i] = new Vec3(vertex[i].r);
    return v;
  }

  /** Set the positions for all the vertices of the mesh. */

  public void setVertexPositions(Vec3 v[])
  {
    for (int i = 0; i < v.length; i++)
      vertex[i].r = v[i];
    bounds = null;
  }

  /** This method rebuilds the mesh based on new lists of vertices and faces.  The smoothness
      values for all edges are lost in the process. */

  public void setShape(Vertex v[], int faces[][])
  {
    Vertex v1, v2;
    int i;

    // Create the vertices and edges.

    vertex = new Vertex [v.length];
    for (i = 0; i < v.length; i++)
      {
        vertex[i] = new Vertex(v[i]);
        vertex[i].firstEdge = -1;
        vertex[i].edges = 0;
      }
    if (faces.length == 0)
      {
        edge = new Edge [0];
        face = new Face [0];
      }
    else
      findEdges(faces);
    bounds = null;

    // Find the edge information for vertices.

    for (i = 0; i < edge.length; i++)
      {
        v1 = vertex[edge[i].v1];
        v2 = vertex[edge[i].v2];
        v1.edges++;
        v2.edges++;
        if (edge[i].f2 == -1)
          v1.firstEdge = v2.firstEdge = i;
        else
          {
            if (v1.firstEdge == -1)
              v1.firstEdge = i;
            if (v2.firstEdge == -1)
              v2.firstEdge = i;
          }
      }
  }

  public boolean isClosed()
  {
    return closed;
  }

  public void setSize(double xsize, double ysize, double zsize)
  {
    Vec3 size = getBounds().getSize();
    double xscale, yscale, zscale;

    if (size.x == 0.0)
      xscale = 1.0;
    else
      xscale = xsize / size.x;
    if (size.y == 0.0)
      yscale = 1.0;
    else
      yscale = ysize / size.y;
    if (size.z == 0.0)
      zscale = 1.0;
    else
      zscale = zsize / size.z;
    for (int i = 0; i < vertex.length; i++)
      {
        vertex[i].r.x *= xscale;
        vertex[i].r.y *= yscale;
        vertex[i].r.z *= zscale;
      }
    if (xscale*yscale*zscale < 0.0)
      reverseNormals();
    bounds = null;
  }

  /** Calculate a set of array representing the boundaries of this mesh.  There is one array
      for each distinct boundary, containing the indices of the edges which form that
      boundary. */

  public int [][] findBoundaryEdges()
  {
    // First, find every edge which is on a boundary.

    Vector<Integer> allEdges = new Vector<Integer>();
    for (int i = 0; i < edge.length; i++)
      if (edge[i].f2 == -1)
        allEdges.addElement(i);

    // Form boundaries one at a time.

    Vector<Vector<Integer>> boundary = new Vector<Vector<Integer>>();
    while (allEdges.size() > 0)
      {
        // Take one edge as a starting point, and follow around.

        Vector<Integer> current = new Vector<Integer>();
        Integer start = allEdges.elementAt(0);
        allEdges.removeElementAt(0);
        current.addElement(start);
        int i = start, j = 0;
        while (j < (allEdges.size()))
          {
            for (j = 0; j < allEdges.size(); j++)
              {
                int k = allEdges.elementAt(j);
                if (edge[i].v1 == edge[k].v1 || edge[i].v1 == edge[k].v2 ||
                    edge[i].v2 == edge[k].v1 || edge[i].v2 == edge[k].v2)
                  {
                    current.addElement(allEdges.elementAt(j));
                    allEdges.removeElementAt(j);
                    i = k;
                    j--;
                    break;
                  }
              }
          }
        boundary.addElement(current);
      }

    // Build the final arrays.

    int index[][] = new int [boundary.size()][];
    for (int i = 0; i < index.length; i++)
      {
        Vector<Integer> current = boundary.elementAt(i);
        index[i] = new int [current.size()];
        for (int j = 0; j < index[i].length; j++)
          index[i][j] = current.elementAt(j);
      }
    return index;
  }

  /** Create a vertex which is a blend of two existing ones. */

  private Vertex blend(Vertex v1, Vertex v2, double w1, double w2)
  {
    return new Vertex (new Vec3(w1*v1.r.x + w2*v2.r.x, w1*v1.r.y + w2*v2.r.y, w1*v1.r.z + w2*v2.r.z));
  }

  /** Create a vertex which is a blend of three existing ones. */

  private Vertex blend(Vertex v1, Vertex v2, Vertex v3, double w1, double w2, double w3)
  {
    return new Vertex (new Vec3(w1*v1.r.x + w2*v2.r.x + w3*v3.r.x, w1*v1.r.y + w2*v2.r.y + w3*v3.r.y, w1*v1.r.z + w2*v2.r.z + w3*v3.r.z));
  }

  /** Set a vertex to be a blend of two other ones. */

  private static void setBlend(Vertex v, Vertex v1, Vertex v2, double w1, double w2)
  {
    v.r.set(w1*v1.r.x + w2*v2.r.x, w1*v1.r.y + w2*v2.r.y, w1*v1.r.z + w2*v2.r.z);
  }

  /** Given a pair of vertices and a new vertex that is to be created between them, find the
      IK binding parameters for the new vertex. */

  private static void blendIKParams(Vertex newvert, Vertex v1, Vertex v2)
  {
    if (v1.ikJoint == v2.ikJoint)
      {
        newvert.ikJoint = v1.ikJoint;
        newvert.ikWeight = 0.5*(v1.ikWeight+v2.ikWeight);
      }
    else if (v1.ikWeight > v2.ikWeight)
      {
        newvert.ikJoint = v1.ikJoint;
        newvert.ikWeight = v1.ikWeight;
      }
    else
      {
        newvert.ikJoint = v2.ikJoint;
        newvert.ikWeight = v2.ikWeight;
      }
  }

  /** When creating a new vertex during subdivision, calculate per-vertex parameter values for the
      new vertex. */

  private static void blendParamValues(double oldValues[][][], double newValues[][][], int paramType[], int v1, int v2, int newv)
  {
    for (int i = 0; i < paramType.length; i++)
      if (paramType[i] == PER_VERTEX)
        newValues[i][0][newv] = 0.5*(oldValues[i][0][v1]+oldValues[i][0][v2]);
  }

  /** Set the per-vertex texture parameters for a newly created vertex. */

  private static void setBlendParams(double newValues[], double vert1Val[], int v2, double w1, double w2, double oldVertValues[][][], int paramType[])
  {
    for (int i = 0; i < paramType.length; i++)
      if (paramType[i] == PER_VERTEX)
        newValues[i] = w1*vert1Val[i]+w2*oldVertValues[i][0][v2];
  }

  /** Set the per-vertex texture parameters for a newly created vertex. */

  private static void setBlendParams(double newValues[], int v1, int v2, double w1, double w2, double oldVertValues[][][], int paramType[])
  {
    for (int i = 0; i < paramType.length; i++)
      if (paramType[i] == PER_VERTEX)
        newValues[i] = w1*oldVertValues[i][0][v1]+w2*oldVertValues[i][0][v2];
  }

  /** Copy the per-vertex texture parameter values for a vertex into an array. */

  private static void recordParamValues(double values[], int v, double vertValues[][][], int paramType[])
  {
    for (int i = 0; i < paramType.length; i++)
      if (paramType[i] == PER_VERTEX)
        values[i] = vertValues[i][0][v];
  }

  /** This method is used for Butterfly subdivision.  Given a face and an edge, it finds the
      other face which is across the edge from the specified one, finds the vertex of that face
      which is opposite the specified edge, and returns its position in pos.  The position of
      this "opposite vertex" can be calculated two different ways.  For smooth edges, it is
      the actual position of the vertex.  For boundary or crease edges, it is a "virtual
      vertex" created by mirroring the specified face across the specified edge.  The
      relative weights of these two are determined by smoothWeight. */

  private static void findOppositeVertex(Vertex pos, int whichFace, int whichEdge, double smoothWeight, Vertex v[], Edge e[], Face f[], double paramVal[], double oldParamValue[][][], int paramType[])
  {
    Face fc;
    Vec3 axis, delta, r;
    double dot;

    // First find the position of the actual vertex.

    if (smoothWeight > 0.0)
      {
        if (e[whichEdge].f1 == whichFace)
          fc = f[e[whichEdge].f2];
        else
          fc = f[e[whichEdge].f1];
        if (fc.e1 == whichEdge)
        {
          pos.copy(v[fc.v3]);
          recordParamValues(paramVal, fc.v3, oldParamValue, paramType);
        }
        else if (fc.e2 == whichEdge)
        {
          pos.copy(v[fc.v1]);
          recordParamValues(paramVal, fc.v1, oldParamValue, paramType);
        }
        else
        {
          pos.copy(v[fc.v2]);
          recordParamValues(paramVal, fc.v2, oldParamValue, paramType);
        }
        pos.scale(smoothWeight);
        for (int i = 0; i < paramVal.length; i++)
          paramVal[i] *= smoothWeight;
      }
    else
    {
      pos.clear();
      for (int i = 0; i < paramVal.length; i++)
        paramVal[i] = 0.0;
    }

    // Next find the position of the virtual vertex.

    if (smoothWeight < 1.0)
    {
      axis = v[e[whichEdge].v1].r.minus(v[e[whichEdge].v2].r);
      axis.normalize();
      fc = f[whichFace];
      if (fc.e1 == whichEdge)
        r = v[fc.v3].r;
      else if (fc.e2 == whichEdge)
        r = v[fc.v1].r;
      else
        r = v[fc.v2].r;
      delta = r.minus(v[e[whichEdge].v2].r);
      dot = delta.dot(axis);
      axis.scale(dot);
      delta.subtract(axis);
      pos.r.x += (1.0-smoothWeight)*(r.x+axis.x-2.0*delta.x);
      pos.r.y += (1.0-smoothWeight)*(r.y+axis.y-2.0*delta.y);
      pos.r.z += (1.0-smoothWeight)*(r.z+axis.z-2.0*delta.z);
    }
  }

  /** This method is used for Butterfly subdivision.  Given the number of edges intersecting
      an extraordinary vertex, it returns an array containing the subdivision coefficients
      to use for that vertex. */

  private static double [] getButterflyCoeff(int numEdges)
  {
    if (numEdges < BUTTERFLY_COEFF.length)
      return BUTTERFLY_COEFF[numEdges];
    double coeff[] = new double [numEdges+1], beta = 2.0*Math.PI/numEdges;
    coeff[numEdges] = 1.0;
    for (int i = 0; i < numEdges; i++)
      {
        coeff[i] = (0.25+Math.cos(beta*i)+0.5*Math.cos(2.0*beta*i))/numEdges;
        coeff[numEdges] -= coeff[i];
      }
    return coeff;
  }

  /** This method is called by the various subdivideXXX() methods to do the actual subdivision.
      The vertex, edge, and face arguments describe the mesh to be subdivided.  newvert
      contains the vertices of the new mesh.  newedge and newface are empty arrays of the
      correct length, into which the new faces and edges will be placed.  mesh is the
      TriangleMesh which the new edges and faces should belong to.  split is an array
      specifying which edges of the old mesh should be split. */

  private static void doSubdivide(TriangleMesh mesh, Vertex vertex[], Edge edge[], Face face[], boolean split[], Vertex newvert[], Edge newedge[], Face newface[], double oldParamValue[][][], double newParamValue[][][], int paramType[])
  {
    Edge tempEdge;
    Face tempFace;
    int i, j, k, n, v1, v2, v3, e1, e2, e3, newEdgeIndex[] = new int [edge.length];

    // First, subdivide edges.

    j = edge.length;
    k = vertex.length;
    for (i = 0; i < edge.length; i++)
      {
        tempEdge = edge[i];
        v1 = tempEdge.v1;
        v2 = tempEdge.v2;
        if (!split[i])
          {
            // This edge does not need to be split, so just copy it over.

            newedge[i] = mesh.new Edge(v1, v2, -1);
            newedge[i].smoothness = tempEdge.smoothness;
            if (vertex[v1].firstEdge == i)
              newvert[v1].firstEdge = i;
            if (vertex[v2].firstEdge == i)
              newvert[v2].firstEdge = i;
            newEdgeIndex[i] = i;
            continue;
          }
        newedge[i] = mesh.new Edge(v1, k, -1);
        newedge[j] = mesh.new Edge(v2, k, -1);
        newedge[i].smoothness = newedge[j].smoothness = tempEdge.smoothness;
        if (vertex[v1].firstEdge == i)
          newvert[v1].firstEdge = i;
        if (vertex[v2].firstEdge == i)
          newvert[v2].firstEdge = j;
        newvert[k].firstEdge = i;
        newEdgeIndex[i] = j++;
        k++;
      }

    // Next, subdivide faces.  For each face in the old mesh, the can be anywhere from
    // one to four faces in the new mesh, depending on how many of its edges were
    // subdivided.

    int addedFace[] = new int [4];
    k = face.length;
    for (i = 0; i < face.length; i++)
      {
        tempFace = face[i];

        // Figure out how to subdivide the face, based on which edges are subdivided.

        if (split[tempFace.e1])
          {
            if (split[tempFace.e2])
              {
                if (split[tempFace.e3])
                  {
                    n = 3;
                    v1 = tempFace.v1;  v2 = tempFace.v2;  v3 = tempFace.v3;
                    e1 = tempFace.e1;  e2 = tempFace.e2;  e3 = tempFace.e3;
                  }
                else
                  {
                    n = 2;
                    v1 = tempFace.v1;  v2 = tempFace.v2;  v3 = tempFace.v3;
                    e1 = tempFace.e1;  e2 = tempFace.e2;  e3 = tempFace.e3;
                  }
              }
            else
              {
                if (split[tempFace.e3])
                  {
                    n = 2;
                    v1 = tempFace.v3;  v2 = tempFace.v1;  v3 = tempFace.v2;
                    e1 = tempFace.e3;  e2 = tempFace.e1;  e3 = tempFace.e2;
                  }
                else
                  {
                    n = 1;
                    v1 = tempFace.v1;  v2 = tempFace.v2;  v3 = tempFace.v3;
                    e1 = tempFace.e1;  e2 = tempFace.e2;  e3 = tempFace.e3;
                  }
              }
          }
        else
          {
            if (split[tempFace.e2])
              {
                if (split[tempFace.e3])
                  {
                    n = 2;
                    v1 = tempFace.v2;  v2 = tempFace.v3;  v3 = tempFace.v1;
                    e1 = tempFace.e2;  e2 = tempFace.e3;  e3 = tempFace.e1;
                  }
                else
                  {
                    n = 1;
                    v1 = tempFace.v2;  v2 = tempFace.v3;  v3 = tempFace.v1;
                    e1 = tempFace.e2;  e2 = tempFace.e3;  e3 = tempFace.e1;
                  }
              }
            else
              {
                if (split[tempFace.e3])
                  {
                    n = 1;
                    v1 = tempFace.v3;  v2 = tempFace.v1;  v3 = tempFace.v2;
                    e1 = tempFace.e3;  e2 = tempFace.e1;  e3 = tempFace.e2;
                  }
                else
                  {
                    n = 0;
                    v1 = tempFace.v1;  v2 = tempFace.v2;  v3 = tempFace.v3;
                    e1 = tempFace.e1;  e2 = tempFace.e2;  e3 = tempFace.e3;
                  }
              }
          }

        // Now subdivide it, and create the new faces and edges.

        switch (n)
        {
          case 0:

            // No edges being split, so simply copy the face over.

            newface[i] = mesh.new Face(v1, v2, v3, e1, e2, e3);
            break;

          case 1:

            // e1 was split.

            newedge[j] = mesh.new Edge(v3, newedge[e1].v2, -1);
            if (edge[e1].v1 == v1)
              {
                newface[i] = mesh.new Face(v1, newedge[e1].v2, v3, e1, j, e3);
                newface[k] = mesh.new Face(v3, newedge[e1].v2, v2, j, newEdgeIndex[e1], e2);
              }
            else
              {
                newface[i] = mesh.new Face(v1, newedge[e1].v2, v3, newEdgeIndex[e1], j, e3);
                newface[k] = mesh.new Face(v3, newedge[e1].v2, v2, j, e1, e2);
              }
            break;

          case 2:

            // e1 and e2 were split.

            newedge[j] = mesh.new Edge(newedge[e1].v2, newedge[e2].v2, -1);
            newedge[j+1] = mesh.new Edge(v3, newedge[e1].v2, -1);
            if (edge[e1].v1 == v1)
              {
                if (edge[e2].v1 == v2)
                  {
                    newface[i] = mesh.new Face(v1, newedge[e1].v2, v3, e1, j+1, e3);
                    newface[k] = mesh.new Face(v3, newedge[e1].v2, newedge[e2].v2, j+1, j, newEdgeIndex[e2]);
                    newface[k+1] = mesh.new Face(newedge[e2].v2, newedge[e1].v2, v2, j, newEdgeIndex[e1], e2);
                  }
                else
                  {
                    newface[i] = mesh.new Face(v1, newedge[e1].v2, v3, e1, j+1, e3);
                    newface[k] = mesh.new Face(v3, newedge[e1].v2, newedge[e2].v2, j+1, j, e2);
                    newface[k+1] = mesh.new Face(newedge[e2].v2, newedge[e1].v2, v2, j, newEdgeIndex[e1], newEdgeIndex[e2]);
                  }
              }
            else
              {
                if (edge[e2].v1 == v2)
                  {
                    newface[i] = mesh.new Face(v1, newedge[e1].v2, v3, newEdgeIndex[e1], j+1, e3);
                    newface[k] = mesh.new Face(v3, newedge[e1].v2, newedge[e2].v2, j+1, j, newEdgeIndex[e2]);
                    newface[k+1] = mesh.new Face(newedge[e2].v2, newedge[e1].v2, v2, j, e1, e2);
                  }
                else
                  {
                    newface[i] = mesh.new Face(v1, newedge[e1].v2, v3, newEdgeIndex[e1], j+1, e3);
                    newface[k] = mesh.new Face(v3, newedge[e1].v2, newedge[e2].v2, j+1, j, e2);
                    newface[k+1] = mesh.new Face(newedge[e2].v2, newedge[e1].v2, v2, j, e1, newEdgeIndex[e2]);
                  }
              }
            break;

          case 3:

            // All edges being split.

            newedge[j] = mesh.new Edge(newedge[e1].v2, newedge[e2].v2, -1);
            newedge[j+1] = mesh.new Edge(newedge[e2].v2, newedge[e3].v2, -1);
            newedge[j+2] = mesh.new Edge(newedge[e3].v2, newedge[e1].v2, -1);
            if (edge[e1].v1 == v1)
              {
                if (edge[e2].v1 == v2)
                  {
                    if (edge[e3].v1 == v3)
                      {
                         newface[i] = mesh.new Face(v1, newedge[e1].v2, newedge[e3].v2, e1, j+2, newEdgeIndex[e3]);
                         newface[k] = mesh.new Face(v2, newedge[e2].v2, newedge[e1].v2, e2, j, newEdgeIndex[e1]);
                         newface[k+1] = mesh.new Face(v3, newedge[e3].v2, newedge[e2].v2, e3, j+1, newEdgeIndex[e2]);
                      }
                    else
                      {
                         newface[i] = mesh.new Face(v1, newedge[e1].v2, newedge[e3].v2, e1, j+2, e3);
                         newface[k] = mesh.new Face(v2, newedge[e2].v2, newedge[e1].v2, e2, j, newEdgeIndex[e1]);
                         newface[k+1] = mesh.new Face(v3, newedge[e3].v2, newedge[e2].v2, newEdgeIndex[e3], j+1, newEdgeIndex[e2]);
                      }
                  }
                else
                  {
                    if (edge[e3].v1 == v3)
                      {
                         newface[i] = mesh.new Face(v1, newedge[e1].v2, newedge[e3].v2, e1, j+2, newEdgeIndex[e3]);
                         newface[k] = mesh.new Face(v2, newedge[e2].v2, newedge[e1].v2, newEdgeIndex[e2], j, newEdgeIndex[e1]);
                         newface[k+1] = mesh.new Face(v3, newedge[e3].v2, newedge[e2].v2, e3, j+1, e2);
                      }
                    else
                      {
                         newface[i] = mesh.new Face(v1, newedge[e1].v2, newedge[e3].v2, e1, j+2, e3);
                         newface[k] = mesh.new Face(v2, newedge[e2].v2, newedge[e1].v2, newEdgeIndex[e2], j, newEdgeIndex[e1]);
                         newface[k+1] = mesh.new Face(v3, newedge[e3].v2, newedge[e2].v2, newEdgeIndex[e3], j+1, e2);
                      }
                  }
              }
            else
              {
                if (edge[e2].v1 == v2)
                  {
                    if (edge[e3].v1 == v3)
                      {
                         newface[i] = mesh.new Face(v1, newedge[e1].v2, newedge[e3].v2, newEdgeIndex[e1], j+2, newEdgeIndex[e3]);
                         newface[k] = mesh.new Face(v2, newedge[e2].v2, newedge[e1].v2, e2, j, e1);
                         newface[k+1] = mesh.new Face(v3, newedge[e3].v2, newedge[e2].v2, e3, j+1, newEdgeIndex[e2]);
                      }
                    else
                      {
                         newface[i] = mesh.new Face(v1, newedge[e1].v2, newedge[e3].v2, newEdgeIndex[e1], j+2, e3);
                         newface[k] = mesh.new Face(v2, newedge[e2].v2, newedge[e1].v2, e2, j, e1);
                         newface[k+1] = mesh.new Face(v3, newedge[e3].v2, newedge[e2].v2, newEdgeIndex[e3], j+1, newEdgeIndex[e2]);
                      }
                  }
                else
                  {
                    if (edge[e3].v1 == v3)
                      {
                         newface[i] = mesh.new Face(v1, newedge[e1].v2, newedge[e3].v2, newEdgeIndex[e1], j+2, newEdgeIndex[e3]);
                         newface[k] = mesh.new Face(v2, newedge[e2].v2, newedge[e1].v2, newEdgeIndex[e2], j, e1);
                         newface[k+1] = mesh.new Face(v3, newedge[e3].v2, newedge[e2].v2, e3, j+1, e2);
                      }
                    else
                      {
                         newface[i] = mesh.new Face(v1, newedge[e1].v2, newedge[e3].v2, newEdgeIndex[e1], j+2, e3);
                         newface[k] = mesh.new Face(v2, newedge[e2].v2, newedge[e1].v2, newEdgeIndex[e2], j, e1);
                         newface[k+1] = mesh.new Face(v3, newedge[e3].v2, newedge[e2].v2, newEdgeIndex[e3], j+1, e2);
                      }
                  }
              }
            newface[k+2] = mesh.new Face(newedge[e1].v2, newedge[e2].v2, newedge[e3].v2, j, j+1, j+2);
        }

        // Copy over per-face and per-face/per-vertex parameter values.

        int numAddedFaces = n+1;
        addedFace[0] = i;
        for (int m = 0; m < n; m++)
          addedFace[m+1] = k+m;
        for (int p = 0; p < paramType.length; p++)
        {
          if (paramType[p] == PER_FACE)
            for (int m = 0; m < numAddedFaces; m++)
              newParamValue[p][0][addedFace[m]] = oldParamValue[p][0][i];
          else if (paramType[p] == PER_FACE_VERTEX)
          {
            int vertInd[] = new int [] {v1, v2, v3, -1, -1, -1};
            if (n > 0)
              vertInd[3] = newedge[e1].v2;
            if (n > 1)
              vertInd[4] = newedge[e2].v2;
            if (n > 2)
              vertInd[5] = newedge[e3].v2;
            double vertVal[];
            if (v1 == tempFace.v1)
              vertVal = new double [] {oldParamValue[p][0][i], oldParamValue[p][1][i], oldParamValue[p][2][i],
                  0.5*(oldParamValue[p][0][i]+oldParamValue[p][1][i]),
                  0.5*(oldParamValue[p][1][i]+oldParamValue[p][2][i]),
                  0.5*(oldParamValue[p][2][i]+oldParamValue[p][0][i])};
            else if (v1 == tempFace.v2)
              vertVal = new double [] {oldParamValue[p][1][i], oldParamValue[p][2][i], oldParamValue[p][0][i],
                  0.5*(oldParamValue[p][1][i]+oldParamValue[p][2][i]),
                  0.5*(oldParamValue[p][2][i]+oldParamValue[p][0][i]),
                  0.5*(oldParamValue[p][0][i]+oldParamValue[p][1][i])};
            else
              vertVal = new double [] {oldParamValue[p][2][i], oldParamValue[p][0][i], oldParamValue[p][1][i],
                  0.5*(oldParamValue[p][2][i]+oldParamValue[p][0][i]),
                  0.5*(oldParamValue[p][0][i]+oldParamValue[p][1][i]),
                  0.5*(oldParamValue[p][1][i]+oldParamValue[p][2][i])};
            for (int m = 0; m < numAddedFaces; m++)
            {
              Face fc = newface[addedFace[m]];
              for (int q = 0; q < 6; q++)
              {
                if (fc.v1 == vertInd[q])
                  newParamValue[p][0][addedFace[m]] = vertVal[q];
                else if (fc.v2 == vertInd[q])
                  newParamValue[p][1][addedFace[m]] = vertVal[q];
                else if (fc.v3 == vertInd[q])
                  newParamValue[p][2][addedFace[m]] = vertVal[q];
              }
            }
          }
        }
        j += n;
        k += n;
      }

    // Record which faces are adjacent to each edge.

    for (i = 0; i < newface.length; i++)
      {
        tempFace = newface[i];
        if (newedge[tempFace.e1].f1 == -1)
          newedge[tempFace.e1].f1 = i;
        else
          newedge[tempFace.e1].f2 = i;
        if (newedge[tempFace.e2].f1 == -1)
          newedge[tempFace.e2].f1 = i;
        else
          newedge[tempFace.e2].f2 = i;
        if (newedge[tempFace.e3].f1 == -1)
          newedge[tempFace.e3].f1 = i;
        else
          newedge[tempFace.e3].f2 = i;
      }

    // Count the number of edges intersecting each vertex.

    for (i = 0; i < newvert.length; i++)
      newvert[i].edges = 0;
    for (i = 0; i < newedge.length; i++)
      {
        newvert[newedge[i].v1].edges++;
        newvert[newedge[i].v2].edges++;
      }
  }

 

  /** If necessary, reorder the points in each face so that the normals will be properly oriented. */

  public void makeRightSideOut()
  {
    Vec3 norm[] = getNormals();
    int maxLenVertex = 0;
    double maxLength = 0.0;
    for (int i = 0; i < vertex.length; i++)
    {
      double length = vertex[i].r.length();
      if (length > maxLength)
      {
        maxLenVertex = i;
        maxLength = length;
      }
    }
    if (vertex[maxLenVertex].r.dot(norm[maxLenVertex]) < 0.0)
      reverseNormals();
  }

  /** Reorder the vertices in each face, so as to reverse all of the normal vectors. */

  public void reverseNormals()
  {
    int i, temp;

    for (i = 0; i < face.length; i++)
      {
        temp = face[i].v2;
        face[i].v2 = face[i].v3;
        face[i].v3 = temp;
        temp = face[i].e1;
        face[i].e1 = face[i].e3;
        face[i].e3 = temp;
      }
  }

  /** Get an array of normal vectors.  This calculates a single normal for each vertex,
      ignoring smoothness values. */


  public Vec3 [] getNormals()
  {
    Vec3 faceNorm, norm[] = new Vec3 [vertex.length];

    // Calculate a normal for each face, and average the face normals for each vertex.

    for (int i = 0; i < norm.length; i++)
      norm[i] = new Vec3();
    for (int i = 0; i < face.length; i++)
      {
        Vec3 edge1 = vertex[face[i].v2].r.minus(vertex[face[i].v1].r);
        Vec3 edge2 = vertex[face[i].v3].r.minus(vertex[face[i].v1].r);
        Vec3 edge3 = vertex[face[i].v3].r.minus(vertex[face[i].v2].r);
        edge1.normalize();
        edge2.normalize();
        edge3.normalize();
        faceNorm = edge1.cross(edge2);
        double length = faceNorm.length();
        if (length == 0.0)
          continue;
        faceNorm.scale(1.0/length);
        double dot1 = edge1.dot(edge2);
        double dot2 = -edge1.dot(edge3);
        double dot3 = edge2.dot(edge3);
        if (dot1 < -1.0)
          dot1 = -1.0; // This can occassionally happen due to roundoff error
        if (dot1 > 1.0)
          dot1 = 1.0;
        if (dot2 < -1.0)
          dot2 = -1.0;
        if (dot2 > 1.0)
          dot2 = 1.0;
        if (dot3 < -1.0)
          dot3 = -1.0;
        if (dot3 > 1.0)
          dot3 = 1.0;
        norm[face[i].v1].add(faceNorm.times(Math.acos(dot1)));
        norm[face[i].v2].add(faceNorm.times(Math.acos(dot2)));
        norm[face[i].v3].add(faceNorm.times(Math.acos(dot3)));
      }
    for (int i = 0; i < norm.length; i++)
      norm[i].normalize();
    return norm;
  }

  public int getFaceCount()
  {
    return face.length;
  }
  
  public int getFaceVertexCount(int faceIndex)
  {
    return 3;
  }

  public int getFaceVertexIndex(int faceIndex, int vertexIndex)
  {
    Face f = face[faceIndex];
    if (vertexIndex == 0)
      return f.v1;
    if (vertexIndex == 1)
      return f.v2;
    return f.v3;
  }

  /** Return a new mesh which is an "optimized" version of the input mesh.  This is done by rearranging edges
      to eliminate very small angles, or vertices where many edges come together.  The resulting mesh will
      generally produce a better looking surface after smoothing is applied to it. */

  public static TriangleMesh optimizeMesh(TriangleMesh mesh)
  {
    Face face[] = mesh.face;
    Edge edge[] = mesh.edge;
    Vertex vertex[] = mesh.vertex;
    boolean candidate[] = new boolean [edge.length];
    TriangleMesh newmesh = null;

    for (int i = 0; i < edge.length; i++)
      candidate[i] = true;
    while (true)
      {
        Vec3 faceNorm[] = new Vec3 [face.length];
        boolean onBoundary[] = new boolean [vertex.length];
        int numEdges[] = new int [vertex.length];

        // Initialize the various arrays, and determine which edges are really candidates for optimization.

        for (int i = 0; i < face.length; i++)
          {
            Face f = face[i];
            if (candidate[f.e1] || candidate[f.e2] || candidate[f.e3])
              {
                Vec3 v1 = vertex[f.v1].r;
                Vec3 v2 = vertex[f.v2].r;
                Vec3 v3 = vertex[f.v3].r;
                Vec3 d1 = v2.minus(v1);
                Vec3 d2 = v3.minus(v1);
                faceNorm[i] = d1.cross(d2);
                double length = faceNorm[i].length();
                if (length > 0.0)
                  faceNorm[i].scale(1.0/length);
                else if (!v1.equals(v2) && !v1.equals(v3) && !v2.equals(v3))
                  faceNorm[i] = null;
              }
          }
        for (int i = 0; i < edge.length; i++)
          {
            Edge e = edge[i];
            numEdges[e.v1]++;
            numEdges[e.v2]++;
            if (e.f2 == -1)
              {
                onBoundary[e.v1] = onBoundary[e.v2] = true;
                candidate[i] = false;
              }
            else if (candidate[i] && faceNorm[e.f1] != null && faceNorm[e.f2] != null && faceNorm[e.f1].dot(faceNorm[e.f2]) < 0.99)
              candidate[i] = false;
          }

        // For each candidate edge, find the list of vertices and angles involved in swapping it.  The vertices
        // are ordered as follows:

        //              <-
        //              /\ 2
        //             /f1\
        //            /    \
        //          0 ------ 1
        //            \    /
        //             \f2/
        //              \/ 3
        //              ->

        int swapVert[][] = new int [edge.length][];
        double minAngle[][] = new double [edge.length][];
        for (int i = 0; i < edge.length; i++)
          if (candidate[i])
            {
              // First find the vertices.

              swapVert[i] = new int [4];
              Edge e = edge[i];
              Face f1 = face[e.f1], f2 = face[e.f2];
              if ((e.v1 == f1.v1 && e.v2 == f1.v2) || (e.v1 == f1.v2 && e.v2 == f1.v3) || (e.v1 == f1.v3 && e.v2 == f1.v1))
                {
                  swapVert[i][0] = e.v1;
                  swapVert[i][1] = e.v2;
                }
              else
                {
                  swapVert[i][0] = e.v2;
                  swapVert[i][1] = e.v1;
                }
              if (e.v1 != f1.v1 && e.v2 != f1.v1)
                swapVert[i][2] = f1.v1;
              else if (e.v1 != f1.v2 && e.v2 != f1.v2)
                swapVert[i][2] = f1.v2;
              else
                swapVert[i][2] = f1.v3;
              if (e.v1 != f2.v1 && e.v2 != f2.v1)
                swapVert[i][3] = f2.v1;
              else if (e.v1 != f2.v2 && e.v2 != f2.v2)
                swapVert[i][3] = f2.v2;
              else
                swapVert[i][3] = f2.v3;

              // Now calculate the angles.

              minAngle[i] = new double [4];
              Vec3 d1 = vertex[swapVert[i][1]].r.minus(vertex[swapVert[i][0]].r);
              Vec3 d2 = vertex[swapVert[i][2]].r.minus(vertex[swapVert[i][0]].r);
              Vec3 d3 = vertex[swapVert[i][3]].r.minus(vertex[swapVert[i][0]].r);
              Vec3 d4 = vertex[swapVert[i][2]].r.minus(vertex[swapVert[i][1]].r);
              Vec3 d5 = vertex[swapVert[i][3]].r.minus(vertex[swapVert[i][1]].r);
              Vec3 d6 = vertex[swapVert[i][3]].r.minus(vertex[swapVert[i][2]].r);
              d1.normalize();
              d2.normalize();
              d3.normalize();
              d4.normalize();
              d5.normalize();
              d6.normalize();
              double a1, a2;
              a1 = Math.acos(d1.dot(d2));
              a2 = Math.acos(d1.dot(d3));
              if (a1+a2 >= Math.PI)
                {
                  candidate[i] = false;
                  continue;
                }
              minAngle[i][0] = (a1 < a2 ? a1 : a2);
              a1 = Math.acos(-d1.dot(d4));
              a2 = Math.acos(-d1.dot(d5));
              if (a1+a2 >= Math.PI)
                {
                  candidate[i] = false;
                  continue;
                }
              minAngle[i][1] = (a1 < a2 ? a1 : a2);
              a1 = Math.acos(-d6.dot(d2));
              a2 = Math.acos(-d6.dot(d4));
              minAngle[i][2] = (a1 < a2 ? a1 : a2);
              a1 = Math.acos(d6.dot(d3));
              a2 = Math.acos(d6.dot(d5));
              minAngle[i][3] = (a1 < a2 ? a1 : a2);
            }

        // Calculate scores for each candidate edge, and decide which ones to swap.

        double score[] = new double [edge.length];
        boolean swap[] = new boolean [edge.length];
        for (int i = 0; i < score.length; i++)
          if (candidate[i])
            score[i] = calcSwapScore(minAngle[i], swapVert[i], numEdges, onBoundary);
        while (true)
          {
            int best = -1;
            double maxScore = 0.01;
            for (int i = 0; i < candidate.length; i++)
              if (candidate[i] && score[i] > maxScore)
                {
                  best = i;
                  maxScore = score[i];
                }
            if (best == -1)
              break;

            // Mark the edge to be swapped.  Remove it and every other edge that shares a face with it
            // from the candidate list.

            swap[best] = true;
            Edge e = edge[best];
            Face f = face[e.f1];
            candidate[f.e1] = candidate[f.e2] = candidate[f.e3] = false;
            f = face[e.f2];
            candidate[f.e1] = candidate[f.e2] = candidate[f.e3] = false;

            // Update the numEdges array, and recalculate scores.

            numEdges[swapVert[best][0]]--;
            numEdges[swapVert[best][1]]--;
            numEdges[swapVert[best][2]]++;
            numEdges[swapVert[best][3]]++;
            for (int i = 0; i < 4; i++)
              {
                int vertEdges[] = vertex[swapVert[best][i]].getEdges();
                for (int j = 0; j < vertEdges.length; j++)
                  if (candidate[vertEdges[j]])
                    score[vertEdges[j]] = calcSwapScore(minAngle[vertEdges[j]], swapVert[vertEdges[j]], numEdges, onBoundary);
              }
          }

        // We now know which edges we want to swap.  Create the new mesh.

        int newface[][] = new int [face.length][];
        int next = 0;
        for (int i = 0; i < face.length; i++)
          {
            Face f = face[i];
            if (!swap[f.e1] && !swap[f.e2] && !swap[f.e3])
              newface[next++] = new int [] {f.v1, f.v2, f.v3};
          }
        int firstSplit = next;
        for (int i = 0; i < edge.length; i++)
          if (swap[i])
            {
              newface[next++] = new int [] {swapVert[i][2], swapVert[i][0], swapVert[i][3]};
              newface[next++] = new int [] {swapVert[i][2], swapVert[i][3], swapVert[i][1]};
            }
        newmesh = new TriangleMesh(vertex, newface);

        // Copy over edge smoothness values.

        Vertex newvert[] = (Vertex []) newmesh.getVertices();
        Edge newedge[] = newmesh.getEdges();
        for (int i = 0; i < edge.length; i++)
          if (!swap[i] && edge[i].smoothness != 1.0)
            {
              int edges2[] = newvert[edge[i].v1].getEdges();
              Edge e1 = edge[i];
              for (int k = 0; k < edges2.length; k++)
                {
                  Edge e2 = newedge[edges2[k]];
                  if ((e1.v1 == e2.v1 && e1.v2 == e2.v2) || (e1.v1 == e2.v2 && e1.v2 == e2.v1))
                    {
                      e2.smoothness = e1.smoothness;
                      break;
                    }
                }
            }

        // Determine which edges are candidates for the next iteration.

        if (firstSplit == next)
          break;
        vertex = newvert;
        edge = newedge;
        face = newmesh.getFaces();
        candidate = new boolean [edge.length];
        for (int i = firstSplit; i < face.length; i++)
          {
            Face f = face[i];
            candidate[f.e1] = candidate[f.e2] = candidate[f.e3] = true;
          }
      }

    // Copy over other mesh properties.
    return newmesh;
  }

  /** This is a utility routine used by optimizeMesh().  It calculates the score for swapping a
      particular edge. */

  private static double calcSwapScore(double minAngle[], int vert[], int numEdges[], boolean onBoundary[])
  {
    double s[] = new double [4];
    for (int i = 0; i < 4; i++)
      {
        int v = vert[i];
        int ideal = (onBoundary[v] ? 4 : 6);
        if (i > 1)
          ideal--;
        s[i] = (numEdges[v] > ideal ? minAngle[i]/(numEdges[v]-ideal+1.5) : minAngle[i]);
      }
    double currentScore = (s[0] < s[1] ? s[0] : s[1]);
    double swapScore = (s[2] < s[3] ? s[2] : s[3]);
    return swapScore-currentScore;
  }

  /** Automatically select smoothness values for all edges in the mesh.  This is done based on the
      angle between the two faces sharing the edge.  If it is greater than the specified cutoff, the
      edge smoothness is set to 0.  Otherwise, it is set to 1.
      @param angle      the cutoff angle, in radians
   */

  public void autosmoothMeshEdges(double angle)
  {
    double cutoff = Math.cos(angle);
    for (int i = 0; i < edge.length; i++)
    {
      if (edge[i].f2 == -1)
      {
        edge[i].smoothness = 1.0f;
        continue;
      }
      Face f1 = face[edge[i].f1];
      Face f2 = face[edge[i].f2];
      Vec3 norm1 = vertex[f1.v1].r.minus(vertex[f1.v2].r).cross(vertex[f1.v1].r.minus(vertex[f1.v3].r));
      Vec3 norm2 = vertex[f2.v1].r.minus(vertex[f2.v2].r).cross(vertex[f2.v1].r.minus(vertex[f2.v3].r));
      norm1.normalize();
      norm2.normalize();
      if (norm1.dot(norm2) < cutoff)
        edge[i].smoothness = 0.0f;
      else
        edge[i].smoothness = 1.0f;
    }
  }
}
