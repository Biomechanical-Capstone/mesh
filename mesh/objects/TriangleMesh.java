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
