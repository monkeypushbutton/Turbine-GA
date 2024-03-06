#!/usr/bin/env python3

"""
author:  Nils Olofsson
email:   nils.olovsson@gmail.com
website: http://nilsolovsson.se
date:    2021-02-28
license: MIT
url: https://bitbucket.org/nils_olovsson/ear_clipping_triangulation/raw/d5c1927262c1a9b2c6a1b7b257bb18689e6e149e/ear_clipping_triangulation.py

"""

import numpy as np
from double_linked_list import DoubleLinkedList

# ==============================================================================
#
# Geometry
#
# ==============================================================================

def angleCCW(a, b):
    """
        Counter clock wise angle (radians) from normalized 2D vectors a to b
    """
    dot = a[0]*b[0] + a[1]*b[1]
    det = a[0]*b[1] - a[1]*b[0]
    angle = np.arctan2(det, dot)
    if angle<0.0 :
        angle = 2.0*np.pi + angle
    return angle

def isConvex(vertex_prev, vertex, vertex_next):
    """
        Determine if vertex lies on the convex hull of the polygon.
    """
    a = vertex_prev - vertex
    b = vertex_next - vertex
    internal_angle = angleCCW(b, a)
    return internal_angle <= np.pi

def insideTriangle(a, b, c, p):
    """
        Determine if a vertex p is inside (or "on") a triangle made of the
        points a->b->c
        http://blackpawn.com/texts/pointinpoly/
    """

    #Compute vectors
    v0 = c - a
    v1 = b - a
    v2 = p - a

    # Compute dot products
    dot00 = np.dot(v0, v0)
    dot01 = np.dot(v0, v1)
    dot02 = np.dot(v0, v2)
    dot11 = np.dot(v1, v1)
    dot12 = np.dot(v1, v2)

    # Compute barycentric coordinates
    denom = dot00*dot11 - dot01*dot01
    if abs(denom) < 1e-20:
        return True
    invDenom = 1.0 / denom
    u = (dot11*dot02 - dot01*dot12) * invDenom
    v = (dot00*dot12 - dot01*dot02) * invDenom

    # Check if point is in triangle
    return (u >= 0) and (v >= 0) and (u + v < 1)

def triangulate(vertices, max_iterations=0):
    """
        Triangulation of a polygon in 2D.
        Assumption that the polygon is simple, i.e has no holes, is closed and
        has no crossings and also that it the vertex order is counter clockwise.
        https://geometrictools.com/Documentation/TriangulationByEarClipping.pdf
    """

    n = vertices.shape[0]
    indices = np.zeros([n-2, 3], dtype=int)

    #print('shape: {}x{}'.format(n,m))

    vertlist = DoubleLinkedList()
    for i in range(0, n):
        vertlist.append(i)

    index_counter = 0
    it_counter = 0

    # Simplest possible algorithm. Create list of indexes.
    # Find first ear vertex. Create triangle. Remove vertex from list
    # Do this while number of vertices > 2.
    node = vertlist.first
    #while vertlist.size > 2 and it_counter < 10:
    while vertlist.size > 2 and (max_iterations<=0 or max_iterations>index_counter):
        #print(it_counter)
        #print('vertlist.size: {}'.format(vertlist.size))
        i = node.prev.data
        j = node.data
        k = node.next.data

        vert_prev = vertices[i,:]
        vert_crnt = vertices[j,:]
        vert_next = vertices[k,:]

        is_convex = isConvex(vert_prev, vert_crnt, vert_next)
        is_ear = True
        if is_convex:
            test_node = node.next.next
            while test_node!=node.prev and is_ear:
                vert = vertices[test_node.data,:]
                is_ear = not insideTriangle(vert_prev, vert_crnt, vert_next, vert)
                test_node = test_node.next
        else:
            is_ear = False

        if is_ear:
            indices[index_counter, :] = np.array([i, j, k], dtype=int)
            index_counter += 1
            vertlist.remove(node.data)
        it_counter += 1
        node = node.next
    indices = indices[0:index_counter, :]
    return indices
