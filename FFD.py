import numpy as np
import matplotlib.pyplot as plt
# NOTE: This library is https://github.com/taxpon/pymesh
# NOT https://github.com/PyMesh/PyMesh
import pymesh
from mpl_toolkits import mplot3d
from scipy.interpolate import RectBivariateSpline


def get_base_curve(segments):
    """Define the baseline unit circle
    
    Parameters
    ----------
    segments : int
        Number of segments curve will be discretised into.
    """
    # 
    theta = np.linspace(0, 2*np.pi, segments)
    x_circle = np.cos(theta)
    y_circle = np.sin(theta)
    return x_circle, y_circle

def get_ruler(ticks):
    """Get an array of points evenly spaced between -1 and 1 with length equal to input ticks
    
    Parameters
    ----------
    ticks : int
        Number of points in returned array.
    """
    return np.linspace(-1, 1, ticks)

def get_control_points(x_ruler, y_ruler):
    """Define the control point grid
    
    Parameters
    ----------
    x_ruler : array_like
        Define spacing of points on x axis
    y_ruler : array_like
        Define spacing of points on y axis

    Returns
    -------
    ndarray
        Mesh points with shape (X*Y, 2) where X = len(x_ruler)
        and Y = len(y_ruler)

    """
    X, Y = np.meshgrid(x_ruler, y_ruler)
    return np.column_stack((X.ravel(), Y.ravel()))

def get_curve_spline_ffd(x_ruler, y_ruler, displacements, x_base_curve, y_base_curve):
    """Perform a B-Spline free form deformation (FFD) of a given curve

    Parameters
    ----------
    x_ruler : array_like
        Define spacing of points on x axis
    y_ruler : array_like
        Define spacing of points on y axis
    displacements : ndarray
        shape (X*Y, 2) Define how control point grid is manipulated to perform deformation
    x_base_curve : array_like
    y_base_curve : array_like

    Returns
    -------
    (array_like, array_like)
        arrays of X, Y co-ordinates of deformed curve
    """
    spline_x = RectBivariateSpline(x_ruler, y_ruler, displacements[:, 0].reshape((len(x_ruler), len(y_ruler))))
    spline_y = RectBivariateSpline(x_ruler, y_ruler, displacements[:, 1].reshape((len(x_ruler), len(y_ruler))))
    x_deformed = x_base_curve + spline_x.ev(x_base_curve, y_base_curve)
    y_deformed = y_base_curve + spline_y.ev(x_base_curve, y_base_curve)
    return x_deformed, y_deformed

def generate_planform_curve(displacements, nx = 5, ny = 5, segments = 100, show_curve = False):
    """
    Create a closed 2d polygon by deforming a circle using a set of displacement information

    Parameters
    ----------
    displacements : ndarray
        shape (nx*ny, 2) Define how control point grid is manipulated to perform deformation
    nx : int
        Number of control points to use to deform curve in x direction
    ny : int
        Number of control points to use to deform curve in y direction
    segments : int
        How many edges resulting polygon should have
    show_curve : boolean
        Should curve be displayed to user in a modal dialog (for debugging)

    Returns
    -------
    ndarray
        Shape (segments, 2) representing points of vertices of resulting polygon
    """
    x_base_curve, y_base_curve = get_base_curve(segments)
    x_ruler = get_ruler(nx)
    y_ruler = get_ruler(ny)
    x_deformed, y_deformed = get_curve_spline_ffd(x_ruler, y_ruler, displacements, x_base_curve, y_base_curve)
    if(show_curve):
        control_points = get_control_points(x_ruler, y_ruler)
        plot_planform_data(x_base_curve, y_base_curve, x_deformed, y_deformed, control_points, displacements)
    return np.column_stack((x_deformed, y_deformed))

def calculate_centroid(triangle, vertices):
    """
    Calculate the centroid of a triangle

    Parameters
    ----------
    triangle : arraylike
        Indices of triangle verices. Expected to be length 3.
    vertices : ndarray
        Vertices that may be referred to by triangle indices. Expected shape (N, 2)

    Returns
    -------
    (float, float)
        X, Y coordinates of centroid of triangle
    """
    return ((vertices[triangle[0]][0] + vertices[triangle[1]][0] + vertices[triangle[2]][0]) / 3, (vertices[triangle[0]][1] + vertices[triangle[1]][1] + vertices[triangle[2]][1]) / 3, )

# Checking if a point is inside a polygon
def point_in_polygon(x, y, polygon):
    """
    Check if a point is inside a given polygon

    Uses ray casting technique. Implementation is from https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/

    Parameters
    ----------
    x : float
        X coordinate of point to check
    y : float
        Y coordinate of point to check
    polygon : ndarray
        Vertices of polygon to check. Shape (N, 2)

    Returns
    -------
    bool
        True if point is inside polygon, False otherwise
    """
    num_vertices = len(polygon)
    inside = False
 
    # Store the first point in the polygon and initialize the second point
    p1 = polygon[0]
 
    # Loop through each edge in the polygon
    for i in range(1, num_vertices + 1):
        # Get the next point in the polygon
        p2 = polygon[i % num_vertices]
        p1x, p1y = p1
        p2x, p2y = p2

        # Check if the point is above the minimum y coordinate of the edge
        if y > min(p1y, p2y):
            # Check if the point is below the maximum y coordinate of the edge
            if y <= max(p1y, p2y):
                # Check if the point is to the left of the maximum x coordinate of the edge
                if x <= max(p1x, p2x):
                    # Calculate the x-intersection of the line connecting the point to the edge
                    x_intersection = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
 
                    # Check if the point is on the same line as the edge or to the left of the x-intersection
                    if p1x == p2x or x <= x_intersection:
                        # Flip the inside flag
                        inside = not inside
 
        # Store the current point as the first point for the next iteration
        p1 = p2
 
    # Return the value of the inside flag
    return inside

def triangulate_poly_pymesh(vertices, flip = False):
    """
    Calculate a triangulation of the given polygon

    Uses the pymesh library to calculate a Delaunay triangulation. This is then
    improved by removing triangles that were exterior to the polygon.

    Parameters
    ----------
    vertices : ndarray
        Array of points in polygon in 2-D or 3-D. If 3-D points must be coplanar.
    flip : boolean
        If true, invert the normal vector of generated triangles

    Returns
    -------
    pymesh.Mesh
        A mesh object representing a triangulation of the interior of the given polygon.
    """
    tri = pymesh.triangle()
    tri.points = vertices
    tri.split_boundary = False
    tri.keep_convex_hull = False
    tri.verbosity = 1
    tri.run()
    top_mesh_vertices = tri.mesh.vertices
    top_mesh_faces = tri.mesh.faces

    face_removals = []
    for face_index in range(len(top_mesh_faces)):
        f = top_mesh_faces[face_index]
        centroid = calculate_centroid(f, top_mesh_vertices)
        if point_in_polygon(centroid[0], centroid[1], planform) == False:
            face_removals.append(face_index)

    top_mesh_faces = np.delete(top_mesh_faces, face_removals, axis=0)

    if flip:
        permutation = [0, 2, 1]
        idx = np.empty_like(permutation)
        idx[permutation] = np.arange(len(permutation))
        top_mesh_faces[:] = top_mesh_faces[:, idx]

    return pymesh.meshio.form_mesh(top_mesh_vertices, top_mesh_faces)

def generate_blade_mesh(planform, depth, show_mesh = False):
    """
    Generate a 3D blade model in STL format by extruding a 2D curve into a prism along the z-axis.

    Parameters
    ----------
    planform : numpy.ndarray
        Array of tuples (x,y) specifying vertices of the planform curve.
    depth : float 
        Z depth of resulting shape

    Returns
    -------
    pymesh.Mesh
        Mesh object representing a candidate turbine blade
    """
    segments = len(planform)
    bottom_vertices = np.column_stack((planform[:,0], planform[:,1], np.zeros(segments)))
    top_vertices    = np.column_stack((planform[:,0][::-1], planform[:,1][::-1], np.full(segments, depth)))
    vertices = np.concatenate((bottom_vertices, top_vertices))
   # Generate the indices for the rectangular walls
    wall_faces = np.empty((segments * 2, 3), int)
    for i in range(segments):
        # Overall quad face is e.g. (0, 1, 6, 7) assuming segments = 4, and counting nodes clockwise
        # Split into two triangles making sure normals match and wrapping 0 = 4 and 4 = 8 vertices
        # 0, 1, 3
        wall_faces[i * 2] = [i, (i + 1) % segments, (2 * segments) - i - 1]
        # 3, 4, 1
        wall_faces[(i * 2) + 1] = [(2 * segments) - i - 1, (i + 1) % segments, (((2 * segments) - i - 2) % segments) + segments]
    wall_mesh = pymesh.meshio.form_mesh(vertices, wall_faces)
    top_mesh = triangulate_poly_pymesh(top_vertices)
    bottom_mesh = triangulate_poly_pymesh(bottom_vertices, True)
    blade_mesh = pymesh.merge_meshes([wall_mesh, top_mesh, bottom_mesh])

    pymesh.remove_degenerated_triangles(blade_mesh)
    pymesh.remove_duplicated_vertices(blade_mesh)
    pymesh.remove_isolated_vertices(blade_mesh)
    pymesh.remove_duplicated_vertices(blade_mesh)

    if(show_mesh):
        plot_mesh(blade_mesh)

    return blade_mesh

def plot_planform_data(x_base_curve, y_base_curve, x_deformed, y_deformed, control_points, displacements):
    """
    Plot data about planform curve in a modal dialog for debugging
    """
    plt.figure(figsize=(8, 6))
    plt.plot(x_base_curve, y_base_curve, 'b-', label='Baseline')
    plt.plot(x_deformed, y_deformed, 'r-', label='Deformed')
    plt.scatter(np.add(control_points[:, 0], displacements[:, 0]), np.add(control_points[:, 1], displacements[:, 1]), c='k', marker='x', label='Control Points')
    plt.axis('equal')
    plt.legend()
    plt.show()
    print(x_base_curve)
    print(y_base_curve)

def plot_mesh(mesh):
    """
    Plot data about 3D mesh in a modal dialog for debugging
    """
    figure = plt.figure()
    axes = mplot3d.Axes3D(figure)

    # Load the STL files and add the vectors to the plot
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(mesh.vectors))

    # Auto scale to the mesh size
    scale = mesh.points.flatten()
    axes.auto_scale_xyz(scale, scale, scale)

    # Show the plot to the screen
    plt.show()

# Driver code.
if __name__ == '__main__':
    # Define the initial displacements of the control points
    nx, ny = 5, 5
    displacements = np.zeros((nx*ny, 2))
    displacements[0] = 1, 0.5
    displacements[10] = -0.2, 0.6
    displacements[20] = 0.5, 1
    displacements[14] = -0.5, -1
    planform = generate_planform_curve(displacements, nx, ny)
    blade_mesh = generate_blade_mesh(planform, 0.3)

    # Save the STL file
    pymesh.save_mesh("blade.stl", blade_mesh, ascii=True)

