import numpy as np
import matplotlib.pyplot as plt
import ear_clipping_triangulation as ear_trian

# NOTE: This library is https://github.com/taxpon/pymesh
# NOT https://github.com/PyMesh/PyMesh
import pymesh

from mpl_toolkits import mplot3d

from stl import mesh
from scipy.interpolate import RectBivariateSpline


def get_base_curve(segments):
    # Define the baseline unit circle
    theta = np.linspace(0, 2*np.pi, segments)
    x_circle = np.cos(theta)
    y_circle = np.sin(theta)
    return x_circle, y_circle

def get_ruler(ticks):
    # Get an array of points evenly spaced between -1 and 1 with length equal to input ticks
    return np.linspace(-1, 1, ticks)

def get_control_points(x_ruler, y_ruler):
    # Define the control point grid
    X, Y = np.meshgrid(x_ruler, y_ruler)
    return np.column_stack((X.ravel(), Y.ravel()))

def get_curve_spline_ffd(x_ruler, y_ruler, displacements, x_base_curve, y_base_curve):
    # Perform B-spline FFD
    spline_x = RectBivariateSpline(x_ruler, y_ruler, displacements[:, 0].reshape((len(x_ruler), len(y_ruler))))
    spline_y = RectBivariateSpline(x_ruler, y_ruler, displacements[:, 1].reshape((len(x_ruler), len(y_ruler))))
    x_deformed = x_base_curve + spline_x.ev(x_base_curve, y_base_curve)
    y_deformed = y_base_curve + spline_y.ev(x_base_curve, y_base_curve)
    return x_deformed, y_deformed

def generate_planform_curve(displacements, nx, ny, segments = 100, show_curve = False):
    # Generate a closed curve using some displacements of a control point matrix over the unit circle
    x_base_curve, y_base_curve = get_base_curve(segments)
    x_ruler = get_ruler(nx)
    y_ruler = get_ruler(ny)
    x_deformed, y_deformed = get_curve_spline_ffd(x_ruler, y_ruler, displacements, x_base_curve, y_base_curve)
    if(show_curve):
        control_points = get_control_points(x_ruler, y_ruler)
        plot_planform_data(x_base_curve, y_base_curve, x_deformed, y_deformed, control_points, displacements)
    return np.column_stack((x_deformed, y_deformed))

def plot_planform_data(x_base_curve, y_base_curve, x_deformed, y_deformed, control_points, displacements):
    # Plot the baseline and deformed curves
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
    figure = plt.figure()
    axes = mplot3d.Axes3D(figure)

    # Load the STL files and add the vectors to the plot
    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(mesh.vectors))

    # Auto scale to the mesh size
    scale = mesh.points.flatten()
    axes.auto_scale_xyz(scale, scale, scale)

    # Show the plot to the screen
    plt.show()

def generate_blade_mesh(planform, depth, show_mesh = False):
    """
        Generate a 3D blade model in STL format by extruding a 2D curve along the z-axis.

    Args:
        planform (numpy.ndarray): array of tuples (x,y) specifying vertices of the planform curve.
        depth (float): z depth of blade

    Returns:
        Mesh object representing a potential turbine blade
    """
    segments = len(planform)
    bottom_vertices = np.column_stack((planform[:,0], planform[:,1], np.zeros(segments)))
    top_vertices    = np.column_stack((planform[:,0][::-1], planform[:,1][::-1], np.full(segments, depth)))
    vertices = np.concatenate((bottom_vertices, top_vertices))
   # Generate the indices for the rectangular walls
    wall_faces = np.empty((segments, 4), int)
    for i in range(segments):
        wall_faces[i] = [i, (i + 1) % segments, i + segments, ((i + 1) % segments) + segments]

    wall_mesh = pymesh.meshio.form_mesh(vertices, wall_faces)

    tri = pymesh.triangle()
    tri.points = top_vertices
    tri.max_area = 0.1
    tri.split_boundary = False
    tri.verbosity = 1
    tri.run()
    top_mesh = tri.mesh

    blade_mesh = pymesh.merge_meshes([wall_mesh, top_mesh])

    #top_faces = ear_trian.triangulate(planform)
    # bottom_faces = top_faces + segments
    # Flip the normals of the bottom face triangles.
    # bottom_faces = bottom_faces[:,[0, 2, 1]]
    # faces = np.concatenate((wall_faces, top_faces, bottom_faces))

    #blade_mesh =

    # blade = mesh.Mesh(np.zeros(faces.shape[0], dtype=mesh.Mesh.dtype))
    # for i, f in enumerate(faces):
    #     for j in range(3):
    #         blade.vectors[i][j] = vertices[f[j],:]


    if(show_mesh):
        plot_mesh(blade_mesh)

    return blade_mesh

# Define the initial displacements of the control points
nx, ny = 5, 5
displacements = np.zeros((nx*ny, 2))
displacements[0] = 0.1, 0.5
displacements[10] = -0.2, 0.6
planform = generate_planform_curve(displacements, nx, ny)
blade_mesh = generate_blade_mesh(planform, 0.3)

# Save the STL file
pymesh.save_mesh("blade.stl", blade_mesh, ascii=True);