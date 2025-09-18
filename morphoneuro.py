bl_info = {
    "name": "MorphoNeuro",
    "author": "Filippo Andrea Sinosi",
    "version": (0, 1),
    "blender": (2, 80, 0),
    "location": "View3D > Sidebar > MorphoNeuro",
    "description": "All in one suite to import, visualize and analyze 3D meshes of neurosurgical corridors.",
    "warning": "",
    "wiki_url": "",
    "category": "Object",
}

import bpy
import bmesh
import math
import mathutils
import json
import numpy as np
import csv
import os
from mathutils import Vector, Matrix
import re
import datetime

# -----------------------------------------------------------
# Scene Properties for global options
# -----------------------------------------------------------

# Define the update function for projection_method in the global scope
def update_projection_method(self, context):
    # Update boolean flags when the projection method changes
    context.scene.use_spherical_geometry = (self.projection_method == 'SPHERICAL')
    context.scene.use_least_squares = (self.projection_method == 'LEAST_SQUARES')
    context.scene.use_centroid_direction = (self.projection_method == 'CENTROID')
    return None

# New: Value Memory item class for collection property
class CalculationValueItem(bpy.types.PropertyGroup):
    value: bpy.props.StringProperty(
        name="Value",
        description="The calculated numerical value"
    )
    unit: bpy.props.StringProperty(
        name="Unit",
        description="Unit of measurement (mm³, mm², °, etc.)"
    )
    label: bpy.props.StringProperty(
        name="Label",
        description="Description of what this value represents"
    )
    timestamp: bpy.props.StringProperty(
        name="Timestamp",
        description="When this calculation was performed"
    )

def init_props():
    # Volume of Surgical Corridor collapsible section
    bpy.types.Scene.show_vsf_options = bpy.props.BoolProperty(
        name="Volume of Surgical Corridor Options",
        description="Show or hide options for creating volumes of surgical corridor",
        default=False
    )
    
    # Volume of Surgical Corridor properties
    bpy.types.Scene.vsf_num_vertices = bpy.props.IntProperty(
        name="Num Vertices",
        description="Number of base vertices after decimation",
        default=1000,
        min=3,
        max=10000
    )
    
    bpy.types.Scene.vsf_distance_from_apex = bpy.props.IntProperty(
        name="Distance from Apex",
        description="Distance from apex to base center or radius of the sphere",
        default=100,
        min=1,
        max=1000
    )
    
    
    # Neuronavigation collapsible section
    bpy.types.Scene.show_neuronavigation_options = bpy.props.BoolProperty(
        name="Show Neuronavigation Options",
        default=False
    )
    bpy.types.Scene.nv_coords_filepath = bpy.props.StringProperty(
         name="File Path",
         description="Path to JSON, DAT or TXT file containing neuronavigation coordinates",
         default="",
         subtype='FILE_PATH'
    )
    bpy.types.Scene.nv_grouping_number = bpy.props.IntProperty(
        name="Grouping Number",
        description="Number of coordinates per group (0 for no grouping)",
        default=9,
        min=0
    )
    bpy.types.Scene.nv_create_without_apex = bpy.props.BoolProperty(
        name="Flat surfaces",
        description="Create base mesh only without connecting to an apex vertex",
        default=False
    )
    
    # Manual selection collapsible section
    bpy.types.Scene.show_manual_selection_options = bpy.props.BoolProperty(
        name="Show Manual Selection Options",
        default=False
    )
    bpy.types.Scene.selection_create_normalized_plane = bpy.props.BoolProperty(
        name="Create Normalized Plane",
        description="Create a normalized plane mesh in addition to the original mesh",
        default=True
    )
    bpy.types.Scene.selection_projection_distance = bpy.props.IntProperty(
        name="Projection Distance",
        description="Distance to project the base vertices from the apex",
        default=100,
        min=1
    )
    bpy.types.Scene.use_spherical_geometry = bpy.props.BoolProperty(
        name="Use Spherical Geometry",
        description="Use spherical geometry for the base",
        default=False
    )
    
    # ICP Registration collapsible section
    bpy.types.Scene.show_icp_options = bpy.props.BoolProperty(
        name="Show ICP Options",
        default=False
    )
    bpy.types.Scene.icp_max_iterations = bpy.props.IntProperty(
        name="Max Iterations", default=50, min=1, description="Maximum number of ICP iterations")
    bpy.types.Scene.icp_convergence_threshold = bpy.props.FloatProperty(
        name="Convergence Threshold", default=0.0001, min=0.00001, max=0.1, description="Stop iterations when error improvement is below this threshold")
    bpy.types.Scene.icp_use_pca_alignment = bpy.props.BoolProperty(
        name="Use PCA Pre-alignment", default=True, description="Use principal component analysis for initial alignment")
    bpy.types.Scene.icp_use_normals = bpy.props.BoolProperty(
        name="Use Normal Compatibility", default=False, description="Consider normal directions when finding corresponding points")
    bpy.types.Scene.icp_normal_weight = bpy.props.FloatProperty(
        name="Normal Weight", default=0.3, min=0.0, max=1.0, description="Weight given to normal compatibility (0 = ignore normals, 1 = only normals)")
    
    # Average Mesh collapsible section
    bpy.types.Scene.show_average_mesh_options = bpy.props.BoolProperty(
        name="Show Average Mesh Options",
        default=False
    )
    bpy.types.Scene.avg_use_icp_first = bpy.props.BoolProperty(
        name="Align meshes", default=True, description="Perform ICP registration to align meshes before averaging")
    
    # Normalize Mesh Height collapsible section
    bpy.types.Scene.show_normalize_height_options = bpy.props.BoolProperty(
        name="Show Normalize Height Options",
        default=False
    )
    bpy.types.Scene.projection_distance = bpy.props.IntProperty(
        name="Projection Distance", default=100, min=1, description="Distance to project the base vertices from the apex")
    bpy.types.Scene.projection_method = bpy.props.EnumProperty(
        name="Projection Method",
        description="Method used for projecting the base vertices",
        items=[
            ('SPHERICAL', "Spherical", "Use spherical geometry for plane fitting"),
            ('LEAST_SQUARES', "Least Squares", "Use least squares method for plane fitting"),
            ('CENTROID', "Centroid Direction", "Project vertices in the direction of the centroid")
        ],
        default='SPHERICAL',
        update=update_projection_method
    )
    
    # For backwards compatibility - these will be set based on projection_method
    bpy.types.Scene.use_spherical_geometry = bpy.props.BoolProperty(
        name="Use Spherical Geometry",
        description="Use spherical geometry for the base",
        default=False
    )
    bpy.types.Scene.use_least_squares = bpy.props.BoolProperty(
        name="Use Least Squares",
        description="Use least squares method for plane fitting",
        default=False
    )
    bpy.types.Scene.use_centroid_direction = bpy.props.BoolProperty(
        name="Use Centroid Direction",
        description="Use apex-centroid direction for plane fitting",
        default=False
    )
    
    # Vertex Distance collapsible section
    bpy.types.Scene.show_vertex_dist_options = bpy.props.BoolProperty(
        name="Show Vertex Distance Options",
        default=False
    )
    bpy.types.Scene.vertex_dist_use_icp = bpy.props.BoolProperty(
        name="Align meshes", default=True, description="Perform ICP registration to align meshes before calculating distances")
    bpy.types.Scene.vertex_dist_export_csv = bpy.props.BoolProperty(
        name="Export to CSV", default=False, description="Export distance measurements to a CSV file")
    bpy.types.Scene.vertex_dist_csv_filepath = bpy.props.StringProperty(
        name="CSV File Path", default="//vertex_distances.csv", subtype='FILE_PATH', description="Path to save the CSV file")
    
    # Pairwise Distance comparison options
    bpy.types.Scene.show_pairwise_options = bpy.props.BoolProperty(
        name="Show Pairwise Distance Options",
        default=False
    )
    bpy.types.Scene.pairwise_export_csv = bpy.props.BoolProperty(
        name="Export to CSV", default=False, description="Export pairwise comparisons to a CSV file")
    bpy.types.Scene.pairwise_csv_filepath = bpy.props.StringProperty(
        name="CSV File Path", default="//{blend_file_name}_{active_mesh_name}_pairwise_distances.csv", subtype='FILE_PATH', description="Path to save the CSV file")
    bpy.types.Scene.pairwise_use_icp = bpy.props.BoolProperty(
        name="Align meshes", default=True, description="Perform ICP registration to align meshes before comparing distances")
    
    # Calculated values memory system
    bpy.types.Scene.calculation_values = bpy.props.CollectionProperty(
        type=CalculationValueItem,
        name="Calculation Values",
        description="Stores multiple calculated values with their units and labels"
    )
    
    bpy.types.Scene.show_calculation_history = bpy.props.BoolProperty(
        name="Show Calculation History",
        default=True,
        description="Show or hide calculation history"
    )

def clear_props():
    # Volume of Surgical Corridor properties
    del bpy.types.Scene.show_vsf_options
    del bpy.types.Scene.vsf_num_vertices
    del bpy.types.Scene.vsf_distance_from_apex
    
    # Neuronavigation properties
    del bpy.types.Scene.show_neuronavigation_options
    del bpy.types.Scene.nv_coords_filepath
    del bpy.types.Scene.nv_grouping_number
    del bpy.types.Scene.nv_create_without_apex
    
    # Manual Selection properties
    del bpy.types.Scene.show_manual_selection_options
    del bpy.types.Scene.selection_create_normalized_plane
    del bpy.types.Scene.selection_projection_distance
    del bpy.types.Scene.use_spherical_geometry

    # ICP properties
    del bpy.types.Scene.show_icp_options
    del bpy.types.Scene.icp_max_iterations
    del bpy.types.Scene.icp_convergence_threshold
    del bpy.types.Scene.icp_use_pca_alignment
    del bpy.types.Scene.icp_use_normals
    del bpy.types.Scene.icp_normal_weight
    
    # Average Mesh properties
    del bpy.types.Scene.show_average_mesh_options
    del bpy.types.Scene.avg_use_icp_first
    
    # Normalize Height properties
    del bpy.types.Scene.show_normalize_height_options
    del bpy.types.Scene.projection_distance
    del bpy.types.Scene.projection_method
    del bpy.types.Scene.use_least_squares
    del bpy.types.Scene.use_centroid_direction

    # Vertex Distance properties
    del bpy.types.Scene.show_vertex_dist_options
    del bpy.types.Scene.vertex_dist_use_icp
    del bpy.types.Scene.vertex_dist_export_csv
    del bpy.types.Scene.vertex_dist_csv_filepath
    
    # Pairwise Distance properties
    del bpy.types.Scene.show_pairwise_options
    del bpy.types.Scene.pairwise_export_csv
    del bpy.types.Scene.pairwise_csv_filepath
    del bpy.types.Scene.pairwise_use_icp
    
    # Calculated values memory system
    del bpy.types.Scene.calculation_values
    del bpy.types.Scene.show_calculation_history

# -----------------------------------------------------------
# Utility functions (shared among operators)
# -----------------------------------------------------------

# New: Utility function to add a value to the calculation memory
def add_calculation_to_memory(context, value, unit, label, timestamp=None):
    """
    Add a new calculation to the memory array
    
    Args:
        context: Blender context
        value: The calculated value (as string)
        unit: Unit of measurement (mm³, mm², °, etc.)
        label: Description of what this value represents
        timestamp: When this calculation was performed (optional)
    
    Returns:
        Index of the newly added item
    """
    if timestamp is None:
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    item = context.scene.calculation_values.add()
    item.value = value
    item.unit = unit
    item.label = label
    item.timestamp = timestamp
    
    return len(context.scene.calculation_values) - 1

# New: Utility function to remove an item from calculation memory
def remove_calculation_from_memory(context, index):
    """
    Remove a calculation from the memory array by index
    
    Args:
        context: Blender context
        index: Index of the item to remove
    
    Returns:
        True if successful, False otherwise
    """
    if 0 <= index < len(context.scene.calculation_values):
        context.scene.calculation_values.remove(index)
        return True
    return False

def create_simple_line(parent, start, end, name="Ray"):
    """Create a simple line object to represent a ray from start to end."""
    mesh = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(name, mesh)

    # Adjust vertices to account for the object's origin
    origin = parent.location
    vertices = [(start - origin).copy(), (end - origin).copy()]

    edges = [(0, 1)]
    mesh.from_pydata(vertices, edges, [])
    bpy.context.collection.objects.link(obj)
    return obj

def calculate_centroid(vectors):
    return sum(vectors, Vector((0, 0, 0))) / len(vectors)

def find_apex_and_base(obj):
    """Find apex (vertex with most connections) and base vertices of a mesh.
    
    The apex is the vertex with the most connections. If no vertex has more than 2 connections,
    there is no apex and all vertices are treated as base vertices.
    Base vertices are only those directly connected to the apex via edges (if apex exists),
    otherwise all vertices are base vertices.
    
    Args:
        obj: Blender mesh object
        
    Returns:
        tuple: (apex_world_coord or None, base_world_coords)
    """
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    
    # Find the vertex with the most connections
    max_connections = max(len(v.link_edges) for v in bm.verts) if bm.verts else 0
    
    # If no vertex has more than 2 connections, there's no apex
    if max_connections <= 2:
        # All vertices are base vertices
        base_verts = list(bm.verts)
        base_verts.sort(key=lambda v: v.index)
        
        print(f"No apex found (max connections: {max_connections}). All {len(base_verts)} vertices treated as base.")
        
        # Convert to world coordinates
        apex_world = None
        base_world = [obj.matrix_world @ v.co for v in base_verts]
        
        bm.free()
        return apex_world, base_world
    
    # Find the vertex with the most connections (apex)
    apex_vert = max(bm.verts, key=lambda v: len(v.link_edges))
    
    # Get only the vertices connected to the apex via edges
    base_verts = []
    for edge in apex_vert.link_edges:
        other_vert = edge.other_vert(apex_vert)
        base_verts.append(other_vert)
    
    # Respect the order of the base vertices
    base_verts.sort(key=lambda v: v.index)
    
    print(f"Apex found with {len(apex_vert.link_edges)} connections. Number of base vertices: {len(base_verts)}")
    
    # Convert to world coordinates
    apex_world = obj.matrix_world @ apex_vert.co
    base_world = [obj.matrix_world @ v.co for v in base_verts]
    
    bm.free()
    return apex_world, base_world

def clear_recursive_children(obj, condition=None):
    for child in obj.children:
        if condition is None or condition in child.name:
            clear_recursive_children(child, condition)
            bpy.data.objects.remove(child, do_unlink=True)

def create_text_object(text, parent_obj, type):
    """Create a text object at the specified location."""

    text = str(text)

    # Delete any existing text similar to this one in the parent object
    for child in parent_obj.children:
        if child.type == 'FONT' and child.name.startswith(str(type)):
            bpy.data.objects.remove(child, do_unlink=True)

    # Get the current 3D view region and space
    area = next(area for area in bpy.context.screen.areas if area.type == 'VIEW_3D')
    region = area.regions[-1]
    rv3d = region.data
    
    # Get the view rotation
    view_rotation = rv3d.view_rotation
    
    # Calculate the centroid of the parent object
    if parent_obj.type == 'MESH' and len(parent_obj.data.vertices) > 0:
        # Calculate centroid from mesh vertices
        co_sum = Vector((0, 0, 0))
        for v in parent_obj.data.vertices:
            co_sum += v.co
        centroid = co_sum / len(parent_obj.data.vertices)
        # Transform to world space
        text_position = parent_obj.matrix_world @ centroid
    else:
        # Fallback to object origin if not a mesh or no vertices
        text_position = parent_obj.matrix_world @ parent_obj.location
    
    # Add a small offset in the direction of the view
    view_vector = view_rotation @ Vector((0.0, 0.0, 5.0))
    text_position += view_vector * 0.5

    # Set the unit of the text object
    unit = ""
    if type == "Volume":
        unit = "mm³"
    elif type == "Area":
        unit = "mm²"
    elif type == "Angle":
        unit = "°"

    # Create the text object at the calculated position
    font_curve = bpy.data.curves.new(type="FONT", name="Font Curve")
    font_curve.body = text + " " + unit
    font_obj = bpy.data.objects.new(f"{type}: {text}", font_curve)
    font_obj.location = text_position
    
    # Make the text bigger
    font_obj.scale = (5.0, 5.0, 5.0)
    
    # Make text face the camera
    font_obj.rotation_mode = 'QUATERNION'
    font_obj.rotation_quaternion = view_rotation

    # Make it stand out more
    font_obj.show_wire = True
    font_obj.show_in_front = True
    
    # Link to scene
    bpy.context.collection.objects.link(font_obj)
    
    # Set parent relationship
    font_obj.parent = parent_obj
    font_obj.matrix_parent_inverse = parent_obj.matrix_world.inverted()
    
    return font_obj

def make_normals_consistent(obj):
    """Ensure that mesh normals are oriented consistently."""
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode='EDIT')
    bpy.ops.mesh.select_all(action="SELECT")
    bpy.ops.mesh.normals_make_consistent(inside=False)
    bpy.ops.object.mode_set(mode='OBJECT')

def create_spherical_cap_mesh(apex_coord, base_coords, name, radius=100, subdivision_level=4):
    """Create a new mesh with base vertices projected onto a sphere centered at the apex.
       All base vertices will be at exactly 'radius' distance from the apex.
       The base is triangulated to create a more spherical-looking cap."""
    bm = bmesh.new()
    apex_vert = bm.verts.new(apex_coord)
    
    # Project base vertices onto a sphere centered at apex with radius
    projected_verts = []
    for co in base_coords:
        direction = (co - apex_coord).normalized()
        projected_point = apex_coord + direction * radius
        projected_verts.append(bm.verts.new(projected_point))
    
    bm.verts.index_update()
    
    # Create the base face
    try:
        base_face = bm.faces.new(projected_verts)
    except Exception as e:
        print("Spherical base face creation failed:", e)
        base_face = None
    
    # Create the side faces connecting apex to base
    n = len(projected_verts)
    for i in range(n):
        v1 = projected_verts[i]
        v2 = projected_verts[(i + 1) % n]
        try:
            bm.faces.new([apex_vert, v1, v2])
        except Exception as e:
            print("Spherical side face creation failed:", e)
    
    # Triangulate the base face to create a more spherical appearance
    if base_face:
        # First, triangulate the face using poke
        bmesh.ops.poke(bm, faces=[base_face])
        
        # Now use a different approach for subdivision:
        # 1. Triangulate all faces (should already be triangulated from poke)
        bmesh.ops.triangulate(bm, faces=bm.faces)
        
        # 2. For each subdivision level, use subdivide_edges on all edges of the base
        for _ in range(subdivision_level):
            # Get all edges that don't connect to the apex
            base_edges = [e for e in bm.edges if apex_vert not in e.verts]
            
            # Subdivide those edges
            bmesh.ops.subdivide_edges(bm, edges=base_edges, cuts=1, use_grid_fill=True)
            
            # Project all vertices (except apex) onto the sphere
            for v in bm.verts:
                if v != apex_vert:
                    direction = (v.co - apex_coord).normalized()
                    v.co = apex_coord + direction * radius
    
    mesh = bpy.data.meshes.new(name)
    bm.to_mesh(mesh)
    bm.free()
    
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    make_normals_consistent(obj)
    
    return obj

def create_normalized_plane_mesh(apex_coord, projected_base_coords, name):
    """Create a new mesh from an apex and projected base vertices.
       The base vertices maintain their original selection order."""
    bm = bmesh.new()
    apex_vert = bm.verts.new(apex_coord)
    projected_verts = [bm.verts.new(co) for co in projected_base_coords]
    bm.verts.index_update()
    
    try:
        bm.faces.new(projected_verts)
    except Exception as e:
        print("Projected base face creation failed:", e)
    
    n = len(projected_verts)
    for i in range(n):
        v1 = projected_verts[i]
        v2 = projected_verts[(i + 1) % n]
        try:
            bm.faces.new([apex_vert, v1, v2])
        except Exception as e:
            print("Projected side face creation failed:", e)
    
    mesh = bpy.data.meshes.new(name)
    bm.to_mesh(mesh)
    bm.free()
    
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    make_normals_consistent(obj)
    
    return obj

# -----------------------------------------------------------
# Functions for Volume of Surgical Corridor
# -----------------------------------------------------------
def calculate_base_plane_normal(base_vertices):
    if len(base_vertices) < 3:
        raise ValueError("At least three vertices required to define a plane.")
    v1, v2, v3 = base_vertices[0], base_vertices[1], base_vertices[2]
    edge1 = v2 - v1
    edge2 = v3 - v1
    return edge1.cross(edge2).normalized()

def find_opposite_vertices(apex, base_vertices, centroid):
    max_angle = 0
    vertex_pair = (None, None)
    for vertex in base_vertices:
        for candidate in base_vertices:
            if candidate == vertex:
                continue
            angle = (apex - vertex).angle(apex - candidate)
            if angle > max_angle:
                max_angle = angle
                vertex_pair = (vertex, candidate)
    print(f"Maximum angle of attack: {math.degrees(max_angle):.2f} degrees")
    return vertex_pair

def create_volume_of_surgical_corridor(apex, view_center, radius, num_vertices, target_obj, distance_from_apex):
    depsgraph = bpy.context.evaluated_depsgraph_get()
    evaluated_obj = target_obj.evaluated_get(depsgraph)
    new_base_center = view_center
    old_base_center = Vector((0, 0, 0))
    vertices_for_centroid_adjustment = 1000
    vertices_precision = 10000
    iteration = 0

    while True:
        iteration += 1
        axis = (new_base_center - apex).normalized()
        new_base_center = apex + axis * distance_from_apex
        is_close_enough = math.isclose(new_base_center.length, old_base_center.length, rel_tol=1e-7)
        # Always use the higher number of vertices for adjustment calculations
        current_vertices_to_generate = vertices_precision if is_close_enough else vertices_for_centroid_adjustment 
        
        perp = axis.orthogonal()
        perp2 = axis.cross(perp)
        base_vertices = []
        # Use the fixed high vertex count for the loop calculations
        angle_step = (2 * math.pi) / current_vertices_to_generate 
        for i in range(current_vertices_to_generate):
            angle = angle_step * i
            direction = (math.cos(angle) * perp + math.sin(angle) * perp2).normalized()
            vertex = new_base_center + direction * radius
            base_vertices.append(vertex)
        
        offset = axis * 0.001
        local_ray_origin = target_obj.matrix_world.inverted() @ (apex + offset)
        adjusted_vertices = []
        for vertex in base_vertices:
            factor = 1
            current_target = new_base_center
            target_closest_to_center = current_target
            target_closest_to_vertex = vertex
            while True:
                current_target = target_closest_to_center.lerp(target_closest_to_vertex, 0.5)
                direction_world = (current_target - (apex + offset)).normalized()
                local_ray_direction = (target_obj.matrix_world.to_quaternion().inverted() @ direction_world).normalized()
                result, loc, norm, idx = evaluated_obj.ray_cast(local_ray_origin, local_ray_direction)
                if not result:
                    target_closest_to_center = current_target
                else:
                    target_closest_to_vertex = current_target
                factor += 1
                if factor >= 30:
                    adjusted_vertices.append(current_target)
                    break
        old_base_center = new_base_center
        new_base_center = calculate_centroid(adjusted_vertices)
        if is_close_enough:
            break

    mesh = bpy.data.meshes.new("Volume of surgical corridor")
    cone_obj = bpy.data.objects.new("Volume of surgical corridor", mesh)
    bpy.context.collection.objects.link(cone_obj)

    vertices = adjusted_vertices.copy()
    faces = []
    base_face = list(range(len(vertices)))
    faces.append(base_face)
    apex_index = len(vertices)
    vertices.append(apex)
    for i in range(len(base_face)):
        faces.append([base_face[i], base_face[(i + 1) % len(base_face)], apex_index])
    
    mesh.from_pydata(vertices, [], faces)
    mesh.update()

    # Decimate the mesh if the target number of vertices is less than the generated number
    if num_vertices < vertices_precision and vertices_precision > 0:
        # Ensure the cone object is active and selected for modifier application
        bpy.context.view_layer.objects.active = cone_obj
        cone_obj.select_set(True)
        
        decimate_mod = cone_obj.modifiers.new(name="DecimateVOF", type='DECIMATE')
        decimate_mod.ratio = num_vertices / vertices_precision
        
        try:
            bpy.ops.object.modifier_apply(modifier=decimate_mod.name)
            print(f"Decimated Volume of Surgical Corridor mesh to approximately {num_vertices} base vertices.")
        except RuntimeError as e:
            print(f"Could not apply decimation modifier: {e}")
            # Remove the modifier if application failed
            cone_obj.modifiers.remove(decimate_mod)
        
        # Deselect the object after operation
        cone_obj.select_set(False)
    elif num_vertices >= vertices_precision:
        print(f"Target vertex count ({num_vertices}) is not less than generated count ({vertices_precision}). No decimation applied.")

    return apex, adjusted_vertices # Note: adjusted_vertices still holds the pre-decimation coordinates

# -----------------------------------------------------------
# Restore active vertex selection after volume creation
# -----------------------------------------------------------
def restore_active_vertex(active_obj, stored_world_coord, tol=1e-5):
    bm = bmesh.from_edit_mesh(active_obj.data)
    for v in bm.verts:
        v.select = False
    for v in bm.verts:
        world_co = active_obj.matrix_world @ v.co
        if (world_co - stored_world_coord).length < tol:
            v.select = True
            bm.select_history.clear()
            bm.select_history.add(v)
            break
    bmesh.update_edit_mesh(active_obj.data)

# -----------------------------------------------------------
# Functions for manual selection to corridors
# -----------------------------------------------------------
def create_cone_mesh(apex_coord, base_coords, name):
    """Create a cone mesh from an apex and base vertices."""
    bm = bmesh.new()
    apex_vert = bm.verts.new(apex_coord)
    base_verts = [bm.verts.new(co) for co in base_coords]
    bm.verts.index_update()
    
    try:
        bm.faces.new(base_verts)
    except Exception as e:
        print("Base face creation failed:", e)
    
    n = len(base_verts)
    for i in range(n):
        v1 = base_verts[i]
        v2 = base_verts[(i + 1) % n]
        try:
            bm.faces.new([apex_vert, v1, v2])
        except Exception as e:
            print("Side face creation failed:", e)
    
    mesh = bpy.data.meshes.new(name)
    bm.to_mesh(mesh)
    bm.free()
    
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)
    make_normals_consistent(obj)
    
    return obj

def restore_selection(active_obj, stored_world_coords, tol=1e-5):
    """Restore vertex selection in Edit Mode based on stored world coordinates."""
    bm = bmesh.from_edit_mesh(active_obj.data)
    for v in bm.verts:
        v.select = False
    for v in bm.verts:
        world_co = active_obj.matrix_world @ v.co
        for stored in stored_world_coords:
            if (world_co - stored).length < tol:
                v.select = True
                break
    bmesh.update_edit_mesh(active_obj.data)

# -----------------------------------------------------------
# Functions for ICP Registration
# -----------------------------------------------------------
def get_mesh_vertices(obj):
    """Extract vertices from a mesh object in world space."""
    mesh = obj.data
    vertices = [obj.matrix_world @ v.co for v in mesh.vertices]
    return vertices

def get_mesh_vertices_and_normals(obj):
    """Extract vertices and normals from a mesh object in world space."""
    mesh = obj.data
    
    # Get world matrix components for transformations
    mat_world = obj.matrix_world
    mat_normal = mat_world.inverted().transposed().to_3x3()
    
    vertices = []
    normals = []
    
    # Create a temporary bmesh to access normals
    bm = bmesh.new()
    bm.from_mesh(mesh)
    bm.faces.ensure_lookup_table()
    
    # Triangulate faces to ensure consistent normal calculation
    bmesh.ops.triangulate(bm, faces=bm.faces)
    
    # Get vertices and normals from faces
    for face in bm.faces:
        face_normal = mat_normal @ face.normal
        face_normal.normalize()
        
        for vert in face.verts:
            vert_co = mat_world @ vert.co
            vertices.append(vert_co)
            normals.append(face_normal.copy())
    
    bm.free()
    
    return vertices, normals

def closest_point_on_mesh(point, target_obj, target_normals=None, normal_weight=0.0):
    """Find the closest point on target_obj to the given point."""
    # Convert point to local space
    local_point = target_obj.matrix_world.inverted() @ point
    
    # Use closest point on mesh operator
    result, closest_point, normal, face_index = target_obj.closest_point_on_mesh(local_point)
    
    # Convert back to world space
    if result:
        world_point = target_obj.matrix_world @ closest_point
        world_normal = target_obj.matrix_world.to_3x3() @ normal
        return world_point, world_normal.normalized()
    return point, Vector((0, 0, 1))  # Return original point if no closest point found

def calculate_spherical_cap_area(apex, base_vertices, radius):
    """Calculate the area of a spherical cap."""

    # Calculate the area of the spherical cap using Girard's theorem
    # Implement Bevis-Cambareri's algorithm
    # Reference: https://link.springer.com/article/10.1007/BF00897843
            
    def vec_to_latlon(v):
        """
        Convert a 3D vector (assumed relative to the sphere center)
        into geographic coordinates (latitude and longitude) in degrees.
        """
        r = v.length
        if r == 0:
            raise ValueError("Zero-length vector cannot be converted to lat/lon.")
        lat = math.degrees(math.asin(v.z / r))
        lon = math.degrees(math.atan2(v.y, v.x))
        return lat, lon

    def trnsfrm_lon(plat, plon, qlat, qlon):
        """
        Transform the geographic coordinates of point Q (qlat, qlon)
        into a coordinate system in which point P (plat, plon) is the north pole.
        
        All angles are in degrees; the result is in radians.
        """
        pi = math.pi
        dtr = pi / 180.0  # degrees to radians conversion factor
        t = math.sin((qlon - plon) * dtr) * math.cos(qlat * dtr)
        b = math.sin(qlat * dtr) * math.cos(plat * dtr) - \
            math.cos(qlat * dtr) * math.sin(plat * dtr) * math.cos((qlon - plon) * dtr)
        return math.atan2(t, b)

    # Convert vertices into geographic coordinates relative to apex
    vlat = []
    vlon = []
    for v in base_vertices:
        rel = v - apex  # shift to sphere-centered coordinates
        lat, lon = vec_to_latlon(rel)
        vlat.append(lat)
        vlon.append(lon)
    
    nv = len(base_vertices)
    two_pi = 2 * math.pi
    total = 0.0
    
    # Loop over each vertex to compute the transformed longitude differences
    for i in range(nv):
        if i == 0:
            # First vertex: use second and last vertices as neighbors
            flat = vlat[1]
            flon = vlon[1]
            blat = vlat[-1]
            blon = vlon[-1]
        elif i < nv - 1:
            # Middle vertices: previous and next
            flat = vlat[i + 1]
            flon = vlon[i + 1]
            blat = vlat[i - 1]
            blon = vlon[i - 1]
        else:
            # Last vertex: use first and second-to-last vertices
            flat = vlat[0]
            flon = vlon[0]
            blat = vlat[-2]
            blon = vlon[-2]
        
        # Compute the transformed longitudes of the two neighbors
        fangl = trnsfrm_lon(vlat[i], vlon[i], flat, flon)
        bang = trnsfrm_lon(vlat[i], vlon[i], blat, blon)
        fvb = bang - fangl

        if fvb < 0.0:
            fvb += two_pi
        total += fvb
    
    # Compute the spherical excess
    spherical_excess = total - math.pi * (nv - 2)
    
    # The area is the excess times the square of the sphere's radius
    spherical_cap_area = spherical_excess * (radius ** 2)
    
    # Compute the area of a full sphere
    sphere_area = 4 * math.pi * (radius ** 2)

    # If the spherical cap area is greater than half the sphere area, we are calculating the area of the inverted spherical cap
    # For our use case, we never obtain more than half the sphere area, so we can use this simple check
    if spherical_cap_area > sphere_area / 2.0:
        spherical_cap_area = sphere_area - spherical_cap_area
    
    return spherical_cap_area
    

def compute_principal_axes(vertices):
    """Compute principal axes of a point cloud using PCA."""
    # Calculate centroid
    centroid = sum(vertices, Vector((0, 0, 0))) / len(vertices)
    
    # Center the points
    centered_points = [v - centroid for v in vertices]
    
    # Compute covariance matrix
    cov_matrix = Matrix.Identity(3)
    for p in centered_points:
        for i in range(3):
            for j in range(3):
                cov_matrix[i][j] += p[i] * p[j]
    
    # Normalize
    for i in range(3):
        for j in range(3):
            cov_matrix[i][j] /= len(vertices)
    
    # Compute eigenvalues and eigenvectors using numpy
    evals, evecs = np.linalg.eigh(np.array([[cov_matrix[i][j] for j in range(3)] for i in range(3)]))
    
    # Sort by eigenvalues in descending order
    idx = evals.argsort()[::-1]
    evals = evals[idx]
    evecs = evecs[:, idx]
    
    # Convert eigenvectors to Blender vectors
    axes = [Vector((evecs[0][i], evecs[1][i], evecs[2][i])) for i in range(3)]
    
    return centroid, axes, evals

def initial_alignment(source_vertices, target_vertices):
    """Perform initial alignment based on principal component analysis."""
    # Compute principal axes
    source_centroid, source_axes, source_evals = compute_principal_axes(source_vertices)
    target_centroid, target_axes, target_evals = compute_principal_axes(target_vertices)
    
    # Create rotation matrices from principal axes
    source_rot = Matrix((source_axes[0], source_axes[1], source_axes[2])).transposed()
    target_rot = Matrix((target_axes[0], target_axes[1], target_axes[2])).transposed()
    
    # Compute rotation from source to target
    rotation = target_rot @ source_rot.inverted()
    
    # Create transformation matrix
    transform = Matrix.Identity(4)
    for i in range(3):
        for j in range(3):
            transform[i][j] = rotation[i][j]
    
    # Add translation to align centroids
    translation = target_centroid - rotation @ source_centroid
    for i in range(3):
        transform[i][3] = translation[i]
    
    # Generate alternative alignments by flipping axes
    alternative_transforms = []
    
    # Try all possible axis flips (8 possibilities)
    for flip_x in [1, -1]:
        for flip_y in [1, -1]:
            for flip_z in [1, -1]:
                if flip_x == 1 and flip_y == 1 and flip_z == 1:
                    continue  # Skip the original orientation
                
                alt_rot = Matrix.Identity(3)
                alt_rot[0][0] = flip_x
                alt_rot[1][1] = flip_y
                alt_rot[2][2] = flip_z
                
                alt_rotation = target_rot @ alt_rot @ source_rot.inverted()
                alt_transform = Matrix.Identity(4)
                
                for i in range(3):
                    for j in range(3):
                        alt_transform[i][j] = alt_rotation[i][j]
                
                alt_translation = target_centroid - alt_rotation @ source_centroid
                for i in range(3):
                    alt_transform[i][3] = alt_translation[i]
                
                alternative_transforms.append(alt_transform)
    
    return transform, alternative_transforms

def icp_iteration(source_vertices, source_normals, target_obj, target_vertices, target_normals, 
                 normal_weight=0.3, use_normals=False):
    """Perform one iteration of ICP algorithm with optional normal weighting."""
    # Find corresponding points
    correspondences = []
    
    for i, source_point in enumerate(source_vertices):
        target_point, target_normal = closest_point_on_mesh(source_point, target_obj)
        
        # If using normals, check normal compatibility
        normal_compatibility = 1.0
        if use_normals and i < len(source_normals):
            source_normal = source_normals[i]
            # Dot product of normals (1 if aligned, -1 if opposite)
            normal_dot = source_normal.dot(target_normal)
            # Only use points with compatible normals (facing same direction)
            if normal_dot < 0:
                continue
            normal_compatibility = (normal_dot + 1) / 2  # Scale to [0,1]
        
        weight = 1.0
        if use_normals:
            weight = (1 - normal_weight) + normal_weight * normal_compatibility
        
        correspondences.append((source_point, target_point, weight))
    
    if len(correspondences) < 3:
        return None, 0, 0  # Not enough correspondences
    
    # Compute weighted centroids
    total_weight = sum(w for _, _, w in correspondences)
    source_centroid = Vector((0, 0, 0))
    target_centroid = Vector((0, 0, 0))
    
    for source, target, weight in correspondences:
        source_centroid += source * weight
        target_centroid += target * weight
    
    source_centroid /= total_weight
    target_centroid /= total_weight
    
    # Compute weighted covariance matrix
    covariance = Matrix.Identity(3)
    for source, target, weight in correspondences:
        source_centered = source - source_centroid
        target_centered = target - target_centroid
        for i in range(3):
            for j in range(3):
                covariance[i][j] += weight * source_centered[i] * target_centered[j]
    
    # Compute SVD
    u, s, v = np.linalg.svd(np.array([[covariance[i][j] for j in range(3)] for i in range(3)]))
    
    # Compute rotation matrix
    rotation = Matrix(v.transpose()) @ Matrix(u.transpose())
    
    # Ensure proper rotation (determinant should be 1)
    if np.linalg.det(rotation) < 0:
        v[-1] *= -1
        rotation = Matrix(v.transpose()) @ Matrix(u.transpose())
    
    # Compute translation
    translation = target_centroid - rotation @ source_centroid
    
    # Create transformation matrix
    transform = Matrix.Identity(4)
    for i in range(3):
        for j in range(3):
            transform[i][j] = rotation[i][j]
        transform[i][3] = translation[i]
    
    # Calculate error
    error = 0
    for source, target, weight in correspondences:
        transformed = transform @ source.to_4d()
        transformed = Vector((transformed[0], transformed[1], transformed[2]))
        error += weight * (transformed - target).length_squared
    
    error /= total_weight
    
    return transform, error, len(correspondences)

def evaluate_alignment(source_vertices, target_obj, transform):
    """Evaluate the quality of an alignment by measuring the average distance."""
    # Apply transformation to source vertices
    transformed_vertices = []
    for v in source_vertices:
        transformed = transform @ v.to_4d()
        transformed_vertices.append(Vector((transformed[0], transformed[1], transformed[2])))
    
    # Find corresponding points
    total_dist = 0
    valid_count = 0
    
    for source_point in transformed_vertices:
        target_point, _ = closest_point_on_mesh(source_point, target_obj)
        dist = (target_point - source_point).length
        
        total_dist += dist
        valid_count += 1
    
    if valid_count == 0:
        return float('inf')
    
    return total_dist / valid_count

# -----------------------------------------------------------
# Functions for Average Mesh Creation
# -----------------------------------------------------------
def find_corresponding_vertices(meshes):
    """Find corresponding vertices across multiple meshes."""
    # Use the first mesh as reference
    reference_mesh = meshes[0]
    reference_vertices = [reference_mesh.matrix_world @ v.co for v in reference_mesh.data.vertices]
    
    # For each vertex in the reference mesh, find closest points in other meshes
    correspondences = []
    
    for i, ref_vert in enumerate(reference_vertices):
        # Start with the reference vertex
        correspondence = [ref_vert]
        
        # Find closest vertices in other meshes
        for mesh in meshes[1:]:
            # Convert to local space
            local_point = mesh.matrix_world.inverted() @ ref_vert
            
            # Find closest vertex
            closest_vert = None
            min_dist = float('inf')
            
            for v in mesh.data.vertices:
                dist = (v.co - local_point).length
                if dist < min_dist:
                    min_dist = dist
                    closest_vert = v
            
            if closest_vert:
                # Convert back to world space
                world_point = mesh.matrix_world @ closest_vert.co
                correspondence.append(world_point)
            else:
                # If no vertex found, use the reference vertex
                correspondence.append(ref_vert)
        
        correspondences.append(correspondence)
    
    return correspondences

def create_average_mesh(correspondences, name="Average Mesh"):
    """Create a new mesh by averaging corresponding vertices."""
    # Average the coordinates
    avg_vertices = []
    for correspondence in correspondences:
        avg_vert = Vector((0, 0, 0))
        for vert in correspondence:
            avg_vert += vert
        avg_vert /= len(correspondence)
        avg_vertices.append(avg_vert)
    
    # Create a new mesh with the same topology as the reference mesh
    mesh = bpy.data.meshes.new(name)
    obj = bpy.data.objects.new(name, mesh)
    
    # Get the reference mesh to copy topology
    reference_mesh = bpy.context.selected_objects[0].data
    
    # Create vertices
    mesh.vertices.add(len(avg_vertices))
    for i, co in enumerate(avg_vertices):
        mesh.vertices[i].co = co
    
    # Copy edges
    edges = [(e.vertices[0], e.vertices[1]) for e in reference_mesh.edges]
    mesh.edges.add(len(edges))
    for i, edge in enumerate(edges):
        mesh.edges[i].vertices = edge
    
    # Copy faces using from_pydata instead of directly setting polygon attributes
    faces = [tuple(p.vertices) for p in reference_mesh.polygons]
    
    # Clear existing mesh data
    mesh.clear_geometry()
    
    # Recreate the mesh with from_pydata
    mesh.from_pydata(avg_vertices, edges, faces)
    
    # Update the mesh
    mesh.update()
    
    # Link the object to the scene
    bpy.context.collection.objects.link(obj)
    
    return obj

# -----------------------------------------------------------
# Operator: Volume of Surgical Corridor
# -----------------------------------------------------------
class OBJECT_OT_create_volume_of_surgical_corridor(bpy.types.Operator):
    bl_idname = "object.create_volume_of_surgical_corridor"
    bl_label = "Create Volume of Surgical Corridor"
    bl_description = "Create a volume of surgical corridor using a selected apex and the current 3D view center"
    bl_options = {'REGISTER', 'UNDO'}
    
    # Using a fixed radius from the original code
    radius = 300
    num_vertices: bpy.props.IntProperty(name="Num Vertices", default=1000)
    distance_from_apex: bpy.props.IntProperty(name="Distance from Apex", default=100)
    calculate_max_angle: bpy.props.BoolProperty(name="Calculate Max Angle", default=True)
    
    def execute(self, context):
        # --- Get 3D View Center ---
        view_center = None
        for area in context.screen.areas:
            if area.type == 'VIEW_3D':
                for region in area.regions:
                    if region.type == 'WINDOW':
                        r3d = region.data
                        view_matrix = r3d.view_matrix
                        # Calculate the center point.  This is the point in world space
                        # that projects to the center of the view.
                        view_center_local = Vector((0, 0, 0))
                        view_center = view_matrix.inverted() @ view_center_local
                        
                        break  # Assuming only one 'WINDOW' region per 'VIEW_3D' area
                if view_center:
                    break  # Stop searching areas if we found the view
        else:
            self.report({'ERROR'}, "No 3D view found.")
            return {'CANCELLED'}

        target_obj = context.object
        if not target_obj or target_obj.type != "MESH":
            self.report({'ERROR'}, "Select a valid mesh object.")
            return {'CANCELLED'}
        
        original_mode = target_obj.mode
        if target_obj.mode != 'EDIT':
            bpy.ops.object.mode_set(mode='EDIT')
        bm = bmesh.from_edit_mesh(target_obj.data)
        selected_verts = [v for v in bm.verts if v.select]
        if len(selected_verts) != 1:
            self.report({'ERROR'}, "Please select exactly one vertex as the apex.")
            return {'CANCELLED'}
        apex_local = selected_verts[0].co.copy()
        apex = target_obj.matrix_world @ apex_local
        stored_apex_world = apex.copy()
        bpy.ops.object.mode_set(mode='OBJECT')
        
        apex_result, base_vertices = create_volume_of_surgical_corridor(apex, view_center, self.radius, self.num_vertices, target_obj, self.distance_from_apex)
        
        if self.calculate_max_angle:
            centroid = calculate_centroid(base_vertices)
            v1, v2 = find_opposite_vertices(apex_result, base_vertices, centroid)
        
        # Get the created Volume of Surgical Corridor object
        vof_obj = None
        for obj in bpy.data.objects:
            if obj.name.startswith("Volume of surgical corridor"):
                vof_obj = obj
                break
        
        # If we found the object and normalization is enabled, perform spherical geometry normalization
        if vof_obj:
            # Get the current apex and base vertices from the decimated mesh
            current_apex, current_base_vertices = find_apex_and_base(vof_obj)
            
            # Only create spherical cap if we have an apex
            if current_apex is not None:
                # Create a new spherical cap mesh based on the decimated mesh
                spherical_obj = create_spherical_cap_mesh(
                    current_apex, 
                    current_base_vertices, 
                    f"Volume of corridor normalized", 
                    radius=self.distance_from_apex
                )
                
                # Set the new normalized mesh as active
                bpy.context.view_layer.objects.active = spherical_obj
                vof_obj.select_set(False)
                spherical_obj.select_set(True)

                # Delete the original Volume of Surgical Corridor object
                bpy.data.objects.remove(vof_obj, do_unlink=True)
            else:
                print(f"Warning: Cannot create spherical cap for '{vof_obj.name}' - no apex found")
        
        self.report({'INFO'}, "Volume of surgical corridor created.")
        bpy.ops.object.mode_set(mode=original_mode)
        if original_mode == 'EDIT':
            restore_active_vertex(target_obj, stored_apex_world)
        return {'FINISHED'}

# -----------------------------------------------------------
# Operator: Neuronavigation importer
# -----------------------------------------------------------
class OBJECT_OT_create_neuronavigation(bpy.types.Operator):
    bl_idname = "object.create_neuronavigation"
    bl_label = "Create neuronavigation surgical corridors"
    bl_description = "Create surgical corridors from coordinates in a JSON, DAT, or TXT file"
    bl_options = {'REGISTER', 'UNDO'}
    
    use_spherical_geometry: bpy.props.BoolProperty(
        name="Use Spherical Geometry",
        description="Use spherical geometry for the base",
        default=False
    )
    
    grouping_number: bpy.props.IntProperty(
        name="Grouping Number",
        description="Number of coordinates per group (0 for no grouping)",
        default=9,
        min=0
    )
    
    create_without_apex: bpy.props.BoolProperty(
        name="Create Without Apex",
        description="Create base mesh only without connecting to an apex vertex",
        default=False
    )
    
    def parse_dat_txt_file(self, filepath):
        """Parse a .dat or .txt file containing neuronavigation coordinates in the specified format."""
        coordinates = []
        points = []
        point_pattern = re.compile(r'point=([-\d\.\,]+)')
        
        try:
            # Read file line by line to find all coordinate points
            with open(filepath, 'r') as f:
                for line in f:
                    match = point_pattern.search(line)
                    if match:
                        point_str = match.group(1)
                        coords = [float(x) for x in point_str.split(',')]
                        points.append(coords)
            
            # Group points if grouping is enabled
            if self.grouping_number > 0 and points:
                for i in range(0, len(points), self.grouping_number):
                    group = points[i:i+self.grouping_number]
                    if group:  # Only add non-empty groups
                        coordinates.append(group)
            else:
                # No grouping, treat all points as one group
                if points:
                    coordinates = [points]
            
            self.report({'INFO'}, f"Parsed {len(points)} coordinates into {len(coordinates)} groups")
            return coordinates
        except Exception as e:
            self.report({'ERROR'}, f"Failed to parse DAT/TXT file: {e}")
    
    def save_json_file(self, coordinates, original_filepath):
        """Save the parsed coordinates to a JSON file with the same base name."""
        try:
            # Create JSON file path with same name but .json extension
            base_path = os.path.splitext(original_filepath)[0]
            json_path = base_path + ".json"
            
            with open(json_path, 'w') as f:
                json.dump(coordinates, f, indent=2)
                
            self.report({'INFO'}, f"Saved coordinates to {json_path}")
            return json_path
            
        except Exception as e:
            self.report({'WARNING'}, f"Failed to save JSON file: {e}")
            return None
    
    def execute(self, context):
        filepath = bpy.path.abspath(context.scene.nv_coords_filepath)
        if not filepath:
            self.report({'ERROR'}, "Please specify a file path for neuronavigation coordinates.")
            return {'CANCELLED'}
        
        # Get grouping number from scene property
        self.grouping_number = context.scene.nv_grouping_number
        self.create_without_apex = context.scene.nv_create_without_apex
        
        # Check file extension to determine how to process it
        file_ext = os.path.splitext(filepath)[1].lower()
        
        if file_ext in ['.dat', '.txt']:
            # Parse DAT/TXT file format
            coordinates = self.parse_dat_txt_file(filepath)
            if coordinates is None:
                return {'CANCELLED'}
                
            # Save parsed coordinates to JSON
            json_path = self.save_json_file(coordinates, filepath)
            if json_path:
                # Update filepath to use the new JSON file
                filepath = json_path
        
        # Load coordinates from JSON file (either original or newly created)
        try:
            with open(filepath, 'r') as f:
                coordinates = json.load(f)
        except Exception as e:
            self.report({'ERROR'}, f"Failed to load coordinates: {e}")
            return {'CANCELLED'}
        
        # Each group in the file should be a list of coordinates,
        # with the first coordinate as the apex and the rest as the base.
        for i, group in enumerate(coordinates):
            name = f"Mesh {i+1}"

            mesh = bpy.data.meshes.new(name)
            obj = bpy.data.objects.new(name, mesh)
            bpy.context.collection.objects.link(obj)
            bm = bmesh.new()
            bm_verts = [bm.verts.new(Vector(coord)) for coord in group]
            bm.verts.index_update()
            
            if self.create_without_apex:
                # Create base-only mesh without apex
                base_verts = bm_verts  # All coordinates are base vertices
                # Create the base face (if there are enough vertices)
                if len(base_verts) >= 3:
                    try:
                        bm.faces.new(base_verts)
                    except Exception as e:
                        print("Base face creation failed:", e)
                apex_coord = None
                base_coords = [v.co.copy() for v in base_verts]
            else:
                # Original behavior: first coordinate is apex, rest are base
                apex = bm_verts[0]
                base_verts = bm_verts[1:]
                try:
                    bm.faces.new(base_verts)
                except Exception as e:
                    print("Base face creation failed:", e)
                
                n = len(base_verts)
                for j in range(n):
                    v1 = base_verts[j]
                    v2 = base_verts[(j + 1) % n]
                    try:
                        bm.faces.new([apex, v1, v2])
                    except Exception as e:
                        print("Side face creation failed:", e)
                
                apex_coord = apex.co.copy()
                base_coords = [v.co.copy() for v in base_verts]
            
            bm.to_mesh(mesh)
            bm.free()
            make_normals_consistent(obj)
            
            # Create a spherical cap mesh if requested (only when apex exists)
            if self.use_spherical_geometry and not self.create_without_apex and apex_coord is not None:
                create_spherical_cap_mesh(apex_coord, base_coords, name + " normalized", radius=100.0)

        self.report({'INFO'}, "Neuronavigation surgical corridors created.")
        return {'FINISHED'}

# -----------------------------------------------------------
# Operator: Manual selection
# -----------------------------------------------------------
class OBJECT_OT_manual_selection(bpy.types.Operator):
    bl_idname = "object.manual_selection"
    bl_label = "Manual selection"
    bl_description = "Create a cone mesh from selected vertices in Edit Mode"
    bl_options = {'REGISTER', 'UNDO'}
    
    create_normalized_plane: bpy.props.BoolProperty(
        name="Create Normalized Plane",
        description="Create a normalized plane mesh in addition to the original mesh",
        default=True
    )
    
    projection_distance: bpy.props.IntProperty(
        name="Projection Distance",
        description="Distance from apex to project the base points",
        default=100
    )

    use_spherical_geometry: bpy.props.BoolProperty(
        name="Use Spherical Geometry",
        description="Use spherical geometry for the base",
        default=False
    )
    
    def invoke(self, context, event):
        # Initialize properties from scene properties
        self.create_normalized_plane = context.scene.selection_create_normalized_plane
        self.projection_distance = context.scene.selection_projection_distance
        self.use_spherical_geometry = context.scene.use_spherical_geometry

        return self.execute(context)
    
    def execute(self, context):
        active_obj = context.active_object
        if not active_obj or active_obj.type != 'MESH':
            self.report({'ERROR'}, "Please select a mesh object.")
            return {'CANCELLED'}
        if context.mode != 'EDIT_MESH':
            self.report({'ERROR'}, "Please be in Edit Mode and select vertices.")
            return {'CANCELLED'}
        
        bm = bmesh.from_edit_mesh(active_obj.data)
        selected_verts = [v for v in bm.verts if v.select]
        if len(selected_verts) < 2:
            self.report({'ERROR'}, "Please select at least 2 vertices (one apex and at least one base vertex).")
            return {'CANCELLED'}
        
        # Get the selection history
        selection_history = list(bm.select_history)
        if not selection_history:
            self.report({'ERROR'}, "Please use sequential selection (e.g., Shift+RMB) to select vertices in order.")
            return {'CANCELLED'}
            
        # The first vertex is the apex
        apex_elem = selection_history[0]

        # Retain the selection
        base_verts = [v for v in selection_history[1:] if isinstance(v, bmesh.types.BMVert)]
        
        stored_world_coords = [active_obj.matrix_world @ v.co.copy() for v in selection_history]
        
        apex_local = apex_elem.co.copy()
        base_local = [v.co.copy() for v in base_verts]
        apex_world = active_obj.matrix_world @ apex_local
        base_world = [active_obj.matrix_world @ co for co in base_local]
        
        bpy.ops.object.mode_set(mode='OBJECT')
        create_cone_mesh(apex_world, base_world, "Surgical corridor")

        if self.create_normalized_plane:
            if self.use_spherical_geometry:
                # --- Begin Projection ---
                # Create a spherical cap mesh with all points at exactly 100 units from apex
                create_spherical_cap_mesh(apex_world, base_world, f"{active_obj.name} normalized", radius=self.projection_distance)
                # --- End Projection ---
            else:
                # --- Begin Projection ---
                centroid = calculate_centroid(base_world)
                axis = (centroid - apex_world).normalized()
                new_base_center = apex_world + axis * self.projection_distance
                
                projected_base_coords = []
                for v in base_world:
                    direction = v - apex_world
                    denom = direction.dot(axis)
                    t = 1.0 if abs(denom) < 1e-6 else ((new_base_center - apex_world).dot(axis)) / denom
                    projected_point = apex_world + t * direction
                    projected_base_coords.append(projected_point)
                
                # Create new mesh with projected vertices
                create_normalized_plane_mesh(apex_world, projected_base_coords, f"{active_obj.name} normalized")
                # --- End Projection ---

        # Restore the selection
        bpy.context.view_layer.objects.active = active_obj
        bpy.ops.object.mode_set(mode='EDIT')
        restore_selection(active_obj, stored_world_coords)
        
        # Update scene properties with the values used
        context.scene.selection_create_normalized_plane = self.create_normalized_plane
        context.scene.selection_projection_distance = self.projection_distance
        context.scene.use_spherical_geometry = self.use_spherical_geometry
        
        self.report({'INFO'}, "surgical corridor mesh created from selection.")
        return {'FINISHED'}

# -----------------------------------------------------------
# Operator: Volume Calculator
# -----------------------------------------------------------
class OBJECT_OT_calculate_volume(bpy.types.Operator):
    bl_idname = "object.calculate_volume"
    bl_label = "Calculate Volume"
    bl_description = "Calculate the volume of selected mesh objects"
    bl_options = {'REGISTER', 'UNDO'}
    
    def calculate_mesh_volume(self, obj):
        """Calculate the volume of a mesh object."""
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.transform(obj.matrix_world)

        # Determine if the mesh is on a sphere; get the apex and base vertices
        use_spherical_geometry = False
        apex, base_vertices = find_apex_and_base(obj)

        # Only check for spherical geometry if we have an apex
        if apex is not None:
            # Calculate the radius of the sphere
            radius = (apex - base_vertices[0]).length

            # If all base vertices are at the same distance from the apex, the mesh is on a sphere
            if all(abs((v - apex).length - radius) < 1e-3 for v in base_vertices):
                use_spherical_geometry = True
        
        if use_spherical_geometry and apex is not None:

            # Calculate the area of the spherical cap
            spherical_cap_area = calculate_spherical_cap_area(apex, base_vertices, radius)
            
            # Calculate the volume of the spherical cap using Girard's theorem
            # Volume = (1/3) * R * A where A is the area of the spherical polygon
            volume = (1.0 / 3.0) * radius * spherical_cap_area
            
            # Set centroid to None as we're using apex as the reference point
            centroid = None

        else:
            # Get the centroid of the mesh
            centroid = Vector((0.0, 0.0, 0.0))
            for v in bm.verts:
                centroid += v.co
            centroid /= len(bm.verts)
            
            # Subtract the centroid from all vertices
            for v in bm.verts:
                v.co -= centroid
            
            # Triangulate the mesh
            bmesh.ops.triangulate(bm, faces=bm.faces)

            # Calculate the volume of the mesh
            volume = bm.calc_volume()
    
        bm.free()
        
        return volume, centroid
    
    def execute(self, context):
        selected_objects = [obj for obj in context.selected_objects if obj.type == "MESH"]
        
        if not selected_objects:
            self.report({'ERROR'}, "Select at least one mesh object.")
            return {'CANCELLED'}
        
        volumes_info = []
        
        for obj in selected_objects:
            volume, centroid = self.calculate_mesh_volume(obj)
            
            # Format volume with 2 decimal places
            volume_str = f"{volume:.2f}"
            
            # Create text object near the mesh and parent it to the mesh
            create_text_object(volume_str, obj, "Volume")
            
            # Store information for reporting
            volumes_info.append((obj.name, volume))
            
            # Print to console
            print(f"Volume of {obj.name}: {volume_str} cubic millimeters")
        
        # Create a summary 
        summary = ", ".join([f"{name}: {vol:.2f} mm³" for name, vol in volumes_info])
        self.report({'INFO'}, f"Volumes calculated: {summary}")
        
        # Store calculated values in memory
        if volumes_info:
            # Store each object's volume separately
                for name, volume in volumes_info:
                    value_str = f"{volume:.2f}"
                    label = f"Volume of {name}"
                    
                # Store in the memory system
                    add_calculation_to_memory(
                        context, 
                        value_str,
                        "mm³", 
                        label
                    )
        
        return {'FINISHED'}


# -----------------------------------------------------------
# Operator: Area Calculator
# -----------------------------------------------------------
class OBJECT_OT_calculate_area(bpy.types.Operator):
    bl_idname = "object.calculate_area"
    bl_label = "Calculate Area"
    bl_description = "Calculate the surface area of selected mesh objects using advanced projection methods"
    bl_options = {'REGISTER', 'UNDO'}

    def normalize_and_apply_shoelace(self, obj, coords):
        """Normalize the points and apply the Shoelace formula."""
        
        # Convert base vertices to numpy array for calculations
        points = np.array([list(v) for v in coords])
        
        # Calculate centroid of base vertices
        centroid = np.mean(points, axis=0)
        
        # Center the points
        centered_points = points - centroid
        
        # Calculate the best-fit plane using SVD
        # The right singular vectors of the centered points matrix give the principal directions
        U, s, Vh = np.linalg.svd(centered_points, full_matrices=False)
        
        # The normal vector is the last right singular vector
        normal = Vh[-1]
        normal = normal / np.linalg.norm(normal)  # Ensure unit length
        
        # Create the plane equation: ax + by + cz + d = 0
        # where (a, b, c) is the normal, and d = -dot(normal, centroid)
        a, b, c = normal
        d = -np.dot(normal, centroid)
        
        # Project base vertices onto the plane
        projected_points = []
        for point in points:
            # Calculate the distance from point to plane
            dist = (a * point[0] + b * point[1] + c * point[2] + d) / np.sqrt(a*a + b*b + c*c)
            
            # Project the point onto the plane
            # Convert point to numpy array for calculation
            point_np = np.array([point[0], point[1], point[2]])
            normal_vec = np.array([a, b, c])
            projected_point_np = point_np - dist * normal_vec
            # Convert result back to Vector for storage
            projected_point = Vector(projected_point_np)
            projected_points.append(projected_point)
        

        # Visualize the projected points as a closed polygon mesh
        # Create a new mesh and object
        mesh = bpy.data.meshes.new(f"{obj.name} normalized area")
        polygon_obj = bpy.data.objects.new(f"{obj.name} normalized area", mesh)
        
        # Link the object to the scene
        bpy.context.collection.objects.link(polygon_obj)
        
        # Create the mesh from projected points
        bm_vis = bmesh.new()
        
        # Add vertices
        for point in projected_points:
            bm_vis.verts.new(point)
        
        # Ensure vertices are indexed
        bm_vis.verts.ensure_lookup_table()
        
        # Create a face from all vertices
        if len(bm_vis.verts) >= 3:  # Need at least 3 vertices for a face
            bm_vis.faces.new(bm_vis.verts)
        
        # Update the mesh
        bm_vis.to_mesh(mesh)
        bm_vis.free()
        
        # Set the object's display type to wire for better visualization
        polygon_obj.display_type = 'WIRE'
        
        # Calculate area using the Shoelace formula
        # We need to project vertices to 2D first
        # We'll use the first two principal components as the 2D plane
        basis_u = Vh[0]
        basis_v = Vh[1]
        
        # Project points onto the 2D plane defined by the principal components
        points_2d = []
        for point in projected_points:
            # Center the point
            centered = point - Vector(centroid)
            # Project onto the 2D basis
            u = np.dot(centered, basis_u)
            v = np.dot(centered, basis_v)
            points_2d.append((u, v))
        
        # Apply Shoelace formula to calculate polygon area
        n = len(points_2d)
        area = 0.0
        for i in range(n):
            j = (i + 1) % n
            area += points_2d[i][0] * points_2d[j][1]
            area -= points_2d[j][0] * points_2d[i][1]
        area = abs(area) / 2.0
        
        return area
    
    def calculate_mesh_area(self, obj):
        """
        Calculate the area of the base of our corridors or of an user-defined polygon if in Edit Mode
        """
        
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        
        # Transform to world space
        for v in bm.verts:
            v.co = obj.matrix_world @ v.co
        
        # Make a copy of the faces before triangulation
        original_faces = []
        for face in bm.faces:
            # Store vertices of each face
            verts = [v.co.copy() for v in face.verts]
            original_faces.append(verts)
        
        # If object mode
        if obj.mode == 'OBJECT':

            # Determine if the base is a spherical cap
            is_sphere = False
            apex, base_vertices = find_apex_and_base(obj)

            # Check if we have an apex (for cone/pyramid shapes) or just a base (for flat polygons)
            if apex is not None:
                # Calculate the radius of the sphere
                radius = (apex - base_vertices[0]).length

                # If all base vertices are at the same distance from the apex, the mesh is on a sphere
                if all(abs((v - apex).length - radius) < 1e-3 for v in base_vertices):
                    is_sphere = True
                
                if is_sphere:
                    # Calculate the area of the spherical cap
                    spherical_cap_area = calculate_spherical_cap_area(apex, base_vertices, radius)
                    # Calculate the area of the mesh
                    total_area = spherical_cap_area
                else:
                    total_area = self.normalize_and_apply_shoelace(obj, base_vertices)
            else:
                # No apex found - treat as flat polygon
                total_area = self.normalize_and_apply_shoelace(obj, base_vertices)
            
            bm.free()
        elif obj.mode == 'EDIT':
            
            # Build a new bmesh 
            bm = bmesh.from_edit_mesh(obj.data)

            # Get the selection history
            selection_history = list(bm.select_history)
            if not selection_history:
                self.report({'ERROR'}, "Please use sequential selection (e.g., Shift+RMB) to select vertices in order.")
                return 0

            # Retain the selection
            coords = [v for v in selection_history if isinstance(v, bmesh.types.BMVert)]

            # Get their coordinates
            base_vertices = [obj.matrix_world @ v.co for v in coords]

            # Ensure they are more than 2
            if len(selection_history) < 3:
                self.report({'ERROR'}, "Please select at least 3 vertices.")
                return 0
            
            total_area = self.normalize_and_apply_shoelace(obj, base_vertices)
        
        return total_area
    
    def execute(self, context):
        selected_objects = [obj for obj in context.selected_objects if obj.type == "MESH"]
        
        if not selected_objects:
            self.report({'ERROR'}, "Select at least one mesh object.")
            return {'CANCELLED'}
        
        areas_info = []
        
        for obj in selected_objects:
            area = self.calculate_mesh_area(obj)
            
            # Convert to square millimeters and round to 2 decimal places
            area_mm2 = round(area, 2)
            
            # Create text object with area information
            create_text_object(area_mm2, obj, "Area")
            
            # Store information for reporting
            areas_info.append((obj.name, area_mm2))
            
            # Print to console
            print(f"Area of {obj.name}: {area_mm2} square millimeters")
        
        # Create a summary 
        summary = ", ".join([f"{name}: {area:.2f} mm²" for name, area in areas_info])
        self.report({'INFO'}, f"Areas calculated: {summary}")
        
        # Store calculated values in memory
        if areas_info:
            # Store each object's area separately
                for name, area in areas_info:
                    value_str = f"{area:.2f}"
                    label = f"Area of {name}"
                    
                # Store in the memory system
                    add_calculation_to_memory(
                        context, 
                        value_str,
                        "mm²", 
                        label
                    )
        
        return {'FINISHED'}


# -----------------------------------------------------------
# Operator: ICP Registration
# -----------------------------------------------------------
class OBJECT_OT_icp_registration(bpy.types.Operator):
    bl_idname = "object.icp_registration"
    bl_label = "ICP Registration"
    bl_description = "Align selected meshes to active mesh using Iterative Closest Point algorithm"
    bl_options = {'REGISTER', 'UNDO'}
    
    max_iterations: bpy.props.IntProperty(
        name="Max Iterations",
        description="Maximum number of ICP iterations",
        default=50,
        min=1
    )
    
    convergence_threshold: bpy.props.FloatProperty(
        name="Convergence Threshold",
        description="Stop iterations when error improvement is below this threshold",
        default=0.0001,
        min=0.00001,
        max=0.1
    )
    
    use_pca_alignment: bpy.props.BoolProperty(
        name="Use PCA Pre-alignment",
        description="Use principal component analysis for initial alignment",
        default=True
    )
    
    use_normals: bpy.props.BoolProperty(
        name="Use Normal Compatibility",
        description="Consider normal directions when finding corresponding points",
        default=True
    )
    
    normal_weight: bpy.props.FloatProperty(
        name="Normal Weight",
        description="Weight given to normal compatibility (0 = ignore normals, 1 = only normals)",
        default=0.3,
        min=0.0,
        max=1.0
    )
    
    def invoke(self, context, event):
        # Use scene properties as defaults
        self.max_iterations = context.scene.icp_max_iterations
        self.convergence_threshold = context.scene.icp_convergence_threshold
        self.use_pca_alignment = context.scene.icp_use_pca_alignment
        self.use_normals = context.scene.icp_use_normals
        self.normal_weight = context.scene.icp_normal_weight
        return self.execute(context)
    
    def execute(self, context):
        # Get target (active) and source (selected) objects
        target_obj = context.active_object
        selected = [obj for obj in context.selected_objects if obj != target_obj and obj.type == 'MESH']
        
        if not target_obj or target_obj.type != 'MESH':
            self.report({'ERROR'}, "Active object must be a mesh.")
            return {'CANCELLED'}
        
        if not selected:
            self.report({'ERROR'}, "Select at least one additional mesh to align to the active mesh.")
            return {'CANCELLED'}
        
        # Store original transforms for all source objects
        original_matrices = {obj: obj.matrix_world.copy() for obj in selected}
        
        # Determine if we should use vertex sampling
        use_sampling = False
        max_verts = 0
        
        # Extract target vertices and normals (with optional sampling)
        target_vertices = [target_obj.matrix_world @ v.co for v in target_obj.data.vertices]
        
        # Create a temporary bmesh to access normals
        bm = bmesh.new()
        bm.from_mesh(target_obj.data)
        bm.faces.ensure_lookup_table()
        
        # Triangulate faces to ensure consistent normal calculation
        bmesh.ops.triangulate(bm, faces=bm.faces)
        
        # Get normals from faces
        target_normals = []
        mat_normal = target_obj.matrix_world.inverted().transposed().to_3x3()
        
        for face in bm.faces:
            face_normal = mat_normal @ face.normal
            face_normal.normalize()
            
            for vert in face.verts:
                target_normals.append(face_normal.copy())
        
        bm.free()
        
        # Process each source object
        for i, source_obj in enumerate(selected):
            self.report({'INFO'}, f"Aligning mesh {i+1} of {len(selected)}...")
            
            # Extract source vertices and normals (with optional sampling)
            source_vertices = [source_obj.matrix_world @ v.co for v in source_obj.data.vertices]
            
            # Perform initial alignment if requested
            best_transform = None
            best_error = float('inf')
            
            if self.use_pca_alignment:
                # Get initial alignment based on PCA
                initial_transform, alternative_transforms = initial_alignment(source_vertices, target_vertices)
                
                # Try all possible initial alignments
                all_transforms = [initial_transform] + alternative_transforms
                
                for transform in all_transforms:
                    # Apply transformation to source vertices for evaluation
                    test_vertices = []
                    for v in source_vertices:
                        transformed = transform @ v.to_4d()
                        test_vertices.append(Vector((transformed[0], transformed[1], transformed[2])))
                    
                    # Evaluate this alignment
                    error = evaluate_alignment(test_vertices, target_obj, Matrix.Identity(4))
                    
                    if error < best_error:
                        best_error = error
                        best_transform = transform
                
                # Apply best initial transform to source object
                source_obj.matrix_world = best_transform @ source_obj.matrix_world
                
                # Update source vertices after initial alignment
                source_vertices = [source_obj.matrix_world @ v.co for v in source_obj.data.vertices]
                
                self.report({'INFO'}, f"Initial alignment complete for mesh {i+1} with error: {best_error:.6f}")
            
            # ICP iterations
            prev_error = float('inf')
            for iteration in range(self.max_iterations):
                # Perform ICP iteration
                transform, error, num_correspondences = icp_iteration(
                    source_vertices, target_normals, target_obj, 
                    target_vertices, target_normals, 
                    self.normal_weight, self.use_normals
                )
                
                if transform is None:
                    self.report({'WARNING'}, f"ICP failed for mesh {i+1}: Not enough corresponding points.")
                    if iteration == 0:
                        # Restore original transform only if we failed on first iteration
                        source_obj.matrix_world = original_matrices[source_obj]
                    break
                
                # Apply transformation to source vertices
                for j in range(len(source_vertices)):
                    transformed = transform @ source_vertices[j].to_4d()
                    source_vertices[j] = Vector((transformed[0], transformed[1], transformed[2]))
                
                
                # Apply transformation to source normals (rotation only)
                rotation = Matrix((transform[0][:3], transform[1][:3], transform[2][:3]))
                for j in range(len(target_normals)):
                    target_normals[j] = rotation @ target_normals[j]
                
                # Apply transformation to object
                source_obj.matrix_world = transform @ source_obj.matrix_world
                
                # Check convergence
                error_improvement = abs(prev_error - error)
                if error_improvement < self.convergence_threshold:
                    self.report({'INFO'}, f"ICP converged for mesh {i+1} after {iteration+1} iterations with error: {error:.6f} ({num_correspondences} points)")
                    break
                
                prev_error = error
            
            else:
                self.report({'INFO'}, f"ICP completed {self.max_iterations} iterations for mesh {i+1} with final error: {error:.6f} ({num_correspondences} points)")
        
        # Update scene properties with the values used
        context.scene.icp_max_iterations = self.max_iterations
        context.scene.icp_convergence_threshold = self.convergence_threshold
        context.scene.icp_use_pca_alignment = self.use_pca_alignment
        context.scene.icp_use_normals = self.use_normals
        context.scene.icp_normal_weight = self.normal_weight
        
        self.report({'INFO'}, f"ICP registration completed for {len(selected)} meshes.")
        return {'FINISHED'}

# -----------------------------------------------------------
# Operator: Create Average Mesh
# -----------------------------------------------------------
class OBJECT_OT_create_average_mesh(bpy.types.Operator):
    bl_idname = "object.create_average_mesh"
    bl_label = "Create Average Mesh"
    bl_description = "Create a new mesh by averaging the coordinates of corresponding vertices across selected meshes"
    bl_options = {'REGISTER', 'UNDO'}
    
    use_icp_first: bpy.props.BoolProperty(
        name="Align meshes",
        description="Perform ICP registration to align meshes before averaging",
        default=True
    )
    
    def invoke(self, context, event):
        # Use scene properties as defaults
        self.use_icp_first = context.scene.avg_use_icp_first
        return self.execute(context)
    
    def execute(self, context):
        selected_objects = [obj for obj in context.selected_objects if obj.type == 'MESH']
        
        if len(selected_objects) < 2:
            self.report({'ERROR'}, "Please select at least two mesh objects.")
            return {'CANCELLED'}
        
        # Store original transforms
        original_matrices = {obj: obj.matrix_world.copy() for obj in selected_objects}
        
        # Optionally align meshes using ICP before averaging
        if self.use_icp_first:
            # Use the first mesh as reference
            reference_obj = selected_objects[0]
            
            # Make the reference mesh active
            context.view_layer.objects.active = reference_obj
            
            # Use the ICP operator directly to align all meshes at once
            # The ICP operator now handles multiple meshes
            bpy.ops.object.icp_registration(
                max_iterations=20,
                convergence_threshold=0.0001,
                use_pca_alignment=True,
                use_normals=False
            )
        
        # Find corresponding vertices across all meshes
        self.report({'INFO'}, "Finding corresponding vertices...")
        correspondences = find_corresponding_vertices(selected_objects)
        
        # Create the average mesh
        self.report({'INFO'}, "Creating average mesh...")
        avg_mesh = create_average_mesh(correspondences, "Average Mesh")
        
        # Restore original transforms if we aligned the meshes
        if self.use_icp_first:
            for obj in selected_objects:
                obj.matrix_world = original_matrices[obj]
        
        # Update scene properties with the values used
        context.scene.avg_use_icp_first = self.use_icp_first
        
        self.report({'INFO'}, "Average mesh created successfully.")
        return {'FINISHED'}

# -----------------------------------------------------------
# Operator: Normalize Mesh Height
# -----------------------------------------------------------
class OBJECT_OT_normalize_mesh_height(bpy.types.Operator):
    bl_idname = "object.normalize_mesh_height"
    bl_label = "Normalize corridor"
    bl_description = "Project base vertices of selected meshes to a normalized plane at a specified distance from the apex"
    bl_options = {'REGISTER', 'UNDO'}
    
    projection_distance: bpy.props.IntProperty(
        name="Projection Distance",
        description="Distance to project the base vertices from the apex",
        default=100,
        min=1
    )

    projection_method: bpy.props.EnumProperty(
        name="Projection Method",
        description="Method used for projecting the base vertices",
        items=[
            ('SPHERICAL', "Spherical", "Use spherical geometry for plane fitting"),
            ('LEAST_SQUARES', "Least Squares", "Use least squares method for plane fitting"),
            ('CENTROID', "Centroid Direction", "Project vertices in the direction of the centroid")
        ],
        default='SPHERICAL'
    )

    def invoke(self, context, event):
        # Use scene property as default
        self.projection_distance = context.scene.projection_distance
        
        # Set projection method based on scene properties
        if context.scene.use_spherical_geometry:
            self.projection_method = 'SPHERICAL'
        elif context.scene.use_least_squares:
            self.projection_method = 'LEAST_SQUARES'
        elif context.scene.use_centroid_direction:
            self.projection_method = 'CENTROID'
            
        return self.execute(context)
    
    def normalize_single_mesh(self, obj):
        """Normalize a single mesh."""
        # Get apex and base vertices
        apex, base_vertices = find_apex_and_base(obj)

        # Check if we have an apex - normalization requires an apex
        if apex is None:
            print(f"Warning: Cannot normalize mesh '{obj.name}' - no apex found (flat polygon)")
            return None

        use_spherical_geometry = self.projection_method == 'SPHERICAL'
        use_least_squares = self.projection_method == 'LEAST_SQUARES'
        use_centroid_direction = self.projection_method == 'CENTROID'
        
        if not use_spherical_geometry:

            if not use_least_squares:
                centroid = calculate_centroid(base_vertices)
                
                # Calculate direction from apex to centroid
                axis = (centroid - apex).normalized()
                
                # Calculate new base center at specified distance
                new_base_center = apex + axis * self.projection_distance
                
                # Project base vertices onto the new plane
                projected_base_coords = []
                for v in base_vertices:
                    direction = v - apex
                    denom = direction.dot(axis)
                    t = 1.0 if abs(denom) < 1e-6 else ((new_base_center - apex).dot(axis)) / denom
                    projected_point = apex + t * direction
                    projected_base_coords.append(projected_point)
            else:
                # Convert base vertices to numpy array for calculations
                points = np.array([list(v) for v in base_vertices])
                
                # Calculate centroid of base vertices
                centroid = np.mean(points, axis=0)
                
                # Center the points
                centered_points = points - centroid
                
                # Calculate the best-fit plane using SVD
                # The right singular vectors of the centered points matrix give the principal directions
                U, s, Vh = np.linalg.svd(centered_points, full_matrices=False)
                
                # The normal vector is the last right singular vector
                normal = Vh[-1]
                normal = normal / np.linalg.norm(normal)  # Ensure unit length
                
                # Create the plane equation: ax + by + cz + d = 0
                # where (a, b, c) is the normal, and d = -dot(normal, centroid)
                a, b, c = normal

                # Calculate the perpendicular axis to the ax + by + cz + d = 0 plane
                perpendicular_axis = Vector((a, b, c)).normalized()
                
                # Ensure the direction is from apex toward the base vertices
                # Calculate average direction from apex to base vertices
                avg_direction = Vector((0, 0, 0))
                for v in base_vertices:
                    avg_direction += (v - apex).normalized()
                avg_direction /= len(base_vertices)
                
                # If the perpendicular axis is pointing in the opposite direction, flip it
                if avg_direction.dot(perpendicular_axis) < 0:
                    perpendicular_axis = -perpendicular_axis
                
                # Calculate new base center at specified distance
                new_base_center = apex + perpendicular_axis * self.projection_distance
                
                # Project base vertices onto the new plane
                projected_base_coords = []
                for v in base_vertices:
                    direction = v - apex
                    denom = direction.dot(perpendicular_axis)
                    t = 1.0 if abs(denom) < 1e-6 else ((new_base_center - apex).dot(perpendicular_axis)) / denom
                    projected_point = apex + t * direction
                    projected_base_coords.append(projected_point)
            
            # Create new mesh with projected vertices
            create_normalized_plane_mesh(apex, projected_base_coords, f"{obj.name} normalized")
        else:
            create_spherical_cap_mesh(apex, base_vertices, f"{obj.name} normalized", radius=self.projection_distance)
        
        return True
    
    def execute(self, context):
        # Get all selected mesh objects
        selected_meshes = [obj for obj in context.selected_objects if obj.type == 'MESH']
        
        if not selected_meshes:
            self.report({'ERROR'}, "Please select at least one mesh object.")
            return {'CANCELLED'}
        
        # Process each selected mesh
        normalized_count = 0
        for obj in selected_meshes:
            try:
                if self.normalize_single_mesh(obj):
                    normalized_count += 1
            except Exception as e:
                self.report({'WARNING'}, f"Failed to normalize {obj.name}: {str(e)}")
                # Add backtrace to the error message
                import traceback
                traceback.print_exc()
        
        # Update scene property with the value used
        context.scene.projection_distance = self.projection_distance
        
        # Update scene properties based on the selected projection method
        context.scene.use_spherical_geometry = (self.projection_method == 'SPHERICAL')
        context.scene.use_least_squares = (self.projection_method == 'LEAST_SQUARES')
        context.scene.use_centroid_direction = (self.projection_method == 'CENTROID')

        self.report({'INFO'}, f"Successfully normalized {normalized_count} of {len(selected_meshes)} meshes.")
        return {'FINISHED'}


# -----------------------------------------------------------
# Operator: Calculate Vertex Distances
# -----------------------------------------------------------
class OBJECT_OT_calculate_vertex_distances(bpy.types.Operator):
    bl_idname = "object.calculate_vertex_distances"
    bl_label = "Calculate Vertex Distances"
    bl_description = "Calculate Euclidean distances between corresponding vertices of two meshes"
    bl_options = {'REGISTER', 'UNDO'}
    
    use_icp_first: bpy.props.BoolProperty(
        name="Align meshes",
        description="Perform ICP registration to align meshes before calculating distances",
        default=True
    )
    
    export_csv: bpy.props.BoolProperty(
        name="Export to CSV",
        description="Export distance measurements to a CSV file",
        default=False
    )
    
    csv_filepath: bpy.props.StringProperty(
        name="CSV File Path",
        description="Path to save the CSV file",
        default="//vertex_distances.csv",
        subtype='FILE_PATH'
    )
    
    def invoke(self, context, event):
        # Use scene properties as defaults
        self.use_icp_first = context.scene.vertex_dist_use_icp
        self.export_csv = context.scene.vertex_dist_export_csv
        self.csv_filepath = context.scene.vertex_dist_csv_filepath
        return self.execute(context)
    
    def execute(self, context):
        # Get reference (active) and target (selected) objects
        reference_obj = context.active_object
        target_objects = [obj for obj in context.selected_objects if obj != reference_obj and obj.type == 'MESH']
        
        if not reference_obj or reference_obj.type != 'MESH':
            self.report({'ERROR'}, "Active object must be a mesh.")
            return {'CANCELLED'}
        
        if not target_objects:
            self.report({'ERROR'}, "Select at least one additional mesh to compare with the active mesh.")
            return {'CANCELLED'}
        
        # Store original transforms
        original_matrices = {obj: obj.matrix_world.copy() for obj in target_objects}
        original_matrices[reference_obj] = reference_obj.matrix_world.copy()
        
        # Get vertices from reference mesh in world space
        reference_vertices = [reference_obj.matrix_world @ v.co for v in reference_obj.data.vertices]
        
        # Prepare data structures for results
        all_distances = {}  # Dictionary to store distances for each target mesh
        vertex_indices = {}  # Dictionary to store corresponding vertex indices
        stats = {}  # Dictionary to store statistics for each target mesh
        
        # If exporting to CSV, prepare the file
        csv_path = None
        csv_writer = None
        if self.export_csv:
            csv_path = bpy.path.abspath(self.csv_filepath)
            try:
                csv_file = open(csv_path, 'w', newline='')
                csv_writer = csv.writer(csv_file)
                
                # Write header row with column titles
                header = ['Reference Vertex']
                for target_obj in target_objects:
                    header.append(f"{target_obj.name} Vertex")
                    header.append(f"{target_obj.name} Distance")
                csv_writer.writerow(header)
            except Exception as e:
                self.report({'ERROR'}, f"Failed to create CSV file: {e}")
                self.export_csv = False
        
        # Process each target mesh one by one
        for target_idx, target_obj in enumerate(target_objects):
            self.report({'INFO'}, f"Processing mesh {target_idx+1}/{len(target_objects)}: {target_obj.name}")
            
            # Optionally align mesh using ICP before calculating distances
            if self.use_icp_first:
                # Make sure the reference mesh is active
                context.view_layer.objects.active = reference_obj
                
                # Deselect all objects and select only the reference and current target
                for obj in context.selected_objects:
                    obj.select_set(False)
                reference_obj.select_set(True)
                target_obj.select_set(True)
                
                # Use the ICP operator to align the target mesh to the reference
                bpy.ops.object.icp_registration(
                    max_iterations=20,
                    convergence_threshold=0.0001,
                    use_pca_alignment=True,
                    use_normals=False,
                    normal_weight=0.3
                )
            
            # Calculate distances between corresponding vertices
            distances = []
            vertex_pairs = []
            max_dist = 0
            min_dist = float('inf')
            total_dist = 0
            
            # For each vertex in the reference mesh, find the closest vertex in the target mesh
            for i, ref_vert in enumerate(reference_vertices):
                # Convert to target's local space
                local_point = target_obj.matrix_world.inverted() @ ref_vert
                
                # Find closest vertex
                closest_vert = None
                min_vertex_dist = float('inf')
                closest_idx = -1
                
                for j, v in enumerate(target_obj.data.vertices):
                    dist = (v.co - local_point).length
                    if dist < min_vertex_dist:
                        min_vertex_dist = dist
                        closest_vert = v
                        closest_idx = j
                
                if closest_vert:
                    # Convert back to world space
                    world_point = target_obj.matrix_world @ closest_vert.co
                    # Calculate Euclidean distance
                    distance = (world_point - ref_vert).length
                    distances.append(distance)
                    vertex_pairs.append((i, closest_idx))
                    
                    # Update statistics
                    max_dist = max(max_dist, distance)
                    min_dist = min(min_dist, distance)
                    total_dist += distance
                else:
                    self.report({'WARNING'}, f"No corresponding vertex found for vertex {i} in {target_obj.name}")
                    distances.append(float('nan'))
                    vertex_pairs.append((i, -1))
            
            # Calculate average distance
            avg_dist = total_dist / len(distances) if distances else 0
            
            # Store results for this target mesh
            all_distances[target_obj.name] = distances
            vertex_indices[target_obj.name] = vertex_pairs
            stats[target_obj.name] = {
                'min': min_dist,
                'max': max_dist,
                'avg': avg_dist,
                'count': len(distances)
            }
            
            # Restore original transform if we aligned the mesh
            if self.use_icp_first:
                target_obj.matrix_world = original_matrices[target_obj]
                
        # If exporting to CSV, write the data rows
        if self.export_csv and csv_writer:
            try:
                # Write data rows
                for i in range(len(reference_vertices)):
                    row = [i]  # Reference vertex index
                    for target_obj in target_objects:
                        if i < len(vertex_indices[target_obj.name]):
                            target_idx = vertex_indices[target_obj.name][i][1]
                            distance = all_distances[target_obj.name][i]
                            row.append(target_idx)
                            row.append(f"{distance:.6f}" if not math.isnan(distance) else "N/A")
                        else:
                            row.append("N/A")
                            row.append("N/A")
                    csv_writer.writerow(row)
                
                # Write summary statistics
                csv_writer.writerow([])
                csv_writer.writerow(['Statistics'])
                
                # Min distances
                row = ['Minimum Distance']
                for target_obj in target_objects:
                    row.append('')  # Empty cell for vertex index
                    row.append(f"{stats[target_obj.name]['min']:.6f}")
                csv_writer.writerow(row)
                
                # Max distances
                row = ['Maximum Distance']
                for target_obj in target_objects:
                    row.append('')  # Empty cell for vertex index
                    row.append(f"{stats[target_obj.name]['max']:.6f}")
                csv_writer.writerow(row)
                
                # Avg distances
                row = ['Average Distance']
                for target_obj in target_objects:
                    row.append('')  # Empty cell for vertex index
                    row.append(f"{stats[target_obj.name]['avg']:.6f}")
                csv_writer.writerow(row)
                
                # Vertex counts
                row = ['Total Vertices']
                for target_obj in target_objects:
                    row.append('')  # Empty cell for vertex index
                    row.append(str(stats[target_obj.name]['count']))
                csv_writer.writerow(row)
                
                # Close the CSV file
                csv_file.close()
                
                self.report({'INFO'}, f"Distances exported to {csv_path}")
            except Exception as e:
                self.report({'ERROR'}, f"Failed to write to CSV: {e}")
        
        # Select all objects again
        for obj in target_objects:
            obj.select_set(True)
        reference_obj.select_set(True)
        
        # Update scene properties with the values used
        context.scene.vertex_dist_use_icp = self.use_icp_first
        context.scene.vertex_dist_export_csv = self.export_csv
        context.scene.vertex_dist_csv_filepath = self.csv_filepath
        
        self.report({'INFO'}, f"Vertex distances calculated for {len(target_objects)} meshes.")
        return {'FINISHED'}

# -----------------------------------------------------------
# Operator: Compare Pairwise Distances
# -----------------------------------------------------------
class OBJECT_OT_compare_pairwise_distances(bpy.types.Operator):
    bl_idname = "object.compare_pairwise_distances"
    bl_label = "Compare Pairwise Distances"
    bl_description = "Compare pairwise distances between vertices of two meshes"
    bl_options = {'REGISTER', 'UNDO'}
    
    export_csv: bpy.props.BoolProperty(
        name="Export to CSV",
        description="Export pairwise distances to a CSV file",
        default=False
    )
    
    csv_filepath: bpy.props.StringProperty(
        name="CSV File Path",
        description="Path to save the CSV file",
        default="//{blend_file_name}_{active_mesh_name}_pairwise_distances.csv",
        subtype='FILE_PATH'
    )
    
    use_icp_first: bpy.props.BoolProperty(
        name="Align meshes",
        description="Perform ICP registration to align meshes before comparing distances",
        default=True
    )
    
    def invoke(self, context, event):
        # Use scene properties as defaults
        self.csv_filepath = context.scene.pairwise_csv_filepath
        self.export_csv = context.scene.pairwise_export_csv
        self.use_icp_first = context.scene.pairwise_use_icp
        return self.execute(context)
    
    def execute(self, context):
        # Get selected mesh objects
        selected_meshes = [obj for obj in context.selected_objects if obj.type == 'MESH']
        
        if len(selected_meshes) != 2:
            self.report({'ERROR'}, "Please select exactly two mesh objects.")
            return {'CANCELLED'}
        
        mesh_a = selected_meshes[0]
        mesh_b = selected_meshes[1]
        
        self.report({'INFO'}, f"Comparing pairwise distances between {mesh_a.name} and {mesh_b.name}")
        
        # Store original transforms
        original_matrix_a = mesh_a.matrix_world.copy()
        original_matrix_b = mesh_b.matrix_world.copy()
        
        # Optionally align meshes using ICP before calculating vertex correspondences
        if self.use_icp_first:
            # Make mesh_a active
            context.view_layer.objects.active = mesh_a
            
            # Deselect all objects and select only the two meshes
            for obj in context.selected_objects:
                obj.select_set(False)
            mesh_a.select_set(True)
            mesh_b.select_set(True)
            
            # Use the ICP operator to align mesh_b to mesh_a
            bpy.ops.object.icp_registration(
                max_iterations=20,
                convergence_threshold=0.0001,
                use_pca_alignment=True,
                use_normals=False,
                normal_weight=0.3
            )
        
        # Get vertices from both meshes in world space
        verts_a = [mesh_a.matrix_world @ v.co for v in mesh_a.data.vertices]
        verts_b = [mesh_b.matrix_world @ v.co for v in mesh_b.data.vertices]
        
        sampled_indices_a = list(range(len(verts_a)))
        
        # Find vertex correspondences: for each vertex in mesh A, find closest vertex in mesh B
        vertex_correspondences = {}  # Maps vertex index in A to closest vertex index in B
        
        for i, vert_a in enumerate(verts_a):
            # Find closest vertex in mesh B
            closest_idx_b = -1
            min_dist = float('inf')
            
            for j, vert_b in enumerate(verts_b):
                dist = (vert_b - vert_a).length
                if dist < min_dist:
                    min_dist = dist
                    closest_idx_b = j
            
            if closest_idx_b >= 0:
                # Map the original vertex index from mesh A to the closest vertex in B
                original_idx_a = sampled_indices_a[i]
                vertex_correspondences[original_idx_a] = closest_idx_b
        
        self.report({'INFO'}, f"Established {len(vertex_correspondences)} vertex correspondences")
        
        # Generate all unique vertex pairs from the correspondences
        # For each pair of vertices in A, calculate the distance and compare with 
        # the corresponding pair in B
        pairwise_comparisons = []
        
        # Get the indices of vertices in A that have a correspondence in B
        valid_indices_a = list(vertex_correspondences.keys())
        
        # Generate all unique pairs of these indices
        pair_count = 0
        
        for i in range(len(valid_indices_a)):
            for j in range(i + 1, len(valid_indices_a)):
                
                idx_a1 = valid_indices_a[i]
                idx_a2 = valid_indices_a[j]
                
                # Get corresponding vertices in B
                idx_b1 = vertex_correspondences[idx_a1]
                idx_b2 = vertex_correspondences[idx_a2]
                
                # Calculate distance between vertices in A
                pos_a1 = mesh_a.matrix_world @ mesh_a.data.vertices[idx_a1].co
                pos_a2 = mesh_a.matrix_world @ mesh_a.data.vertices[idx_a2].co
                dist_a = (pos_a2 - pos_a1).length
                
                # Calculate distance between corresponding vertices in B
                pos_b1 = mesh_b.matrix_world @ mesh_b.data.vertices[idx_b1].co
                pos_b2 = mesh_b.matrix_world @ mesh_b.data.vertices[idx_b2].co
                dist_b = (pos_b2 - pos_b1).length
                
                # Store the comparison
                pairwise_comparisons.append({
                    'pair_a': (idx_a1, idx_a2),
                    'pair_b': (idx_b1, idx_b2),
                    'dist_a': dist_a,
                    'dist_b': dist_b,
                    'diff': abs(dist_a - dist_b)
                })
                
                pair_count += 1
        
        # Sort pairwise comparisons by difference (closest matches first)
        pairwise_comparisons.sort(key=lambda x: x['diff'])
        
        # Restore original transforms if we aligned the meshes
        if self.use_icp_first:
            mesh_a.matrix_world = original_matrix_a
            mesh_b.matrix_world = original_matrix_b
        
        # Create a text object to display the summary
        if pairwise_comparisons:
            # Calculate some statistics
            diffs = [data['diff'] for data in pairwise_comparisons]
            min_diff = min(diffs)
            max_diff = max(diffs)
            avg_diff = sum(diffs) / len(diffs)
        
        # Write to CSV
        if self.export_csv:
            csv_path = bpy.path.abspath(self.csv_filepath)
            try:
                with open(csv_path, 'w', newline='') as csvfile:
                    writer = csv.writer(csvfile)
                    
                    # Write header row with column titles as requested in the format
                    writer.writerow(['', mesh_a.name, mesh_b.name])
                    
                    # Write data rows
                    for data in pairwise_comparisons:
                        pair_a_str = f"Pair {data['pair_a'][0]}-{data['pair_a'][1]}"
                        writer.writerow([
                            pair_a_str,
                            f"{data['dist_a']:.6f}",
                            f"{data['dist_b']:.6f}"
                        ])
                    
                    # Add summary statistics
                    writer.writerow([])
                    writer.writerow(['Statistics', '', ''])
                    writer.writerow(['Total Pairs', len(pairwise_comparisons), ''])
                    
                    if pairwise_comparisons:
                        writer.writerow(['Min Difference', f"{min_diff:.6f}", ''])
                        writer.writerow(['Max Difference', f"{max_diff:.6f}", ''])
                        writer.writerow(['Avg Difference', f"{avg_diff:.6f}", ''])
                        
                        # Add vertex correspondence information
                        writer.writerow([])
                        writer.writerow(['Vertex Correspondences', '', ''])
                        writer.writerow([f'Vertex {mesh_a.name}', f'Vertex {mesh_b.name}', ''])
                        
                        # Sort correspondences by vertex index in A
                        for idx_a, idx_b in sorted(vertex_correspondences.items()):
                            writer.writerow([idx_a, idx_b, ''])
                
                self.report({'INFO'}, f"Pairwise distances exported to {csv_path}")
                
            except Exception as e:
                self.report({'ERROR'}, f"Failed to write to CSV: {e}")
                return {'CANCELLED'}
        
        # Update scene properties with the values used
        context.scene.pairwise_export_csv = self.export_csv
        context.scene.pairwise_csv_filepath = self.csv_filepath
        context.scene.pairwise_use_icp = self.use_icp_first
        
        return {'FINISHED'}

# -----------------------------------------------------------
# Operator: Calculate Angles
# -----------------------------------------------------------
class OBJECT_OT_calculate_angles(bpy.types.Operator):
    bl_idname = "object.calculate_angles"
    bl_label = "Calculate Angles"
    bl_description = "Calculate angles between apex and selected vertices"
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        import math
        import bmesh
        
        # Different behavior based on mode
        if context.mode == 'EDIT_MESH':
            # Get the active object
            obj = context.active_object
            if not obj or obj.type != 'MESH':
                self.report({'ERROR'}, "Active object must be a mesh")
                return {'CANCELLED'}
            
            # Create a bmesh from the edit mesh
            bm = bmesh.from_edit_mesh(obj.data)
            
            # Get selected vertices
            selected_verts = [v for v in bm.verts if v.select]
            
            if len(selected_verts) < 2 or len(selected_verts) > 3:
                self.report({'ERROR'}, "Select exactly two vertices in addition to the apex")
                return {'CANCELLED'}
            
            # If only two vertices, we infer the apex from the base
            if len(selected_verts) == 2:
                # Find apex and base
                apex, base = find_apex_and_base(obj)
            
                if apex is None:
                    self.report({'ERROR'}, "Could not determine apex vertex")
                    return {'CANCELLED'}
            
                # apex is already in world coordinates from find_apex_and_base
                apex_global = apex

                v1 = selected_verts[0]
                v2 = selected_verts[1]
            # If three vertices, the second is the apex
            elif len(selected_verts) == 3:
                apex = selected_verts[1]
                apex_global = obj.matrix_world @ apex.co

                v1 = selected_verts[0]
                v2 = selected_verts[2]

            v1_global = obj.matrix_world @ v1.co
            v2_global = obj.matrix_world @ v2.co
            
            vec1 = (v1_global - apex_global).normalized()
            vec2 = (v2_global - apex_global).normalized()
            
            # Calculate angle in radians
            dot_product = vec1.dot(vec2)
            # Clamp dot product to [-1, 1] to handle floating point errors
            dot_product = max(min(dot_product, 1.0), -1.0)
            angle_rad = math.acos(dot_product)
            angle = math.degrees(angle_rad)
            
            # Log the angle
            self.report({'INFO'}, f"Calculated angle: {angle:.2f} degrees between vertices {v1.index} and {v2.index}")
            
            # Store the angle in the memory system
            value_str = f"{angle:.2f}"
            label = f"Angle between vertices {v1.index} and {v2.index}"
            add_calculation_to_memory(context, value_str, "°", label)
            

        else:  # OBJECT mode
            # Get selected mesh objects
            selected_objects = [obj for obj in context.selected_objects if obj.type == 'MESH']
            
            if not selected_objects:
                self.report({'ERROR'}, "No mesh objects selected")
                return {'CANCELLED'}
            
            all_results = []
            
            for obj in selected_objects:
                # Find apex and base for this object
                apex, base = find_apex_and_base(obj)
                
                if apex is None or not base:
                    self.report({'WARNING'}, f"Could not determine apex and base vertices for {obj.name}")
                    continue
                
                # apex and base are already in world coordinates
                apex_global = apex
                
                # Number of base vertices
                n = len(base)
                
                # Check if there are enough vertices to create pairs (at least 4)
                if n < 4:
                    self.report({'WARNING'}, f"At least 4 base vertices are needed to calculate perpendicular angles for {obj.name}")
                    continue
                
                # Calculate primary and perpendicular angles using the fixed pattern
                angle_pairs = []
                
                # For primary pairs: vertex i is paired with vertex i+n/2
                for i in range(n):
                    opposite_idx = (i + n//2) % n
                    
                    # Calculate angle between primary vertex pair
                    v1 = base[i]
                    v2 = base[opposite_idx]
                    vec1 = (v1 - apex_global).normalized()
                    vec2 = (v2 - apex_global).normalized()
                    
                    # Calculate angle in radians
                    dot_product = vec1.dot(vec2)
                    # Clamp dot product to [-1, 1] to handle floating point errors
                    dot_product = max(min(dot_product, 1.0), -1.0)
                    angle_rad = math.acos(dot_product)
                    primary_angle_deg = math.degrees(angle_rad)
                    
                    # For perpendicular pairs: (i+n/4, i+3n/4) to the primary pair (i, i+n/2)
                    perp1_idx = (i + n//4) % n
                    perp2_idx = (i + 3*n//4) % n
                    
                    # Calculate angle between perpendicular vertex pair
                    v3 = base[perp1_idx]
                    v4 = base[perp2_idx]
                    vec3 = (v3 - apex_global).normalized()
                    vec4 = (v4 - apex_global).normalized()
                    
                    # Calculate angle in radians
                    dot_product = vec3.dot(vec4)
                    # Clamp dot product to [-1, 1] to handle floating point errors
                    dot_product = max(min(dot_product, 1.0), -1.0)
                    angle_rad = math.acos(dot_product)
                    perp_angle_deg = math.degrees(angle_rad)
                    
                    # Store both angles along with indices
                    combined_angle = primary_angle_deg + perp_angle_deg

                    primary = (i, opposite_idx, primary_angle_deg)
                    perpendicular = (perp1_idx, perp2_idx, perp_angle_deg)

                    # The primary angle is always the bigger one
                    if primary_angle_deg < perp_angle_deg:
                        # Invert the primary and perpendicular angles
                        perpendicular = (i, opposite_idx, primary_angle_deg)
                        primary = (perp1_idx, perp2_idx, perp_angle_deg)

                    angle_pairs.append({
                        'primary': primary,
                        'perpendicular': perpendicular,
                        'combined': combined_angle,
                        'object': obj
                    })
                
                # Find the pair with maximum combined angle for this object
                if angle_pairs:
                    max_pair = max(angle_pairs, key=lambda x: x['combined'])
                    all_results.append(max_pair)
            
            # Process and report results
            if all_results:

                # Report all objects
                for result in all_results:
                    obj = result['object']
                    primary = result['primary']
                    perpendicular = result['perpendicular']

                    # Remove any existing angle visualization for the max object
                    clear_recursive_children(obj, "angle")

                    # Get apex and base for the max object to create visualizations
                    apex, base = find_apex_and_base(obj)
                    
                    # Get coordinates for primary vertices
                    p_v1_idx, p_v2_idx, primary_angle = primary
                    v1 = base[p_v1_idx if p_v1_idx < len(base) else 0]
                    v2 = base[p_v2_idx if p_v2_idx < len(base) else 1]
                    
                    # Get coordinates for perpendicular vertices
                    perp_v1_idx, perp_v2_idx, perp_angle = perpendicular
                    v3 = base[perp_v1_idx if perp_v1_idx < len(base) else 2]
                    v4 = base[perp_v2_idx if perp_v2_idx < len(base) else 3]
                    
                    # Create ray visualizations
                    ray1 = create_simple_line(obj, v1, v2, "Primary angle")
                    ray2 = create_simple_line(obj, v3, v4, "Perpendicular angle")

                    # Create text object
                    create_text_object(f"{primary[2]:.2f}", ray1, "Angle")
                    create_text_object(f"{perpendicular[2]:.2f}", ray2, "Angle")

                    # Store in memory system
                    primary_angle_str = f"{primary_angle:.2f}"
                    label = f"Primary Angle of {obj.name}"
                    add_calculation_to_memory(context, primary_angle_str, "°", label)
                    
                    perp_angle_str = f"{perp_angle:.2f}"
                    perp_label = f"Perpendicular Angle of {obj.name}"
                    add_calculation_to_memory(context, perp_angle_str, "°", perp_label)

                    # Set ray properties
                    for ray in [ray1, ray2]:
                        
                        # Parent ray to the original mesh
                        ray.parent = obj
                        
                        # Make the ray stand out more
                        ray.show_wire = True
                        ray.show_in_front = True
                    
                    self.report({'INFO'}, f"Object {obj.name}: Primary angle: {primary[2]:.2f}°, Perpendicular angle: {perpendicular[2]:.2f}°, Combined: {result['combined']:.2f}°")
                    
                self.report({'INFO'}, "Created ray visualizations for the angle pairs")
            else:
                self.report({'ERROR'}, "No angles to calculate across selected objects")

        return {'FINISHED'}

# -----------------------------------------------------------
# Operator: Copy Value to Clipboard
# -----------------------------------------------------------
class OBJECT_OT_copy_value_to_clipboard(bpy.types.Operator):
    bl_idname = "object.copy_value_to_clipboard"
    bl_label = "Copy to Clipboard"
    bl_description = "Copy the calculated value to clipboard"
    bl_options = {'REGISTER'}
    
    collection_index: bpy.props.IntProperty(
        name="Collection Index",
        default=-1,
        description="Index of the value in the collection to copy"
    )
    
    def execute(self, context):
        if self.collection_index >= 0 and self.collection_index < len(context.scene.calculation_values):
                value = context.scene.calculation_values[self.collection_index].value
                context.window_manager.clipboard = value
                self.report({'INFO'}, f"Value '{value}' copied to clipboard")
        else:
            self.report({'WARNING'}, "Invalid collection index")
        return {'FINISHED'}

# -----------------------------------------------------------
# UI Panel
# -----------------------------------------------------------
class VIEW3D_PT_morphoneuro(bpy.types.Panel):
    bl_label = "MorphoNeuro"
    bl_idname = "VIEW3D_PT_morphoneuro"
    bl_space_type = 'VIEW_3D'
    bl_region_type = 'UI'
    bl_category = 'MorphoNeuro'
    
    def draw(self, context):
        layout = self.layout
        scene = context.scene
        
        # 1 - CREATION SECTION
        layout.label(text="Modeling", icon='MATERIAL')
        layout.separator()
        
        # Neuronavigation Section (collapsible)
        box = layout.box()
        row = box.row()
        row.prop(scene, "show_neuronavigation_options", icon="TRIA_DOWN" if scene.show_neuronavigation_options else "TRIA_RIGHT", 
                 icon_only=True, emboss=False)
        # Take advantage of the full row width for the label
        label_split = row.split(factor=0.6)
        label_split.label(text="Neuronavigation")
        # Button in the remaining space
        button_split = label_split.split(factor=1)
        op = button_split.operator("object.create_neuronavigation", text="Create")
        op.grouping_number = scene.nv_grouping_number
        op.create_without_apex = scene.nv_create_without_apex

        # Show Neuronavigation options if expanded
        if scene.show_neuronavigation_options:
            col = box.column(align=True)
            col.prop(scene, "nv_coords_filepath")
            col.prop(scene, "nv_grouping_number")
            col.prop(scene, "nv_create_without_apex")
        
        # Manual selection Section (collapsible)
        box = layout.box()
        row = box.row()
        row.prop(scene, "show_manual_selection_options", icon="TRIA_DOWN" if scene.show_manual_selection_options else "TRIA_RIGHT", 
                 icon_only=True, emboss=False)
        # Take advantage of the full row width for the label
        label_split = row.split(factor=0.6)
        label_split.label(text="Manual selection")
        # Button in the remaining space
        button_split = label_split.split(factor=1)
        op = button_split.operator("object.manual_selection", text="Create")
        op.create_normalized_plane = scene.selection_create_normalized_plane
        op.use_spherical_geometry = scene.use_spherical_geometry
        
        # Show Manual selection options if expanded
        if scene.show_manual_selection_options:
            col = box.column(align=True)
            col.prop(scene, "selection_create_normalized_plane")

            if scene.selection_create_normalized_plane:
                col.prop(scene, "selection_projection_distance")
                col.prop(scene, "use_spherical_geometry")

        # Volume of Surgical Corridor Section (collapsible)
        box = layout.box()
        row = box.row()
        row.prop(scene, "show_vsf_options", icon="TRIA_DOWN" if scene.show_vsf_options else "TRIA_RIGHT", 
                 icon_only=True, emboss=False)
        # Take advantage of the full row width for the label
        label_split = row.split(factor=0.6)
        label_split.label(text="Automatic corridor")
        # Button in the remaining space
        button_split = label_split.split(factor=1)
        op = button_split.operator("object.create_volume_of_surgical_corridor", text="Create")
        op.num_vertices = scene.vsf_num_vertices
        op.distance_from_apex = scene.vsf_distance_from_apex
        
        # Show VSF options if expanded
        if scene.show_vsf_options:
            col = box.column(align=True)
            col.prop(scene, "vsf_num_vertices")
            col.prop(scene, "vsf_distance_from_apex")
        
        # 2 - POLISHING SECTION
        layout.separator()
        layout.label(text="Polishing", icon='BRUSH_DATA')
        layout.separator()
        
        # Normalize Mesh Height Section (collapsible)
        box = layout.box()
        row = box.row()
        row.prop(scene, "show_normalize_height_options", icon="TRIA_DOWN" if scene.show_normalize_height_options else "TRIA_RIGHT", 
                 icon_only=True, emboss=False)
        # Take advantage of the full row width for the label
        label_split = row.split(factor=0.6)
        label_split.label(text="Normalize meshes")
        # Button in the remaining space
        button_split = label_split.split(factor=1)
        op = button_split.operator("object.normalize_mesh_height", text="Run")
        op.projection_distance = scene.projection_distance
        op.projection_method = scene.projection_method
        
        # Show Normalize Height options if expanded
        if scene.show_normalize_height_options:
            col = box.column(align=True)
            col.prop(scene, "projection_distance")
            col.prop(scene, "projection_method")
            
        # Average Mesh Section (collapsible)
        box = layout.box()
        row = box.row()
        row.prop(scene, "show_average_mesh_options", icon="TRIA_DOWN" if scene.show_average_mesh_options else "TRIA_RIGHT", 
                 icon_only=True, emboss=False)
        # Take advantage of the full row width for the label
        label_split = row.split(factor=0.6)
        label_split.label(text="Average coordinates")
        # Button in the remaining space
        button_split = label_split.split(factor=1)
        op = button_split.operator("object.create_average_mesh", text="Run")
        op.use_icp_first = scene.avg_use_icp_first
        
        # Show Average Mesh options if expanded
        if scene.show_average_mesh_options:
            col = box.column(align=True)
            col.prop(scene, "avg_use_icp_first")
        
        # 3 - ANALYSIS SECTION
        layout.separator()
        layout.label(text="Analysis", icon='DRIVER_DISTANCE')
        layout.separator()
        
        # Calculate Volume, Area, Angle
        box = layout.box()
        row = box.row()
        row.operator("object.calculate_volume", text="Volume", icon='CONE')
        row.operator("object.calculate_area", text="Area", icon='MESH_CIRCLE')
        row.operator("object.calculate_angles", text="Angle", icon='LINCURVE')
        
        # Display the calculation history in a collapsible section
        box = layout.box()
        row = box.row()
        row.prop(scene, "show_calculation_history", icon="TRIA_DOWN" if scene.show_calculation_history else "TRIA_RIGHT", 
                 icon_only=True, emboss=False)
        
        # Split to put label and clear button on same row
        split = row.split(factor=0.6)
        split.label(text="Calculation History")
        
        # Only show clear button if there are values to clear
        if len(scene.calculation_values) > 0:
            split.operator("object.remove_all_calculations", text="Clear All", icon='TRASH')
        
        # Show calculation history if expanded
        if scene.show_calculation_history:
            # Display all stored calculation values
            if scene.calculation_values:
                # Add a small header
                column = box.column()
                
                # Loop through all stored values with most recent first
                for i, item in enumerate(reversed(scene.calculation_values)):
                    idx = len(scene.calculation_values) - 1 - i
                    
                    row = column.row(align=True)
                    split = row.split(factor=0.5)
                    
                    # Label column
                    label_col = split.column()
                    label_col.label(text=f"{item.label}:")
                    
                    # Value column with unit
                    value_col = split.row(align=True)
                    value_col.label(text=f"{item.value} {item.unit}")
                    
                    # Copy button
                    copy_op = value_col.operator("object.copy_value_to_clipboard", text="", icon='COPYDOWN')
                    copy_op.collection_index = idx
                    
                    # Remove button
                    value_col.operator("object.remove_calculation", text="", icon='X').index = idx
            else:
            # If no calculations are stored, show a message
                box.label(text="No calculations stored yet.")
        
        # Vertex Distance Section (collapsible)
        box = layout.box()
        row = box.row()
        row.prop(scene, "show_vertex_dist_options", icon="TRIA_DOWN" if scene.show_vertex_dist_options else "TRIA_RIGHT", 
                 icon_only=True, emboss=False)
        # Take advantage of the full row width for the label
        label_split = row.split(factor=0.6)
        label_split.label(text="Vertex Distances")
        # Button in the remaining space
        button_split = label_split.split(factor=1)
        op = button_split.operator("object.calculate_vertex_distances", text="Calculate")
        op.use_icp_first = scene.vertex_dist_use_icp
        op.export_csv = scene.vertex_dist_export_csv
        op.csv_filepath = scene.vertex_dist_csv_filepath
        
        # Show Vertex Distance options if expanded
        if scene.show_vertex_dist_options:
            col = box.column(align=True)
            col.prop(scene, "vertex_dist_use_icp")
            col.prop(scene, "vertex_dist_export_csv")
            if scene.vertex_dist_export_csv:
                col.prop(scene, "vertex_dist_csv_filepath")
        
        # Pairwise Distance Comparison Section (collapsible)
        box = layout.box()
        row = box.row()
        row.prop(scene, "show_pairwise_options", icon="TRIA_DOWN" if scene.show_pairwise_options else "TRIA_RIGHT", 
                 icon_only=True, emboss=False)
        # Take advantage of the full row width for the label
        label_split = row.split(factor=0.6)
        label_split.label(text="Pairwise Distances")
        # Button in the remaining space
        button_split = label_split.split(factor=1)
        op = button_split.operator("object.compare_pairwise_distances", text="Compare")
        op.use_icp_first = scene.pairwise_use_icp
        
        # Show Pairwise Distance options if expanded
        if scene.show_pairwise_options:
            col = box.column(align=True)
            col.prop(scene, "pairwise_use_icp")
            col.prop(scene, "pairwise_export_csv")
            if scene.pairwise_export_csv:
                col.prop(scene, "pairwise_csv_filepath")
        
        # 4 - OTHER TOOLS SECTION
        layout.separator()
        layout.label(text="Other Tools", icon='TOOL_SETTINGS')
        layout.separator()
        
        # ICP Registration Section (collapsible)
        box = layout.box()
        row = box.row()
        row.prop(scene, "show_icp_options", icon="TRIA_DOWN" if scene.show_icp_options else "TRIA_RIGHT", 
                 icon_only=True, emboss=False)
        # Take advantage of the full row width for the label
        label_split = row.split(factor=0.6)
        label_split.label(text="ICP Registration")
        # Button in the remaining space
        button_split = label_split.split(factor=1)
        op = button_split.operator("object.icp_registration", text="Run")
        op.max_iterations = scene.icp_max_iterations
        op.convergence_threshold = scene.icp_convergence_threshold
        op.use_pca_alignment = scene.icp_use_pca_alignment
        op.use_normals = scene.icp_use_normals
        op.normal_weight = scene.icp_normal_weight
        
        # Show ICP options if expanded
        if scene.show_icp_options:
            col = box.column(align=True)
            col.prop(scene, "icp_max_iterations")
            col.prop(scene, "icp_convergence_threshold")
            col.prop(scene, "icp_use_pca_alignment")
            col.prop(scene, "icp_use_normals")
            if scene.icp_use_normals:
                col.prop(scene, "icp_normal_weight")

# -----------------------------------------------------------
# Operator: Remove Calculation from Memory
# -----------------------------------------------------------
class OBJECT_OT_remove_calculation(bpy.types.Operator):
    bl_idname = "object.remove_calculation"
    bl_label = "Remove Calculation"
    bl_description = "Remove a calculation from memory"
    bl_options = {'REGISTER', 'UNDO'}
    
    index: bpy.props.IntProperty(
        name="Index",
        description="Index of the calculation to remove",
        default=0
    )
    
    def execute(self, context):
        if remove_calculation_from_memory(context, self.index):
            self.report({'INFO'}, "Calculation removed from memory")
        else:
            self.report({'ERROR'}, "Invalid calculation index")
        return {'FINISHED'}

# -----------------------------------------------------------
# Operator: Remove All Calculations from Memory
# -----------------------------------------------------------
class OBJECT_OT_remove_all_calculations(bpy.types.Operator):
    bl_idname = "object.remove_all_calculations"
    bl_label = "Clear History"
    bl_description = "Remove all calculations from memory"
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        if len(context.scene.calculation_values) > 0:
            # Remove all items by repeatedly removing the first item
            while len(context.scene.calculation_values) > 0:
                context.scene.calculation_values.remove(0)
            self.report({'INFO'}, "All calculations removed from memory")
        else:
            self.report({'INFO'}, "No calculations to remove")
        return {'FINISHED'}

# -----------------------------------------------------------
# Registration
# -----------------------------------------------------------
classes = [
    CalculationValueItem,
    OBJECT_OT_create_neuronavigation,
    OBJECT_OT_manual_selection,
    OBJECT_OT_create_volume_of_surgical_corridor,
    OBJECT_OT_normalize_mesh_height,
    OBJECT_OT_create_average_mesh,
    OBJECT_OT_calculate_volume,
    OBJECT_OT_calculate_area,
    OBJECT_OT_calculate_angles,
    OBJECT_OT_calculate_vertex_distances,
    OBJECT_OT_compare_pairwise_distances,
    OBJECT_OT_icp_registration,
    OBJECT_OT_copy_value_to_clipboard,
    OBJECT_OT_remove_calculation,
    OBJECT_OT_remove_all_calculations,
    VIEW3D_PT_morphoneuro,
]

def register():
    for cls in classes:
        bpy.utils.register_class(cls)
    init_props()

def unregister():
    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)
    clear_props()

if __name__ == "__main__":
    register()
