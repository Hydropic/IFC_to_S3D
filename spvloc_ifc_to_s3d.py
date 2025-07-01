import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.shape_builder
import json
import numpy as np
import sys
from scipy.spatial import ConvexHull

# File paths
ifc_file_path = "/raid/USERDATA/ovrobi4y/SPVLOC/spvloc/2025-03-10_Stand_AA/amogus_together_moved.ifc"
output_json_path = "/raid/USERDATA/ovrobi4y/SPVLOC/spvloc/2025-03-10_Stand_AA/zind/scene_00100/annotation_3d.json"
junctions = []
junctions_dict = {}
door_window_coordinates = []
lines = []
planes = []
semantics = []
cuboids = []
total_vertices = 0
manhattan_structures = []
ifc_model = ifcopenshell.open(ifc_file_path)
specific_ifc = True

# HELPER FUNCTIONS

def get_global_transform(placement):
    """Recursively compute global transform for a given IfcLocalPlacement"""
    transform = np.identity(4)

    while placement:
        if placement.is_a("IfcLocalPlacement"):
            rel = placement.RelativePlacement
            loc = rel.Location.Coordinates if rel else [0, 0, 0]
            loc = np.array(loc, dtype=float)
            rotation = np.identity(3)

            if hasattr(rel, "RefDirection") and hasattr(rel, "Axis"):
                # Build rotation matrix from RefDirection and Axis
                ref_dir = np.array(rel.RefDirection.DirectionRatios)
                axis = np.array(rel.Axis.DirectionRatios)
                z_axis = axis / np.linalg.norm(axis)
                x_axis = ref_dir / np.linalg.norm(ref_dir)
                y_axis = np.cross(z_axis, x_axis)
                rotation = np.vstack((x_axis, y_axis, z_axis)).T

            # Combine rotation and translation into a 4x4 matrix
            local_tf = np.identity(4)
            local_tf[:3, :3] = rotation
            local_tf[:3, 3] = loc
            transform = np.dot(local_tf, transform)

            placement = getattr(placement, "PlacementRelTo", None)
        else:
            break

    return transform

    
def _compute_normal(p1, p2, p3):
    v1 = np.array(p2) - np.array(p1)
    v2 = np.array(p3) - np.array(p1)
    normal = np.cross(v1, v2)
    norm_length = np.linalg.norm(normal)
    EPSILON = 1e-8
    if norm_length > EPSILON:
        return (normal / norm_length).tolist()
    else:
        return [0.0, 0.0, 0.0]


def _process_mesh(verts, edges, element_id, obj_type, obj_name):
    scale = 1000 # meter to millimeter
    points = [(verts[i]*scale, verts[i+1]*scale, verts[i+2]*scale) for i in range(0, len(verts), 3)]
    global junctions
    global junctions_dict
    global door_window_coordinates
    if obj_type == "IfcDoor" or obj_type == "IfcWindow":
        # get the coordinates of the door/window with the corresponding obj_type (used in extract_semantics assign obj_type to opening)
        door_window_coordinates.append({
            "coordinates": points,
            "obj_type": obj_type
        })
        return
    for i in range(0, len(edges), 2):
        start = points[edges[i]]
        end = points[edges[i + 1]]

        if start not in junctions_dict:
            junctions_dict[start] = len(junctions)
            junctions.append({
                "ID": len(junctions),
                "coordinate": list(start),
                "element_id": element_id,
                "obj_type": obj_type,
                "obj_name": obj_name
            })
        
        direction = tuple(e - s for s, e in zip(start, end))
        #threshold = 0.001 
        #direction = tuple(0.0 if abs(x) < threshold else x for x in direction)
        #if np.linalg.norm(direction) <= 1:  # Avoid zero-length lines
        #    print(f"⚠ Warning: Zero-length line detected between {start} and {end}. Skipping.")
        #    continue

        if end not in junctions_dict:
            junctions_dict[end] = len(junctions)
            junctions.append({
                "ID": len(junctions),
                "coordinate": list(end),
                "element_id": element_id,
                "obj_type": obj_type,
                "obj_name": obj_name
            })

        lines.append({
            "ID": len(lines),
            "point": start,
            "direction": direction,
            "element_id": element_id,
            "obj_type": obj_type,
            "obj_name": obj_name
        })


# END HELPER FUNCTIONS

def extract_junctions_and_lines():
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True) # global coordinates
    global total_vertices
    global junctions_dict
    object_types = ["IfcWall", "IfcDoor", "IfcWindow", "IfcColumn", "IfcSpace"]
    for obj_type in object_types:
        try:
            objects = ifc_model.by_type(obj_type)
            for obj in objects:
                if obj.Representation:
                    try:
                        # Create geometry for the object
                        shape = ifcopenshell.geom.create_shape(settings, obj)
                        mesh = shape.geometry

                        # Extract vertices
                        vertices = mesh.verts  # Flat list of coordinates [x1, y1, z1, x2, y2, z2, ...]
                        edges = mesh.edges
                        _process_mesh(vertices,edges, obj.id(), obj_type, obj.Name)
                        #vertices = [(vertices[i], vertices[i + 1], vertices[i + 2]) for i in range(0, len(vertices), 3)]

                        #print(f"{obj_type} {obj.GlobalId} has {len(vertices)} vertices:")
                        #num_vertices += len(vertices)
                        #for vertex in vertices:
                            #print(f"  - {vertex}")
                    except Exception as e:
                        print(f"⚠ Could not process {obj_type} {obj.GlobalId}: {e}")
                else:
                    print(f"{obj_type} {obj.GlobalId} has no representation.")
        except Exception as e:
            # Print the exception message
            print(f"⚠ Exception occurred: {e}")
        
    print(f"✅ Extracted {len(junctions)} junctions and {len(lines)} lines.")
    print(f"total vertices = {total_vertices}")


def extract_planes():
    print("\n🔍 🏛️ Extracting planes from IfcWall, IfcSlab, and IfcDoor...")

    type_map = {
        "IfcWall": "wall"
    }

    found_planes = 0

    for entity_type in type_map.keys():
        for element in ifc_model.by_type(entity_type):

            placement = getattr(element, "ObjectPlacement", None)
            element_id = element.id()
            element_name = getattr(element, "Name", None)
            print(f"🔹 Processing {entity_type} ID {element_id}")
            shape = element.Representation
            if not shape:
                print(f"⚠ No representation for {entity_type} ID {element_id}")
                continue

            candidate_faces = []

            for rep in shape.Representations:
                for item in rep.Items:
                    def process_faces(faces, coords=None):
                        print("len faces: ", len(faces))
                        print("WALL")
                        for face in faces:
                            # Check for duplicate faces
                            
                            if face in candidate_faces:
                                print("Duplicate face found, skipping...")
                                print(face)
                                continue

                            print(face)
                            # Handle IfcIndexedPolygonalFace (IFC4)
                            if face.is_a("IfcIndexedPolygonalFace") and coords is not None:
                                # face.CoordIndex is 1-based
                                indices = [idx - 1 for idx in face.CoordIndex]
                                face_points = [coords[i] for i in indices]
                                print("face_points:")
                                for i, pt in enumerate(face_points):
                                    print(f"  {i}: [{pt[0]:10.3f}, {pt[1]:10.3f}, {pt[2]:10.3f}]")
                                global_tf = get_global_transform(placement)
                                face_points = [
                                    (global_tf @ np.append(pt, 1))[:3]
                                    for pt in face_points
                                ]
                                if len(face_points) >= 3:
                                    normal = _compute_normal(face_points[0], face_points[1], face_points[2])
                                    offset = np.dot(normal, face_points[0])
                                    area = compute_face_area(face_points)
                                    candidate_faces.append({
                                        "normal": normal,
                                        "offset": offset,
                                        "area": area,
                                        "points": face_points,
                                    })
                                continue

                            # Handle faces with Bounds (e.g., IfcFace)
                            if hasattr(face, "Bounds"):
                                for bound in face.Bounds:
                                    loop = bound.Bound
                                    if loop.is_a("IfcPolyLoop"):
                                        global_tf = get_global_transform(placement)
                                        face_points = [
                                            (global_tf @ np.append(pt.Coordinates, 1))[:3]
                                            for pt in loop.Polygon
                                        ]
                                        if len(face_points) >= 3:
                                            normal = _compute_normal(face_points[0], face_points[1], face_points[2])
                                            offset = np.dot(normal, face_points[0])
                                            area = compute_face_area(face_points)
                                            candidate_faces.append({
                                                "normal": normal,
                                                "offset": offset,
                                                "area": area,
                                                "points": face_points,
                                            })


                    if item.is_a("IfcFacetedBrep") and hasattr(item, "Outer"):
                        process_faces(item.Outer.CfsFaces)

                    elif item.is_a("IfcSurfaceModel") and hasattr(item, "Faces"):
                        process_faces(item.Faces)

                    elif item.is_a("IfcIndexedPolyCurve") and hasattr(item, "Faces"):
                        process_faces(item.Faces)

                    elif item.is_a("IfcPolygonalFaceSet") and hasattr(item, "Faces"): # For IFC4
                        # Get the coordinates array (CoordList)
                        coords = [np.array(coord) for coord in item.Coordinates.CoordList]
                        process_faces(item.Faces, coords=coords)

                    elif item.is_a("IfcMappedItem"):
                        mapped_items = item.MappingSource.MappedRepresentation.Items
                        for mapped in mapped_items:
                            if mapped.is_a("IfcFacetedBrep") and hasattr(mapped, "Outer"):
                                process_faces(mapped.Outer.CfsFaces)
                            elif mapped.is_a("IfcSurfaceModel") and hasattr(mapped, "Faces"):
                                process_faces(mapped.Faces)

            # Filter and select the vertical faces as wall planes
            if entity_type == "IfcWall":
                vertical_faces = [
                    face for face in candidate_faces if is_vertical(face["normal"])
                ]
                """best_face = max(vertical_faces, key=lambda f: f["area"], default=None)
                if best_face:
                    plane_data = {
                        "ID": len(planes),
                        "element_id": element_id,
                        "type": type_map.get(entity_type, "unknown"),
                        "normal": best_face["normal"],
                        "offset": best_face["offset"],
                    }
                    planes.append(plane_data)
                    found_planes += 1"""
                for face in vertical_faces:
                    plane_data = {
                        "ID": len(planes),
                        "element_id": element_id,
                        "type": type_map.get(entity_type, "unknown"),
                        "normal": face["normal"],
                        "offset": face["offset"],
                        "element_name": element_name,
                    }
                    planes.append(plane_data)
                    found_planes += 1
    print(f"✅ Extracted {found_planes} planes.")

def is_vertical(normal, tolerance=1e-2):
    """Check if a face is vertical (normal vector is perpendicular to Z-axis)."""
    normal = np.array(normal) / np.linalg.norm(normal)  # Normalize the vector
    print(abs(normal[2]))
    return abs(normal[2]) < tolerance  # Z-component should be close to 0

def is_horizontal(normal, tolerance=1e-6):
    """Check if a face is horizontal (normal vector is aligned with Z-axis)."""
    normal = np.array(normal) / np.linalg.norm(normal)  # Normalize the vector
    return abs(abs(normal[2]) - 1) < tolerance  # Z-component should be close to ±1

def compute_face_area(points):
    """Compute the area of a 3D polygon using triangulation."""
    if len(points) < 3:
        return 0.0
    area = 0.0
    origin = points[0]
    for i in range(1, len(points) - 1):
        v1 = points[i] - origin
        v2 = points[i + 1] - origin
        cross = np.cross(v1, v2)
        area += np.linalg.norm(cross) / 2.0
    return area  # Z-component of the cross product

def is_line_on_plane(line, plane, tol=0.05):
    start_point = np.array(line["point"])
    end_point = start_point + np.array(line["direction"])
    normal = np.array(plane["normal"])
    offset = plane["offset"]
    print("st:", start_point, "end:", end_point, "dir:", line["direction"])
    print("start_point distance: ", np.abs(np.dot(normal, start_point) - offset))
    print("end_point distance: ", np.abs(np.dot(normal, end_point) - offset))
    return (
        np.isclose(np.abs(np.dot(normal, start_point) - offset), 0, atol=tol) and
        np.isclose(np.abs(np.dot(normal, end_point) - offset), 0, atol=tol)
    )

def find_matching_wall_plane(lines_for_opening, wall_planes, tol=0.05):
    for plane in wall_planes:
        if all(is_line_on_plane(line, plane, tol) for line in lines_for_opening):
            return plane
    return None
def extract_semantics():
    print("\n🔍 Extracting semantics and assigning planes...")

    # Step 1: Map IFC Wall IDs to Plane IDs
    wall_to_plane_map = {}
    for plane in planes:
        if "element_id" in plane:
            wall_to_plane_map.setdefault(plane["element_id"], []).append(plane["ID"])

    print("wp_map: ",wall_to_plane_map)
    #raise ValueError("DEBUG: Wall to Plane Map created, now extracting semantics.")

    count = 0

    # Step 2: Extract Room Assignments via IfcRelSpaceBoundary
    space_plane_map = {}
    for boundary in ifc_model.by_type("IfcRelSpaceBoundary"):
        if boundary.RelatingSpace and boundary.RelatedBuildingElement:
            room_id = boundary.RelatingSpace.id()
            wall_id = boundary.RelatedBuildingElement.id()
            if wall_id in wall_to_plane_map:
                space_plane_map.setdefault(room_id, []).extend(wall_to_plane_map[wall_id]) # flatten the list of planes
    # Make entries unique per room_id
    for room_id in space_plane_map:
        space_plane_map[room_id] = list(set(space_plane_map[room_id]))
    # After building space_plane_map
    room_plane_counts = []
    for room_id, plane_lists in space_plane_map.items():
        # Flatten the list in case wall_to_plane_map[wall_id] is a list of plane IDs
        flat_planes = [item for sublist in plane_lists for item in (sublist if isinstance(sublist, list) else [sublist])]
        room_plane_counts.append((room_id, len(flat_planes)))

    # Sort by count descending
    room_plane_counts.sort(key=lambda x: x[1], reverse=True)

    for room_id, count in room_plane_counts:
        print(f"Room ID {room_id}: {count} entries")
    print("1272: ",space_plane_map.get(1272, []))
    print("space_plane_map: ", space_plane_map)
    
    # Add Rooms to semantics
    for space in ifc_model.by_type("IfcSpace"):
        room_id = space.id()
        assigned_planes = space_plane_map.get(room_id, [])
        print("space.name: ", space.Name)
        semantics.append({
            "ID": count,
            "planeID": assigned_planes,
            "type": "office",
            "element_id": room_id
        })
        if(space.Name == "R-143"):
            print("semantics: ", semantics[-1])
            #raise ValueError("DEBUG: Found R-143 room, check assigned planes.")
        count += 1
        print(f"🏠 Room {room_id} → Assigned Planes: {assigned_planes}")
    #raise ValueError("DEBUG: Space Plane Map created, now assigning planes to rooms.")

    # Step 3: Map IfcOpeningElement to Walls via IfcRelVoidsElement
    opening_to_wall_map = {}  # {opening_id: wall_id}
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    # Get all wall openings
    for rel_void in ifc_model.by_type("IfcRelVoidsElement"):
        opening = rel_void.RelatedOpeningElement
        wall = rel_void.RelatingBuildingElement
        if opening and wall and wall.is_a("IfcWall"):
            opening_to_wall_map[opening.id()] = wall.id() # map opening to wall
            print(f"🛠 DEBUG: IfcOpeningElement {opening.id()} is voiding IfcWall {wall.id()}")

        # Get the wall planes from the opening's wall ID
        if opening.id() in opening_to_wall_map:
            wall_id = opening_to_wall_map[opening.id()]
            if wall_id in wall_to_plane_map:
                base_plane_ids = wall_to_plane_map[wall_id] # list of plane IDs of the wall
            else:
                raise ValueError(f"Wall ID {wall_id} has no associated planes.")
        else:
            raise ValueError(f"Opening ID {opening_id} not found in opening_to_wall_map.")
        print("element_id of WALL: ", planes[base_plane_ids[0]]["element_id"])
        print("name of WALL: ", planes[base_plane_ids[0]]["element_name"])
        # Get geometry of the void (cut in the wall)
        shape = ifcopenshell.geom.create_shape(settings, opening)
        coords = np.array(shape.geometry.verts).reshape(-1, 3) * 1000  # Convert to millimeters

        print(f"Opening (in wall {wall.GlobalId}): {opening.GlobalId}")
        print("Opening vertex sample:", coords)
        # idea: the shortest axis of the opening will be the space between the two doors/windows
        # Get min and max for each axis
        min_x, min_y, min_z = np.min(coords, axis=0)
        max_x, max_y, max_z = np.max(coords, axis=0)

        x_diff = max_x - min_x
        y_diff = max_y - min_y
        z_diff = max_z - min_z
        # Assume x_diff, y_diff, z_diff are already defined
        diffs = {'x': x_diff, 'y': y_diff, 'z': z_diff}
        min_axis = min(diffs, key=diffs.get)
        min_value = diffs[min_axis]
        print(f"Minimum difference is {min_value} along {min_axis}-axis")
        # Find the index for the min_axis
        axis_idx = {'x': 0, 'y': 1, 'z': 2}[min_axis]
        other_axes = [i for i in range(3) if i != axis_idx]
        # Create the 8 bounding box (hull) coordinates (will be the coordinates of the two doors)
        hull_coords = np.array([
            [min_x, min_y, min_z],
            [min_x, min_y, max_z],
            [min_x, max_y, min_z],
            [min_x, max_y, max_z],
            [max_x, min_y, min_z],
            [max_x, min_y, max_z],
            [max_x, max_y, min_z],
            [max_x, max_y, max_z],
        ])
        print("hull_coords: ", hull_coords)
        # Get the two coordinates: one with min value, one with max value along min_axis
        min_axis_val = np.min(hull_coords[:, axis_idx])
        max_axis_val = np.max(hull_coords[:, axis_idx])

        # Get all coordinates of the two doors/windows
        door_window_coords_min = hull_coords[np.isclose(hull_coords[:, axis_idx], min_axis_val)]
        door_window_coords_max = hull_coords[np.isclose(hull_coords[:, axis_idx], max_axis_val)]

        # Sort by the other two axes for a closed loop order
        def loop_order(door_window, other_axes):
            # Project to the two main axes
            coords_2d = door_window[:, other_axes]
            hull = ConvexHull(coords_2d)
            loop_indices = list(hull.vertices)
            print("Convex Hull loop indices:", loop_indices)
            for i in loop_indices:
                print(f"Coordinate {i}: {door_window_coords_min[i]}")
            return loop_indices
        # Get the loop order for the two doors/windows
        loop_indices_min = loop_order(door_window_coords_min, other_axes)
        loop_indices_max = loop_order(door_window_coords_max, other_axes)
        # Sort the coordinates by the loop order
        door_window_coords_min = door_window_coords_min[loop_indices_min]
        door_window_coords_max = door_window_coords_max[loop_indices_max]
        direction_min = []
        direction_max = []
        for i, door_window in enumerate([door_window_coords_min, door_window_coords_max]):
            for j, coord in enumerate(door_window):
                # calculate for each coordinate the direction to the next coordinate
                if j < len(door_window) - 1:
                    next_coord = door_window[j + 1]
                    direction = np.array(next_coord) - np.array(coord)
                else:
                    # For the last coordinate, calculate direction to the first coordinate
                    next_coord = door_window[0]
                    direction = np.array(next_coord) - np.array(coord)
                if i == 0:
                    direction_min.append(direction.tolist())
                else:
                    direction_max.append(direction.tolist())
        # search wall plane for the coordinates of the two doors/windows
        wall_plane_min = None
        wall_plane_max = None
        for candidate_plane_id in base_plane_ids:
            print(" \ncandidate_plane_id: ", candidate_plane_id)
            candidate_plane = planes[candidate_plane_id]
            print("info: ", candidate_plane)
            correct_plane_min = True
            correct_plane_max = True
            for vert_idx in range(len(door_window_coords_min)):
                print("check min")
                if not is_line_on_plane({"point": door_window_coords_min[vert_idx], "direction": direction_min[vert_idx]}, candidate_plane, tol=10):
                    correct_plane_min = False
                print("check max")
                if not is_line_on_plane({"point": door_window_coords_max[vert_idx], "direction": direction_max[vert_idx]}, candidate_plane, tol=10):
                    correct_plane_max = False
            if correct_plane_min:
                wall_plane_min = candidate_plane
            if correct_plane_max:
                wall_plane_max = candidate_plane
        print("wall_plane_min: ", wall_plane_min)
        print("wall_plane_max: ", wall_plane_max )
        if wall_plane_min is None or wall_plane_max is None:
            raise ValueError(f"Could not find wall plane for opening {opening.GlobalId} with coordinates {door_window_coords_min} and {door_window_coords_max}.")
                
        print("door_window_coords_min:", door_window_coords_min)
        print("door_window_coords_max:", door_window_coords_max)
        # Create junctions for the two doors/windows
        global junctions_dict
        global junctions
        global lines
        plane_id_max = None
        plane_id_min = None
        for i, door_window in enumerate([door_window_coords_min, door_window_coords_max]):

            # Assign a different obj_name and element_id for min and max
            if i == 0:
                obj_name = f"{opening.Name}_min"
                element_id = int(str(opening.id()) + "00001")
            else:
                obj_name = f"{opening.Name}_max"
                element_id = int(str(opening.id()) + "00002")
            for j, coord in enumerate(door_window):
                coord_tuple = tuple(coord)
                if coord_tuple not in junctions_dict:
                    # Default obj_type
                    obj_type = "IfcOpeningElement"
                    # Try to find the correct obj_type from door_window_coordinates by closest door/window coordinate
                    closest_dw = None
                    min_dist = float('inf')
                    # assign obj_type based on the closest door/window coordinate
                    for dw in door_window_coordinates:
                        for dw_coord in dw["coordinates"]:
                            dist = np.linalg.norm(np.array(coord) - np.array(dw_coord))
                            if dist < min_dist:
                                min_dist = dist
                                closest_dw = dw
                    if closest_dw is not None:
                        obj_type = closest_dw["obj_type"]
                    junctions_dict[coord_tuple] = len(junctions)
                    junctions.append({
                        "ID": len(junctions),
                        "coordinate": list(coord),
                        "element_id": element_id,
                        "obj_type": obj_type,
                        "obj_name": obj_name
                    })
                    lines.append({
                        "ID": len(lines),
                        "point": list(coord),
                        "direction": direction_min[j] if i == 0 else direction_max[j],
                        "element_id": element_id,
                        "obj_type": obj_type,
                        "obj_name": obj_name
                    })
                    #print the added junction
                    print("Newly added junction:", junctions[-1])
                    print("Newly added line:", lines[-1])

            if obj_type == "IfcDoor":
                type_value = "door"
            elif obj_type == "IfcWindow":
                type_value = "window"
            else:
                raise ValueError(f"Unknown object type: {obj_type}")

            if i == 0:
                plane_id_min = len(planes)
            elif i == 1:
                plane_id_max = len(planes)
            
            planes.append({
            "ID": len(planes),
            "element_id": element_id,  # Associate this unique plane with the door/window element
            "type": type_value,
            "normal": wall_plane_min["normal"] if i == 0 else wall_plane_max["normal"],
            "offset": wall_plane_min["offset"] if i == 0 else wall_plane_max["offset"],
            "element_name": obj_name,
            "wall_id" : wall_id,
            "wall_plane_id": wall_plane_min["ID"] if i == 0 else wall_plane_max["ID"],
            })
            print("Newly added plane:", planes[-1])
        # find corresponding room by going trough each semantic element and checking if the wall_id is in the planeID list
        list_of_room_ids = []
        print("base_plane_ids: ", base_plane_ids)
        for semantic in semantics:
            print("semantic: ", semantic)
            if all(plane_id in semantic["planeID"] for plane_id in base_plane_ids) and semantic["type"] == "office":
                room_id = semantic["element_id"]
                list_of_room_ids.append(room_id)
        print("list_of_room_ids: ", list_of_room_ids)

        # take junctions of each ifcspace
        from collections import defaultdict
        junctions_for_each_space = defaultdict(list)
        
        for junction in junctions:
            if junction["element_id"] in list_of_room_ids:
                junctions_for_each_space[junction["element_id"]].append(junction)

        # Take ifcspace that has the closest distance to door
        closest_space_max = None
        closest_space_min = None
        closest_distance_door_max = float('inf')
        closest_distance_door_min = float('inf')
        for room_id, room_junctions in junctions_for_each_space.items():
            for junction in room_junctions:
                for i, door_window in enumerate([door_window_coords_max, door_window_coords_min]):
                    for coord in door_window:
                        dist = np.linalg.norm(np.array(junction["coordinate"]) - np.array(coord))
                        if dist < closest_distance_door_max and i == 0:  # For max door/window
                            closest_distance_door_max = dist
                            closest_space_max = room_id

                        if dist < closest_distance_door_min and i == 1:
                            temp = closest_space_min
                            closest_space_min = room_id  # For min door/window
                            if closest_space_min == closest_space_max:
                                closest_space_min = temp
                                continue
                            closest_distance_door_min = dist
        
        # Add Doors to semantics

        if obj_type == "IfcDoor":
            type_value = "door"
        elif obj_type == "IfcWindow":
            type_value = "window"
        else:
            raise ValueError(f"Unknown object type: {obj_type}")
        semantics.append({
            "ID": closest_space_max,
            "planeID": [plane_id_max],
            "type": type_value
        })

        semantics.append({
            "ID": closest_space_min,
            "planeID": [plane_id_min],
            "type": type_value
        })

        print("semantics: ", semantics[-2])
        print("semantics: ", semantics[-1])

    return


def generate_matrices():
    print("\n🔍 Generating matrices...")

    if not planes:
        print("❌ No planes found, skipping matrix generation.")
        return [], []

    planeLineMatrix = np.zeros((len(planes), len(lines)), dtype=int)
    lineJunctionMatrix = np.zeros((len(lines), len(junctions)), dtype=int)
    for p_idx, plane in enumerate(planes):
        show = True
        for l_idx, line in enumerate(lines):
            # Skip lines that do not belong to the same object as the plane
            wall_check = False

            dont_skip = False
            # Dont skip if plane is from door and line is from the wall that belongs to the door
            if wall_check:
                if plane["type"] == "wall" and line["obj_type"] == "IfcDoor":
                    if plane["ID"] == line["wall_plane_id"]:
                        dont_skip = True
                        #wall_check = True

            if "element_id" in line:
                if isinstance(line["element_id"], list):
                    match = plane["element_id"] in line["element_id"]
                else:
                    match = line["element_id"] == plane["element_id"]
                if  not match and not dont_skip:
                    continue

            # Check if the line is on the plane
            start_point = np.array(line["point"])
            end_point = start_point + np.array(line["direction"])
            normal = np.array(plane["normal"])
            offset = plane["offset"]

            normal_magnitude = np.linalg.norm(normal)

            
            relative_tolerance = None
            if (line["obj_type"] == "IfcDoor")  and specific_ifc:
                #if "second_plane" in line and "second_plane" in plane:
                    #pass
                    #print("second plane found")
                    #print("distances for door: ", np.abs(np.dot(normal, start_point) - offset), "end:", np.abs(np.dot(normal, end_point) - offset))
                relative_tolerance = 10 
            elif line["obj_type"] == "IfcWindow" and specific_ifc:
                relative_tolerance = 0.002
            else:
                relative_tolerance = 1e-8
            if show and "wall_id" in plane and wall_check:
                #print("start_point: ", start_point)
                #print("end_point: ", end_point)
                #print("normal: ", normal)
                #print("offset: ", offset)
                #print("np.dot(normal, start_point): ", np.dot(normal, start_point- offset))
                #print("np.isclose(np.dot(normal, end_point):", np.dot(normal, end_point- offset))
                #show = False
                print(np.abs(np.dot(normal, start_point) - offset))
                print(np.abs(np.dot(normal, end_point) - offset))
            if np.isclose(np.abs(np.dot(normal, start_point) - offset), 0, atol=relative_tolerance) and np.isclose(np.abs(np.dot(normal, end_point) - offset), 0, atol=relative_tolerance): # THIS MIGHT BE PROBLEMATIC TO USE CERTAIN DISTANCES TO CHECK FOR ADJACENT LINES TO PLANES
                #raise ValueError
                planeLineMatrix[p_idx, l_idx] = 1 # line belongs to this plane
            """elif line["obj_type"] == "IfcDoor" and plane["ID"] == 90:
                print("Door distances")
                print("planeID: ", plane["ID"])
                print("obj_name: ", line["obj_name"])
                print(np.abs(np.dot(normal, start_point) - offset))
                print(np.abs(np.dot(normal, end_point) - offset))"""
            """
            value = np.dot(normal, start_point) - offset
            if 0.1 < value < 10:
                print("value between 0 and 10:", value)
            """

    # Iterate through each line and find the closest junctions 
    for l_idx, line in enumerate(lines):
        point = np.array(line["point"])
        end_point = point + np.array(line["direction"])

        # Initialize variables to store the closest junctions and their distances
        closest_junction_to_point = None
        closest_junction_to_end_point = None
        min_distance_to_point = float('inf')
        min_distance_to_end_point = float('inf')

        only_junction_with_second = False
        if "second_plane" in line:
            only_junction_with_second = True

        for j_idx, junction in enumerate(junctions):
            # Skip junctions that do not belong to the same object as the line
            if "second_plane" not in junction and only_junction_with_second:
                #print("Skipping junction without second_plane")
                continue

            junction_coord = np.array(junction["coordinate"])

            # Calculate distances to point and end_point
            distance_to_point = np.linalg.norm(junction_coord - point)
            distance_to_end_point = np.linalg.norm(junction_coord - end_point)

            # Update the closest junction to point           
            if distance_to_point < min_distance_to_point:
                min_distance_to_point = distance_to_point
                closest_junction_to_point = j_idx
            
            # Update the closest junction to end_point
            if distance_to_end_point < min_distance_to_end_point:
                min_distance_to_end_point = distance_to_end_point
                closest_junction_to_end_point = j_idx
        #print("min_distance_to_point:",min_distance_to_point) # these should be very small
        #print("min_distance_to_end_point:",min_distance_to_end_point)
        # Overwrite the start and end coordinates of the line with the closest junction coordinates
        if closest_junction_to_point is not None:
            line["point"] = junctions[closest_junction_to_point]["coordinate"]
        if closest_junction_to_end_point is not None:
            line["direction"] = (
                np.array(junctions[closest_junction_to_end_point]["coordinate"]) -
                np.array(line["point"])
            ).tolist()
        # Set the matrix entries for the closest junctions
        if closest_junction_to_point is not None:
            lineJunctionMatrix[l_idx, closest_junction_to_point] = 1
        if closest_junction_to_end_point is not None:
            lineJunctionMatrix[l_idx, closest_junction_to_end_point] = 1
    return planeLineMatrix.tolist(), lineJunctionMatrix.tolist()
def add_doors_to_wall_planematrix(planeLineMatrix):
    #for each door plane, add their lines to the entries of the wall plane they are associated with in the planeLineMatrix
    global semantics
    global lines
    global planes
    for p_idx, plane in enumerate(planes):
        if plane["type"] == "door":
            wall_plane_id = plane["wall_plane_id"]
            door_plm = planeLineMatrix[plane["ID"]]
            wall_plm = planeLineMatrix[wall_plane_id]
            planeLineMatrix[wall_plane_id] = [max(door_plm[i], wall_plm[i]) for i in range(len(door_plm))]
    return planeLineMatrix

def add_floor_planes_for_offices(planeLineMatrix): # CEILING IS SOMEHOW NOT VISIBLE IN THE 3D MODEL
    """
    For each office, create a new horizontal floor plane and ceiling plane bounded by its space lines.
    Update semantics, planes, and planeLineMatrix accordingly.
    """
    global semantics
    global lines
    global planes
    next_plane_id = max(p["ID"] for p in planes) + 1 if planes else 0
    for semantic in semantics:
        if semantic["type"] == "office" and semantic["planeID"]:
            # Find all lines associated with this office semantic
            space_line_ids = [line["ID"] for line in lines if line["obj_type"] == "IfcSpace" and line["element_id"] == semantic["element_id"]]
            # Debug: print all space lines for this semantic
            print(f"\n--- Space lines for semantic (office) ID {semantic['ID']} (element_id={semantic['element_id']}):")
            for idx in space_line_ids:
                print(lines[idx])
            floor_height = float('inf')  # Initialize to a very high value
            ceiling_height = -float('inf')  # Initialize to a very low value
            for boundary_line_idx in space_line_ids: # get the lowest/highest Z coordinate of the boundary lines
                l = lines[boundary_line_idx]
                z1 = l["point"][2]
                z2 = l["point"][2] + l["direction"][2]
                floor_height = min(floor_height, z1, z2)
                ceiling_height = max(ceiling_height, z1, z2)
            # Add new floor and ceiling plane
            floor_plane = {
                "ID": next_plane_id,
                "type": "floor",
                "normal": [0, 0, 1],
                "offset": floor_height
            }
            print(f"floor_height: {floor_height}")
            ceiling_plane = {
                "ID": next_plane_id + 1,
                "type": "ceiling",
                "normal": [0, 0, -1],
                "offset": ceiling_height
            }
            # append to planes
            planes.append(floor_plane)
            planes.append(ceiling_plane)
            # Add a new row to planeLineMatrix for each plane
            planeLineMatrix.append([0 for _ in range(len(lines))])
            planeLineMatrix.append([0 for _ in range(len(lines))])
            # Assign lines to floor plane
            for i in space_line_ids:
                start_point = np.array(lines[i]["point"])
                end_point = start_point + np.array(lines[i]["direction"])
                # floor
                normal = np.abs(np.array(floor_plane["normal"]))
                offset = floor_plane["offset"]
                if np.isclose(np.abs(np.dot(normal, start_point) - offset), 0, atol=0.1) and np.isclose(np.abs(np.dot(normal, end_point) - offset), 0, atol=0.1):
                    planeLineMatrix[next_plane_id][lines[i]["ID"]] = 1
                #ceiling
                normal = np.abs(np.array(ceiling_plane["normal"]))
                offset = ceiling_plane["offset"]
                if np.isclose(np.abs(np.dot(normal, start_point) - offset), 0, atol=0.1) and np.isclose(np.abs(np.dot(normal, end_point) - offset), 0, atol=0.1):
                    planeLineMatrix[next_plane_id+1][lines[i]["ID"]] = 1

            # Update semantics to reference both new planes
            semantic["planeID"].append(next_plane_id)
            semantic["planeID"].append(next_plane_id+1)
            next_plane_id += 2
            
    return planeLineMatrix


def adjust_plane_normals():
    """
    For each plane, adjust its normal to point toward the closest IfcSpace (office).
    """
    global semantics, lines, planes

    import numpy as np

    # Step 1: Build a map from element_id to centroid of its IfcSpace lines
    space_centroids = {}
    for semantic in semantics:
        if semantic["type"] == "office":
            space_line_coords = []
            for line in lines:
                if line["obj_type"] == "IfcSpace" and line["element_id"] == semantic["element_id"]:
                    start = np.array(line["point"])
                    end = start + np.array(line["direction"])
                    space_line_coords.extend([start, end])
            if space_line_coords:
                centroid = np.mean(space_line_coords, axis=0)
                space_centroids[semantic["element_id"]] = centroid

    # Step 2: Adjust each plane's normal
    for plane in planes:
        if "normal" not in plane or "offset" not in plane:
            continue

        normal = np.array(plane["normal"])
        offset = plane["offset"]

        if np.linalg.norm(normal) == 0:
            continue

        # Compute a point on the plane
        point_on_plane = normal * (offset / np.linalg.norm(normal))

        # Step 3: Find closest centroid
        closest_centroid = None
        min_distance = float('inf')
        for centroid in space_centroids.values():
            distance = np.linalg.norm(centroid - point_on_plane)
            if distance < min_distance:
                min_distance = distance
                closest_centroid = centroid

        if closest_centroid is not None:
            to_centroid = closest_centroid - point_on_plane

            # Step 4: Flip normal if it points away
            if np.dot(normal, to_centroid) < 0:
                plane["normal"] = (-normal).tolist()

def fix_wall_normals_postprocess():
    """
    Post-process: for each wall element, ensure its two largest planes have opposite normals.
    """
    global planes
    import numpy as np
    from collections import defaultdict

    # Step 1: Group wall planes by element_id
    wall_plane_groups = defaultdict(list)
    for plane in planes:
        if plane.get("type") == "wall" and "element_id" in plane and "normal" in plane:
            wall_plane_groups[plane["element_id"]].append(plane)

    # Step 2: For each group, fix normals of top 2 planes
    for element_id, group in wall_plane_groups.items():
        if len(group) < 2:
            continue  # Need at least 2 to compare

        # Sort by area descending
        sorted_planes = sorted(group, key=lambda p: p.get("area", 0), reverse=True)
        p1, p2 = sorted_planes[0], sorted_planes[1]

        

        n1 = np.array(p1["normal"])
        n2 = np.array(p2["normal"])
        print(np.dot(n1,n2))

        # If normals are too similar, flip the second one
        if np.dot(n1, n2) >= 0.1:
            p1["normal"] = (-n2).tolist()


import numpy as np
from shapely.geometry import Point, Polygon

def generate_grid_samples(junction_coords, grid_spacing):
    """
    Generate a grid of sample points within the 2D footprint of the room.

    Args:
        junction_coords (list of tuples): List of 3D coordinates (x, y, z) of the room's junctions.
        grid_spacing (float): The spacing between grid points.

    Returns:
        list of tuples: List of 3D sample points (x, y, z).
    """
    # Step 1: Project 3D junction points to 2D (XY) to define room footprint
    polygon_2D = [(x, y) for (x, y, z) in junction_coords]

    # Step 2: Compute 2D bounding box of the room
    x_min = min(x for x, y in polygon_2D)
    x_max = max(x for x, y in polygon_2D)
    y_min = min(y for x, y in polygon_2D)
    y_max = max(y for x, y in polygon_2D)

    # Step 3: Estimate average height (z) for sample points
    avg_z = np.mean([z for x, y, z in junction_coords])

    # Step 4: Initialize empty list of sample points
    samples = []

    # Step 5: Create a Shapely polygon for point-in-polygon testing
    room_polygon = Polygon(polygon_2D)

    # Step 6: Loop through 2D grid inside the bounding box
    x = x_min
    while x <= x_max:
        y = y_min
        while y <= y_max:
            point_2D = Point(x, y)

            # Step 7: Check if point is inside the room polygon
            if room_polygon.contains(point_2D):
                sample_point = (x, y, avg_z)
                samples.append(sample_point)

            y += grid_spacing
        x += grid_spacing

    return samples

def adjust_plane_normals_works():
    """
    For each plane, adjust its normal to point toward the closest IfcSpace (office).
    """
    global semantics, lines, planes

    # For each plane, adjust normal to point toward closest IfcSpace centroid
    for plane in planes:
        plane_normal = np.array(plane["normal"])
        plane_offset = plane["offset"]
        plane_type = plane["type"]

        shortest_distance = float('inf')
        closest_element_id = None
        closest_junction = None
        for semantic in semantics:
            if semantic["type"] == "office":
                for junc in junctions:
                    if junc["obj_type"] == "IfcSpace" and junc["element_id"] == semantic["element_id"]:
                        point = np.array(junc["coordinate"])
                        distance = np.abs(np.dot(plane_normal, point) - plane_offset)
                        if distance < shortest_distance:
                            shortest_distance = distance
                            closest_element_id = semantic["element_id"]
                            closest_junction = junc


        junction_coords = np.array([j["coordinate"] for j in junctions if j["element_id"] == closest_element_id])

        """grid_spacing = 10
        grid_samples = generate_grid_samples(junction_coords, grid_spacing)

        closest_grid_sample = min(grid_samples, key=lambda sample: np.abs(np.dot(plane_normal, sample) - plane_offset))"""

        centroid = np.mean(junction_coords, axis=0)
        point_on_plane = plane_normal * plane_offset
        vector_to_centroid = centroid - point_on_plane
        plane_point = plane_normal * (plane_offset / np.linalg.norm(plane_normal))
        #vector_to_sample = np.array(closest_grid_sample) - point_on_plane
        # make sure the normals points in the same direction of the room centroid
        if np.dot(plane_normal, vector_to_centroid) < 0:
            plane["normal"] = (-plane_normal).tolist()
            plane["offset"] = (-plane_offset)
        else:
            if plane_type == "door" or plane_type == "window":
                # For doors and windows, we want the normal to point outward
                plane["normal"] = (-plane_normal).tolist()
                plane["offset"] = (-plane_offset)

def adjust_plane_normals_works2():
    """
    For each plane, adjust its normal to point toward the closest IfcSpace (office).
    """
    global semantics, lines, planes

    # For each plane, adjust normal to point toward closest IfcSpace centroid
    for plane in planes:
        plane_normal = np.array(plane["normal"])
        plane_offset = plane["offset"]


        shortest_distance = float('inf')
        second_shortest_distance = float('inf')
        second_closest_junction = None
        closest_element_id = None
        closest_junction = None
        for semantic in semantics:
            if semantic["type"] == "office":
                for junc in junctions:
                    if junc["obj_type"] == "IfcSpace" and junc["element_id"] == semantic["element_id"]:
                        point = np.array(junc["coordinate"])
                        distance = np.abs(np.dot(plane_normal, point) - plane_offset)
                        if distance < shortest_distance:
                            shortest_distance = distance
                            closest_element_id = semantic["element_id"]
                            closest_junction = junc
        for semantic in semantics:
            if semantic["type"] == "office":
                for junc in junctions:
                    if junc.get("element_id") == closest_element_id:
                        # Your code here, e.g.:
                        point = np.array(junc["coordinate"])
                        distance = np.abs(np.dot(plane_normal, point) - plane_offset)
                        if distance < second_shortest_distance and distance > 0.0:
                            second_shortest_distance = distance
                            second_closest_junction = junc
                        # ... do something with point ...
        print("shortest_distance: ", shortest_distance)

        #junction_coords = np.array([j["coordinate"] for j in junctions if j["element_id"] == closest_element_id])
        #centroid = np.mean(junction_coords, axis=0)
        point_on_plane = plane_normal * plane_offset
        vector_to_centroid = second_closest_junction["coordinate"] - point_on_plane
        plane_point = plane_normal * (plane_offset / np.linalg.norm(plane_normal))
        # make sure the normals points in the same direction of the room centroid
        if np.dot(plane_normal, vector_to_centroid) < 0:
            plane["normal"] = (-plane_normal).tolist()    
            plane["offset"] = (-plane_offset).tolist()  # Also flip the offset to match the new normal direction
def convert_linejunctionmatrix_to_BIMKIT(lnm):
    """
    Converts the lineJunctionMatrix to a BIMKIT-compatible format.

    Args:
        lnm (list): The lineJunctionMatrix to convert.

    Returns:
        list: The converted matrix in BIMKIT format.
    """
    # Initialize the BIMKIT matrix
    bimkit_matrix = []
    all_rows_valid = True

    # Track all junctions that are referenced by any line
    referenced_junctions = set()

    # Iterate through each row of the lineJunctionMatrix
    for row_idx, row in enumerate(lnm):
        # Find the indices of the junctions that are connected
        connected_junctions = [idx for idx, value in enumerate(row) if value == 1]
        # Track referenced junctions
        referenced_junctions.update(connected_junctions)
        # Check if this row has exactly two junctions
        if len(connected_junctions) != 2:
            print(f"❌ Row {row_idx} does not have exactly two junctions: {connected_junctions}")
            all_rows_valid = False
        # Append the row index and connected junctions to the BIMKIT matrix
        bimkit_matrix.append({"line": row_idx, "connected_junctions": connected_junctions})

    # Check for junctions with no lines
    total_junctions = len(lnm[0]) if lnm else 0
    all_junctions = set(range(total_junctions))
    unreferenced_junctions = all_junctions - referenced_junctions
    if unreferenced_junctions:
        print(f"⚠ Junctions with no lines: {sorted(unreferenced_junctions)}")
    else:
        print("✅ All junctions are referenced by at least one line.")

    if all_rows_valid:
        print("✅ All rows have exactly two junctions (two values of 1 per row).")
    else:
        print("❌ Not all rows have exactly two junctions.")

    # Save the BIMKIT matrix to a JSON file with proper indentation
    output_json_path = "/raid/USERDATA/ovrobi4y/SPVLOC/spvloc/2025-03-10_Stand_AA/spvloc_ifc_to_json/linejunctionMatrix.json"
    with open(output_json_path, 'w') as f:
        json.dump(bimkit_matrix, f, indent=4, separators=(',', ': '))  # Proper formatting

    return bimkit_matrix

def print_door_lines():
    """Print all lines associated with doors."""
    door_ids = set(door.id() for door in ifc_model.by_type("IfcDoor"))
    print("Door IDs:", door_ids)
    print("Lines associated with doors:")
    for line in lines:
        # If you store element association in the line dict, use that.
        # Otherwise, you may need to relate lines to doors by geometry or semantics.
        # Here, we check if the line's point or direction matches any door junction.
        # If you have a better mapping, adjust accordingly.
        for semantic in semantics:
            if semantic["type"] == "office":
                for plane_id in semantic["planeID"]:
                    # Find lines associated with this plane
                    for p_idx, plane in enumerate(planes):
                        if plane["ID"] == plane_id:
                            # Print all lines that belong to this plane
                            for l_idx, val in enumerate(planeLineMatrix[p_idx]):
                                if val == 1:
                                    print(lines[l_idx])
    print("Done printing door lines.")
def print_top_planes_by_line_count(planeLineMatrix, planes, top_n=130):
    """Prints the top N planes with the most and fewest lines on them."""
    # Count lines per plane
    line_counts = [(i, sum(row)) for i, row in enumerate(planeLineMatrix)]
    # Sort by line count
    line_counts_sorted = sorted(line_counts, key=lambda x: x[1], reverse=True)

    print(f"\nTop {top_n} planes with the MOST lines:")
    for idx, count in line_counts_sorted[:top_n]:
        plane_id = planes[idx]["ID"] if "ID" in planes[idx] else idx
        print(f"Plane {plane_id}: {count} lines")

    print(f"\nTop {top_n} planes with the FEWEST lines:")
    for idx, count in line_counts_sorted[-top_n:]:
        plane_id = planes[idx]["ID"] if "ID" in planes[idx] else idx
        print(f"Plane {plane_id}: {count} lines")
def print_lines_for_plane(planeLineMatrix, planes, lines, plane_id):
    """
    Prints the line IDs associated with the given plane_id.
    """
    # Find the index of the plane in the planes list
    plane_idx = next((i for i, p in enumerate(planes) if p["ID"] == plane_id), None)
    if plane_idx is None:
        print(f"Plane ID {plane_id} not found.")
        return

    line_ids = [lines[l_idx]["ID"] for l_idx, val in enumerate(planeLineMatrix[plane_idx]) if val == 1]
    print(f"Lines for plane {plane_id}: {line_ids}")
if __name__ == "__main__":
    
    extract_junctions_and_lines()
    extract_planes()
    extract_semantics()
    
    old_planeLineMatrix, lineJunctionMatrix = generate_matrices()
    planeLineMatrix = add_floor_planes_for_offices(old_planeLineMatrix) # Expand the floor plane by 500mm in X and Y
    planeLineMatrix = add_doors_to_wall_planematrix(planeLineMatrix) # Add door lines to the wall planes they are associated with
    adjust_plane_normals_works() # Adjust plane normals to point towards the closest IfcSpace (office)
    #fix_wall_normals_postprocess()

    # Fix `type` to remove "ifc" prefix
    for plane in planes:
        if "type" in plane:
            plane["type"] = plane["type"].replace("ifc", "")

    # Fix `type` for semantics elements (doors/windows)
    for semantic in semantics:
        if "type" in semantic:
            semantic["type"] = semantic["type"].replace("ifc", "")

    output_planes = [
        {k: v for k, v in plane.items() if k not in ["edges"]} for plane in planes
    ]
    output_lines_first = [
        {k: v for k, v in line.items() if k not in ["point_index", "end_point", "end_point_index", "obj_type"]} for line in lines
    ]
    output_lines = [
        {k: v for k, v in line.items()} for line in output_lines_first
    ]
    output_data = {
        "junctions": junctions,
        "lines": output_lines,
        "planes": output_planes, # Excluding element_id
        "semantics": semantics,
        "planeLineMatrix": planeLineMatrix,
        "lineJunctionMatrix": lineJunctionMatrix,
        "cuboids": cuboids,
        "manhattan": manhattan_structures,
    }
    rows = len(planeLineMatrix)  # Number of rows
    cols = len(planeLineMatrix[0]) if rows > 0 else 0  # Number of columns
    print(f"Shape of the list: ({rows}, {cols})")
    with open(output_json_path, 'w') as f:
        json.dump(output_data, f)

    # Extract edges from lineJunctionMatrix
    edges = []
    for row in lineJunctionMatrix:
        connected_junctions = [idx for idx, value in enumerate(row) if value == 1]
        if len(connected_junctions) == 2:  # Ensure the row represents a valid edge
            edges.append((connected_junctions[0], connected_junctions[1]))
    print(f"Extracted {len(edges)} edges from lineJunctionMatrix.")
    print_top_planes_by_line_count(planeLineMatrix, planes)
    convert_linejunctionmatrix_to_BIMKIT(lineJunctionMatrix) # TEMPORARY
    print_lines_for_plane(planeLineMatrix, planes, lines, 331)
    print(f"✅ JSON export complete: {output_json_path}")