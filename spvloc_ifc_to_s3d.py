import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.shape_builder
import json
import numpy as np
import sys

# File paths
ifc_file_path = "/raid/USERDATA/ovrobi4y/SPVLOC/spvloc/2025-03-10_Stand_AA/centerediIFCmm.ifc"
output_json_path = "/raid/USERDATA/ovrobi4y/SPVLOC/spvloc/2025-03-10_Stand_AA/zind/scene_00100/annotation_3d.json"
junctions = []
junctions_dict = {}
lines = []
planes = []
semantics = []
cuboids = []
total_vertices = 0
manhattan_structures = []
ifc_model = ifcopenshell.open(ifc_file_path)

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
        return [0.0, 0.0, 1.0]


def _process_mesh(verts, edges):
    scale = 1000
    points = [(verts[i]*scale, verts[i+1]*scale, verts[i+2]*scale) for i in range(0, len(verts), 3)]
    print(points)
    global junctions
    global junctions_dict

    for i in range(0, len(edges), 2):
        start = points[edges[i]]
        end = points[edges[i + 1]]

        if start not in junctions_dict:
            junctions_dict[start] = len(junctions)
            junctions.append({
                "ID": len(junctions),
                "coordinate": list(start)
            })

        if end not in junctions_dict:
            junctions_dict[end] = len(junctions)
            junctions.append({
                "ID": len(junctions),
                "coordinate": list(end)
            })

        direction = tuple(e - s for s, e in zip(start, end))

        lines.append({
            "ID": len(lines),
            "point": start,
            "direction": direction
        })


# EMD HELPER FUNCTIONS

def extract_junctions_and_lines():
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True) 
    print("junctions should be empty: ", junctions)
    print("lines should be empty: ", lines)
    global total_vertices
    global junctions_dict
    object_types = ["IfcWall", "IfcDoor", "IfcWindow", "IfcSlab","IfcColumn", "IfcSpace"]
    for obj_type in object_types:
        try:
            objects = ifc_model.by_type(obj_type)
            for obj in objects:
                print("obj_name", obj.Name)
                if obj.Representation:
                    try:
                        # Create geometry for the object
                        shape = ifcopenshell.geom.create_shape(settings, obj)
                        mesh = shape.geometry

                        # Extract vertices
                        vertices = mesh.verts  # Flat list of coordinates [x1, y1, z1, x2, y2, z2, ...]
                        edges = mesh.edges
                        _process_mesh(vertices,edges)
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
    print("\n🔍 Extracting planes from IfcWall, IfcSlab, and IfcDoor...")

    type_map = {
        "IfcWall": "wall",
        "IfcSlab": "floor",
    }

    found_planes = 0

    for entity_type in type_map.keys():
        for element in ifc_model.by_type(entity_type):
            placement = getattr(element, "ObjectPlacement", None)
            element_id = element.id()
            print(f"🔹 Processing {entity_type} ID {element_id}")
            shape = element.Representation
            if not shape:
                print(f"⚠ No representation for {entity_type} ID {element_id}")
                continue

            candidate_faces = []

            for rep in shape.Representations:
                print(f"  🔹 Representation Type: {rep.RepresentationType}")

                for item in rep.Items:
                    print(item)
                    
                    def process_faces(faces):
                        
                        for face in faces:
                            for bound in face.Bounds:
                                loop = bound.Bound
                                if loop.is_a("IfcPolyLoop"):
                                    global_tf = get_global_transform(placement)
                                    face_points = [
                                        (global_tf @ np.append(pt.Coordinates, 1))[:3]
                                        for pt in loop.Polygon
                                    ]
                                    print("##### FACE_POINTS: ",face_points)
                                    if len(face_points) >= 3:
                                        normal = _compute_normal(face_points[0], face_points[1], face_points[2])
                                        print("##### NORMAL: ",normal)
                                        offset = np.dot(normal, face_points[0])
                                        print("##### OFFSET: ",offset)
                                        area = compute_face_area(face_points)
                                        print("##### AREA: ",area)
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
                    
                    elif item.is_a("IfcPolygonalFaceSet") and hasattr(item, "Faces"):
                        process_faces(item.Faces)

                    elif item.is_a("IfcMappedItem"):
                        mapped_items = item.MappingSource.MappedRepresentation.Items
                        for mapped in mapped_items:
                            if mapped.is_a("IfcFacetedBrep") and hasattr(mapped, "Outer"):
                                process_faces(mapped.Outer.CfsFaces)
                            elif mapped.is_a("IfcSurfaceModel") and hasattr(mapped, "Faces"):
                                process_faces(mapped.Faces)

            # Filter and select the best face
            if entity_type == "IfcWall":
                # Select the largest vertical face
                vertical_faces = [
                    face for face in candidate_faces if is_vertical(face["normal"])
                ]
                print(vertical_faces)
                best_face = max(vertical_faces, key=lambda f: f["area"], default=None)

            elif entity_type == "IfcSlab":
                # Select the largest horizontal face
                horizontal_faces = [
                    face for face in candidate_faces if is_horizontal(face["normal"])
                ]
                best_face = max(horizontal_faces, key=lambda f: f["area"], default=None)

            else:
                best_face = None

            if best_face:
                plane_data = {
                    "ID": len(planes),
                    "element_id": element_id,
                    "type": type_map.get(entity_type, "unknown"),
                    "normal": best_face["normal"],
                    "offset": best_face["offset"],
                }
                planes.append(plane_data)
                found_planes += 1

    print(f"✅ Extracted {found_planes} planes.")

def is_vertical(normal, tolerance=1e-6):
    """Check if a face is vertical (normal vector is perpendicular to Z-axis)."""
    normal = np.array(normal) / np.linalg.norm(normal)  # Normalize the vector
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
    
def extract_semantics():
    print("\n🔍 Extracting semantics and assigning planes...")

    # Step 1: Map IFC Wall IDs to Plane IDs
    wall_to_plane_map = {plane["element_id"]: plane["ID"] for plane in planes if "element_id" in plane}
    print("wall_to_plane_map : ", wall_to_plane_map)
    print(f"🛠 DEBUG: Mapping IFC Wall IDs to Plane IDs: {wall_to_plane_map}")

    count = 0

    # Step 2: Extract Room Assignments via IfcRelSpaceBoundary
    space_plane_map = {}
    for boundary in ifc_model.by_type("IfcRelSpaceBoundary"):
        if boundary.RelatingSpace and boundary.RelatedBuildingElement:
            room_id = boundary.RelatingSpace.id()
            wall_id = boundary.RelatedBuildingElement.id()      
            if wall_id in wall_to_plane_map:   
                plane = wall_to_plane_map[wall_id]
                if plane not in space_plane_map.setdefault(room_id, []):
                    space_plane_map[room_id].append(plane)

    # Add Rooms to semantics
    for space in ifc_model.by_type("IfcSpace"):
        room_id = space.id()
        assigned_planes = space_plane_map.get(room_id, [])
        semantics.append({
            "ID": count,
            "planeID": assigned_planes,
            "type": "office"
        })
        count += 1
        print(f"🏠 Room {room_id} → Assigned Planes: {assigned_planes}")

    # Step 3: Map IfcOpeningElement to Walls via IfcRelVoidsElement
    opening_to_wall_map = {}  # {opening_id: wall_id}
    for rel_void in ifc_model.by_type("IfcRelVoidsElement"):
        opening = rel_void.RelatedOpeningElement
        wall = rel_void.RelatingBuildingElement
        if opening and wall and wall.is_a("IfcWall"):
            opening_to_wall_map[opening.id()] = wall.id()
            print(f"🛠 DEBUG: IfcOpeningElement {opening.id()} is voiding IfcWall {wall.id()}")

    # Step 4: Assign Doors and Windows to unique planes for each opening
    # We'll create a mapping keyed by (element_id, opening_id) so that multiple openings on the same wall get separate plane assignments.
    door_window_to_plane_map = {}  # {(door/window_id, opening_id): plane_id}

    for fills in ifc_model.by_type("IfcRelFillsElement"):
        if not hasattr(fills, "RelatingOpeningElement"):
            print(f"⚠ Warning: IfcRelFillsElement {fills.id()} has no RelatingOpeningElement!")
            continue
        opening = fills.RelatingOpeningElement
        if opening is None:
            continue
        opening_id = opening.id()
        related_element = fills.RelatedBuildingElement  # Should be IfcDoor or IfcWindow
        if related_element is None:
            continue
        print(f"🔍 Checking IfcRelFillsElement: {fills.id()} → Opening: {opening_id} → Related Element: {related_element.id()} ({related_element.is_a()})")
        if opening_id in opening_to_wall_map:
            wall_id = opening_to_wall_map[opening_id]
            if wall_id in wall_to_plane_map:
                base_plane_id = wall_to_plane_map[wall_id]
                # Create a unique plane for this door/window based on the opening.
                new_plane_id = len(planes)
                # We copy the normal and offset from the wall's plane.
                planes.append({
                    "ID": new_plane_id,
                    "element_id": related_element.id(),  # Associate this unique plane with the door/window element
                    "type": related_element.is_a().lower(),
                    "normal": planes[base_plane_id]["normal"],
                    "offset": planes[base_plane_id]["offset"]
                })
                door_window_to_plane_map[(related_element.id(), opening_id)] = new_plane_id
                print(f"✅ {related_element.is_a()} {related_element.id()} → Opening {opening_id} → Wall {wall_id} → Unique Plane {new_plane_id}")
            else:
                print(f"⚠ Warning: Wall {wall_id} has no mapped plane!")
        else:
            print(f"⚠ Warning: Opening {opening_id} is not linked to any wall!")

    print(f"🛠 DEBUG: Final Door and Window Plane Map: {door_window_to_plane_map}")

    # Step 5: Group unique plane assignments by door/window ID
    door_to_planes = {}  # {door_id: [plane_ids]}
    window_to_planes = {}  # {window_id: [plane_ids]}
    for (elem_id, opening_id), plane_id in door_window_to_plane_map.items():
        if any(door.id() == elem_id for door in ifc_model.by_type("IfcDoor")):
            door_to_planes.setdefault(elem_id, []).append(plane_id)
        if any(window.id() == elem_id for window in ifc_model.by_type("IfcWindow")):
            window_to_planes.setdefault(elem_id, []).append(plane_id)

    # Add Doors to semantics
    for door in ifc_model.by_type("IfcDoor"):
        door_id = door.id()
        assigned_planes = door_to_planes.get(door_id, [])
        semantics.append({
            "ID": count,
            "planeID": assigned_planes,
            "type": "door"
        })
        print(f"🚪 Door {door_id} → Assigned Planes: {assigned_planes}")
        count += 1

    # Add Windows to semantics
    for window in ifc_model.by_type("IfcWindow"):
        window_id = window.id()
        assigned_planes = window_to_planes.get(window_id, [])
        semantics.append({
            "ID": count,
            "planeID": assigned_planes,            
            "type": "window"

        })
        print(f"🪟 Window {window_id} → Assigned Planes: {assigned_planes}")
        count += 1

    print(f"✅ Extracted {len(semantics)} semantic elements.")


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
            start_point = np.array(line["point"])
            end_point = start_point + np.array(line["direction"])
            normal = np.array(plane["normal"])
            offset = plane["offset"]

            normal_magnitude = np.linalg.norm(normal)
            print(np.dot(normal, start_point) - offset)

            relative_tolerance = 1e-4 #0.01 * normal_magnitude
            if show:
                #print("start_point: ", start_point)
                #print("end_point: ", end_point)
                #print("normal: ", normal)
                print("offset: ", offset)
                #print("np.dot(normal, start_point): ", np.dot(normal, start_point- offset))
                #print("np.isclose(np.dot(normal, end_point):", np.dot(normal, end_point- offset))
                show = False
            if np.isclose(np.dot(normal, start_point) - offset, 0, atol=relative_tolerance) and np.isclose(np.dot(normal, end_point) - offset, 0, atol=relative_tolerance): # THIS MIGHT BE PROBLEMATIC TO USE CERTAIN DISTANCES TO CHECK FOR ADJACENT LINES TO PLANES
                print(plane)
                print(line)
                #raise ValueError
                planeLineMatrix[p_idx, l_idx] = 1

    # Iterate through each line and find the closest junctions (unused due to distances being possibly unreliable)
    for l_idx, line in enumerate(lines):
        point = np.array(line["point"])
        end_point = point + np.array(line["direction"])

        # Initialize variables to store the closest junctions and their distances
        closest_junction_to_point = None
        closest_junction_to_end_point = None
        min_distance_to_point = float('inf')
        min_distance_to_end_point = float('inf')

        for j_idx, junction in enumerate(junctions):
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

    # Iterate through each row of the lineJunctionMatrix
    for row_idx, row in enumerate(lnm):
        # Find the indices of the junctions that are connected
        connected_junctions = [idx for idx, value in enumerate(row) if value == 1]
        # Check if this row has exactly two junctions
        if len(connected_junctions) != 2:
            print(f"❌ Row {row_idx} does not have exactly two junctions: {connected_junctions}")
            all_rows_valid = False
        # Append the row index and connected junctions to the BIMKIT matrix
        bimkit_matrix.append({"line": row_idx, "connected_junctions": connected_junctions})
    if all_rows_valid:
        print("✅ All rows have exactly two junctions (two values of 1 per row).")
    else:
        print("❌ Not all rows have exactly two junctions.")

    # Save the BIMKIT matrix to a JSON file with proper indentation
    output_json_path = "/raid/USERDATA/ovrobi4y/SPVLOC/spvloc/2025-03-10_Stand_AA/spvloc_ifc_to_json/linejunctionMatrix.json"
    with open(output_json_path, 'w') as f:
        json.dump(bimkit_matrix, f, indent=4, separators=(',', ': '))  # Proper formatting

    return bimkit_matrix
if __name__ == "__main__":
    
    extract_junctions_and_lines()
    extract_planes()
    extract_semantics()

    planeLineMatrix, lineJunctionMatrix = generate_matrices()

    # Fix `type` to remove "ifc" prefix
    for plane in planes:
        if "type" in plane:
            plane["type"] = plane["type"].replace("ifc", "")

    # Fix `type` for semantics elements (doors/windows)
    for semantic in semantics:
        if "type" in semantic:
            semantic["type"] = semantic["type"].replace("ifc", "")

    output_planes = [
        {k: v for k, v in plane.items() if k not in ["element_id", "edges"]} for plane in planes
    ]
    output_lines = [
        {k: v for k, v in line.items() if k not in ["point_index", "end_point", "end_point_index", ]} for line in lines
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
            
    convert_linejunctionmatrix_to_BIMKIT(lineJunctionMatrix) # TEMPORARY
    print(f"✅ JSON export complete: {output_json_path}")