import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.shape_builder
import json
import numpy as np

# File paths
ifc_file_path = "/raid/USERDATA/ovrobi4y/SPVLOC/spvloc/2025-03-10_Stand_AA/IC6_s3d_IFC2x3.ifc"
output_json_path = "/raid/USERDATA/ovrobi4y/SPVLOC/spvloc/2025-03-10_Stand_AA/annotation_3d.json"
junctions = []
junctions_dict = {}
lines = []
planes = []
semantics = []
cuboids = []
manhattan_structures = []
ifc_model = ifcopenshell.open(ifc_file_path)

# HELPER FUNCTIONS

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

def _save_junction_lines(points):
    if len(points) >= 2:
        for i in range(len(points) - 1):
            start, end = points[i], points[i + 1]

            start = tuple(float(coord) for coord in (start + (0.0,) * (3 - len(start))))
            end = tuple(float(coord) for coord in (end + (0.0,) * (3 - len(end))))

            direction = np.array(end) - np.array(start)

            if np.linalg.norm(direction) < 1e-6:
                print(f"⚠ Skipping degenerate line (ID {len(lines)}) - start and end points are the same.")
                continue
            
            if start not in junctions_dict:
                junctions_dict[start] = len(junctions)
                junctions.append({"ID": len(junctions), "coordinate": list(start)})

            if end not in junctions_dict:
                junctions_dict[end] = len(junctions)
                junctions.append({"ID": len(junctions), "coordinate": list(end)})

            new_point = list(start)
            new_direction = list(direction)

            exists = any(line["point"] == new_point and line["direction"] == new_direction for line in lines)

            if exists:
                continue

            lines.append({
                "ID": len(lines),
                "point": new_point,
                "direction": new_direction
            })


def _process_mesh(verts, edges):

        points = [(verts[i], verts[i+1], verts[i+2]) for i in range(0, len(verts), 3)]
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
    for product in ifc_model.by_type("IfcProduct"):
        if not product.Representation:
            continue
        try:
            settings = ifcopenshell.geom.settings()
            settings.set(settings.USE_WORLD_COORDS, True)
            shape = ifcopenshell.geom.create_shape(settings, product)
            mesh = shape.geometry
            _process_mesh(mesh.verts, mesh.edges)
        except Exception as e:
            # Print the exception message
            print(f"⚠ Exception occurred: {e}")
        
    print(f"✅ Extracted {len(junctions)} junctions and {len(lines)} lines.")


def extract_planes():
    print("\n🔍 Extracting planes from IfcWall, IfcSlab, and IfcDoor...")

    type_map = {
        "IfcWall": "wall",
        "IfcSlab": "floor",
    }

    found_planes = 0

    for entity_type in type_map.keys():
        for element in ifc_model.by_type(entity_type):
            element_id = element.id()
            print(f"🔹 Processing {entity_type} ID {element_id}")
            shape = element.Representation
            if not shape:
                print(f"⚠ No representation for {entity_type} ID {element_id}")
                continue

            for rep in shape.Representations:
                print(f"  🔹 Representation Type: {rep.RepresentationType}")

                for item in rep.Items:
                    def process_faces(faces):
                        for face in faces:
                            for bound in face.Bounds:
                                loop = bound.Bound
                                if loop.is_a("IfcPolyLoop"):
                                    face_points = [np.array(pt.Coordinates) for pt in loop.Polygon]
                                    if len(face_points) >= 3:
                                        normal = _compute_normal(face_points[0], face_points[1], face_points[2])
                                        offset = np.dot(normal, face_points[0])
                                        plane_data = {
                                            "ID": len(planes),
                                            "element_id": element_id,
                                            "type": type_map.get(entity_type, "unknown"),
                                            "normal": normal,
                                            "offset": offset
                                        }
                                        planes.append(plane_data)
                                        nonlocal found_planes
                                        found_planes += 1

                    if item.is_a("IfcFacetedBrep") and hasattr(item, "Outer"):
                        process_faces(item.Outer.CfsFaces)

                    elif item.is_a("IfcSurfaceModel") and hasattr(item, "Faces"):
                        process_faces(item.Faces)

                    elif item.is_a("IfcMappedItem"):
                        mapped_items = item.MappingSource.MappedRepresentation.Items
                        for mapped in mapped_items:
                            if mapped.is_a("IfcFacetedBrep") and hasattr(mapped, "Outer"):
                                process_faces(mapped.Outer.CfsFaces)
                            elif mapped.is_a("IfcSurfaceModel") and hasattr(mapped, "Faces"):
                                process_faces(mapped.Faces)

    print(f"✅ Extracted {found_planes} planes.")


def extract_semantics():
    print("\n🔍 Extracting semantics and assigning planes...")

    # Step 1: Map IFC Wall IDs to Plane IDs
    wall_to_plane_map = {plane["element_id"]: plane["ID"] for plane in planes if "element_id" in plane}
    print(f"🛠 DEBUG: Mapping IFC Wall IDs to Plane IDs: {wall_to_plane_map}")

    count = 0

    # Step 2: Extract Room Assignments via IfcRelSpaceBoundary
    space_plane_map = {}
    for boundary in ifc_model.by_type("IfcRelSpaceBoundary"):
        if boundary.RelatingSpace and boundary.RelatedBuildingElement:
            room_id = boundary.RelatingSpace.id()
            wall_id = boundary.RelatedBuildingElement.id()
            if wall_id in wall_to_plane_map:
                space_plane_map.setdefault(room_id, []).append(wall_to_plane_map[wall_id])

    # Add Rooms to semantics
    for space in ifc_model.by_type("IfcSpace"):
        room_id = space.id()
        assigned_planes = space_plane_map.get(room_id, [])
        semantics.append({
            "ID": count,
            "planeID": assigned_planes,
            "type": "room"
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
        for l_idx, line in enumerate(lines):
            point = np.array(line["point"])
            normal = np.array(plane["normal"])
            offset = plane["offset"]

            # Check if the point already lies on a plane
            #if np.any(planeLineMatrix[:, l_idx] == 1):
            #   continue

            normal_magnitude = np.linalg.norm(normal)

            relative_tolerance = 0.05 * normal_magnitude

            if np.isclose(np.dot(normal, point) - offset, 0, atol=10000000):
                planeLineMatrix[p_idx, l_idx] = 1
            
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

        # Set the matrix entries for the closest junctions
        if closest_junction_to_point is not None:
            lineJunctionMatrix[l_idx, closest_junction_to_point] = 1
        if closest_junction_to_end_point is not None:
            lineJunctionMatrix[l_idx, closest_junction_to_end_point] = 1

    return planeLineMatrix.tolist(), lineJunctionMatrix.tolist()

def check_empty_rows_in_planeLineMatrix(planeLineMatrix):
    """
    Checks each row of planeLineMatrix to see if it contains only 0 elements.
    Prints the row numbers of empty rows.

    Args:
        planeLineMatrix (list of lists or np.ndarray): The matrix to check.
    """
    for row_idx, row in enumerate(planeLineMatrix):
        if all(value == 0 for value in row):  # Check if all elements in the row are 0
            print(f"Row {row_idx} contains only 0 elements.")
def check_linejunctionMatrix(lineJunctionMatrix):
    """
    Checks each row of lineJunctionMatrix to see if it contains only 0 elements.
    Prints the row numbers of empty rows.

    Args:
        lineJunctionMatrix (list of lists or np.ndarray): The matrix to check.
    """
    for row_idx, row in enumerate(lineJunctionMatrix):
        next_junc = row[1]
        for elem in lineJunctionMatrix:
            if elem[0] == next_junc:
                next_junc = elem[1]
                break

def find_cycles_or_noncycles(data):
    from collections import defaultdict
    from tqdm import tqdm
    # Build undirected graph
    graph = defaultdict(list)
    for i, (a, b) in enumerate(data):
        graph[a].append((i, b))
        graph[b].append((i, a))

    visited_paths = []

    def dfs(path, visited_indices, visited_nodes, current_node, start_node):
        if current_node in visited_nodes:
            return

        if not graph[current_node]:
            # Dead end
            visited_paths.append(path[:])
            print("Dead end:", path)
            return

        dead_end = True
        for i, neighbor in graph[current_node]:
            if i not in visited_indices:
                visited_indices.add(i)
                visited_nodes.add(current_node)
                path.append((min(current_node, neighbor), max(current_node, neighbor)))

                dfs(path, visited_indices, visited_nodes, neighbor, start_node)

                path.pop()
                visited_indices.remove(i)
                visited_nodes.remove(current_node)
                dead_end = False

        if dead_end:
            visited_paths.append(path[:])
            print("Dead end:", path)

    # Add a progress bar for the outer loop
    print("len(data):", len(data))
    with tqdm(total=len(data), desc="Processing edges") as progress_bar:
        for i, (a, b) in enumerate(data):
            dfs([(a, b)], {i}, {a}, b, a)
            progress_bar.update(1)  # Update progress bar only for the outer loop


    print("\nNon-cyclic paths:")
    for path in visited_paths:
        print(path)
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

    output_data = {
        "junctions": junctions,
        "lines": lines,
        "planes": planes,
        "semantics": semantics,
        "planeLineMatrix": planeLineMatrix,
        "lineJunctionMatrix": lineJunctionMatrix,
        "cuboids": cuboids,
        "manhattan": manhattan_structures,
    }
    check_empty_rows_in_planeLineMatrix(planeLineMatrix)
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
    #edges = [(0,1),(1,2),(2,3),(0,3),(1,3),(0,2),(0,4), (4,5),(5,6)]
    find_cycles_or_noncycles(edges)

    print(f"✅ JSON export complete: {output_json_path}")