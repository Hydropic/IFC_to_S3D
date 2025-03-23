# --------------------------
# Extract Planes from Walls and Slabs
import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util.shape_builder
import json
import numpy as np

# File paths
ifc_file_path = "/raid/USERDATA/ovrobi4y/SPVLOC/spvloc/2025-03-10_Stand_AA/IC6_s3d_2x3.ifc"
output_json_path = "/raid/USERDATA/ovrobi4y/SPVLOC/spvloc/2025-03-10_Stand_AA/annotation_3d.json"

# Load IFC file
print("🔍 Loading IFC file...")
try:
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    ifc_model = ifcopenshell.open(ifc_file_path)
    vertices_list = []
    for element in ifc_model.by_type('IfcProduct'):
        if element.Representation:
            shape = ifcopenshell.geom.create_shape(settings, element)
            geometry = shape.geometry
            
            # Extract vertices and convert to world coordinates
            vertices = np.array(geometry.verts).reshape(-1, 3)
            print(f"Element {element.GlobalId}:")
            print(vertices)
            vertices_list.append(vertices)
    print("Ifcproduct: ",len(ifc_model.by_type('IfcProduct')))
    print("IfcPolyline: ",len(ifc_model.by_type('IfcPolyline')))
    #breakpoint()
    #settings.set(settings.USE_PYTHON_OPENCASCADE, True)
    print("✅ IFC file loaded successfully.")
    print("📌 IFC Schema:", ifc_model.wrapped_data.schema)
except Exception as e:
    print("❌ Error loading IFC file:", e)
    exit()

# Print all entity types in the IFC file
entity_types = set(entity.is_a() for entity in ifc_model)
print("\n📋 IFC Entity Types Found:", sorted(entity_types))

# Data structures
junctions = []
junctions_dict = {}
lines = []
planes = []
semantics = []
cuboids = []
manhattan_structures = []
# --------------------------
# Compute normal vector
def compute_normal(p1, p2, p3):
    v1 = np.array(p2) - np.array(p1)
    v2 = np.array(p3) - np.array(p1)
    normal = np.cross(v1, v2)
    norm_length = np.linalg.norm(normal)
    return (normal / norm_length).tolist() if norm_length != 0 else [0, 0, 1]

# --------------------------
# Extract Junctions and Lines
def extract_junctions_and_lines(): # CHANGES WERE MADE HERE TO CREATE MAKE UNIQUE LINES
    print("\n🔍 Extracting junctions and lines...")

    if "IfcPolyline" not in entity_types:
        print("⚠ Warning: No IfcPolyline found in IFC model!")
        return

    # Use a set to track unique lines
    unique_lines = set()
    idxpoly = 0
    shape_builder = ifcopenshell.util.shape_builder.ShapeBuilder(ifc_model)
    #print("polylines: ",ifc_model.by_type("IfcPolyline"))
    for polyline in vertices_list:
        """ 
        print(f"Processing polyline of type: {type(polyline)}")
        t = list(shape_builder.get_polyline_coords(polyline))
        print("t: ",t)
        """
        points = [tuple(pt) for pt in polyline]
        idxpoly += 1
        if len(points) >= 2:
            for i in range(len(points) - 1):
                start, end = points[i], points[i + 1]
                print("poly: ",idxpoly," start coordinate: ",start)
                # Ensure 3D coordinates (fill missing dimensions with 0.0)
                start = tuple(float(coord) for coord in (start + (0.0,) * (3 - len(start))))
                end = tuple(float(coord) for coord in (end + (0.0,) * (3 - len(end))))

                # Compute direction vector
                direction = np.array(end) - np.array(start)

                # Validate the line has a valid direction (not zero vector)
                if np.linalg.norm(direction) < 1e-6:
                    print(f"⚠ Skipping degenerate line (ID {len(lines)}) - start and end points are the same.")
                    continue

                # Ensure each junction is stored uniquely
                if start not in junctions_dict:
                    junctions_dict[start] = len(junctions)
                    junctions.append({"ID": len(junctions), "coordinate": list(start)})

                if end not in junctions_dict:
                    junctions_dict[end] = len(junctions)
                    junctions.append({"ID": len(junctions), "coordinate": list(end)})

                # Create a tuple representing the line (sorted to handle both directions)
                line_tuple = tuple(sorted((start, end)))

                # Add the line only if it's not already in the set of unique lines
                if line_tuple not in unique_lines:
                    unique_lines.add(line_tuple)
                    lines.append({
                        "ID": len(lines),
                        "point": list(start),
                        "direction": list(direction)  # Convert to Python list
                    })

    print(f"✅ Extracted {len(junctions)} junctions and {len(lines)} lines.")


# --------------------------
# Extract Planes from Walls and Slabs
def extract_planes():
    print("\n🔍 Extracting planes from IfcWall and IfcSlab...")

    found_planes = 0
    for entity_type in ["IfcWall", "IfcSlab"]:
        if entity_type in entity_types:
            for element in ifc_model.by_type(entity_type):
                element_id = element.id()  # Internal use only
                print(f"🔹 Processing {entity_type} ID {element_id}")
                shape = element.Representation
                if not shape:
                    print(f"⚠ No representation for {entity_type} ID {element_id}")
                    continue
                for rep in shape.Representations:
                    print(f"  🔹 Representation Type: {rep.RepresentationType}")
                    for item in rep.Items:
                        if item.is_a("IfcFacetedBrep") or item.is_a("IfcMappedItem") or item.is_a("IfcSurfaceModel"):
                            for shell in item.Outer.CfsFaces:
                                face_points = []
                                for bound in shell.Bounds:
                                    loop = bound.Bound
                                    if loop.is_a("IfcPolyLoop"):
                                        face_points = [np.array(pt.Coordinates) for pt in loop.Polygon]
                                        if len(face_points) >= 3:
                                            normal = compute_normal(face_points[0], face_points[1], face_points[2])
                                            offset = np.dot(normal, face_points[0])

                                            # Store element_id internally, remove in final output
                                            plane_data = {
                                                "ID": len(planes),
                                                "element_id": element_id,  # Internal use only
                                                "type": "wall" if entity_type == "IfcWall" else "floor",
                                                "normal": normal,
                                                "offset": offset
                                            }

                                            planes.append(plane_data)
                                            found_planes += 1

    print(f"✅ Extracted {found_planes} planes.")


# --------------------------
# Extract Semantics (Rooms, Doors, Windows)
def extract_semantics():
    print("\n🔍 Extracting semantics and assigning planes...")

    # Step 1: Map IFC Wall IDs to Plane IDs
    wall_to_plane_map = {plane["element_id"]: plane["ID"] for plane in planes if "element_id" in plane}
    print(f"🛠 DEBUG: Mapping IFC Wall IDs to Plane IDs: {wall_to_plane_map}")

    # Step 2: Extract Room Assignments via IfcRelSpaceBoundary
    space_plane_map = {}  # {room_id: [plane_ids]}
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
            "ID": room_id,
            "planeID": assigned_planes,
            "type": "room"
        })
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
                    #"element_id": related_element.id(),  # Associate this unique plane with the door/window element
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
            "ID": door_id,
            "planeID": assigned_planes,
            "type": "door"
        })
        print(f"🚪 Door {door_id} → Assigned Planes: {assigned_planes}")

    # Add Windows to semantics
    for window in ifc_model.by_type("IfcWindow"):
        window_id = window.id()
        assigned_planes = window_to_planes.get(window_id, [])
        semantics.append({
            "ID": window_id,
            "planeID": assigned_planes,
            "type": "window"
        })
        print(f"🪟 Window {window_id} → Assigned Planes: {assigned_planes}")

    print(f"✅ Extracted {len(semantics)} semantic elements.")

# --------------------------
# Run Extractions
extract_junctions_and_lines()
extract_planes()
extract_semantics()

# --------------------------
# Generate Matrices
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
            if np.isclose(np.dot(normal, point) - offset, 0, atol=1e-3): # This might also include more or less lines than in actuality
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
        #print(f"Line {l_idx}: Start Junction {closest_junction_to_point} (Distance: {min_distance_to_point}, Coordinates: {point}), End Junction {closest_junction_to_end_point} (Distance: {min_distance_to_end_point}, Coordinates: {end_point})")
    # Verification step #1: Check if each plane is assigned to exactly 1 line start point and end point
    """ 
    for j_idx in range(len(junctions)):
        start_points = np.sum(lineJunctionMatrix[:, j_idx] == 1)
        end_points = np.sum(lineJunctionMatrix[:, j_idx] == 1)
        if start_points != 1 or end_points != 1:
            print(f"❌ Junction {j_idx} is assigned to {start_points} start points and {end_points} end points.")
        else:
            print(f"✅ Junction {j_idx} is correctly assigned to 1 start point and 1 end point.")
    """
    #test = np.array([[1,1,0,0],[0,1,1,0],[0,0,1,1],[1,0,0,1]])
    #verify_closed_loops(test)
    #find_unconnected_junctions(lineJunctionMatrix)

    return planeLineMatrix.tolist(), lineJunctionMatrix.tolist()

def find_unconnected_junctions(lineJunctionMatrix):
    unconnected_junctions = []
    for j_idx in range(lineJunctionMatrix.shape[1]):
        connections = np.sum(lineJunctionMatrix[:, j_idx])
        if connections == 0:
            unconnected_junctions.append(j_idx)
    if unconnected_junctions:
        print(f"❌ Unconnected junctions found: {unconnected_junctions}")
    else:
        print("✅ All junctions are connected.")

def verify_closed_loops(lineJunctionMatrix):
    from collections import defaultdict, deque

    # Create an adjacency list for the graph
    adjacency_list = defaultdict(list)
    for l_idx, row in enumerate(lineJunctionMatrix):
        junction_indices = np.where(row == 1)[0]
        if len(junction_indices) == 2:
            j1, j2 = junction_indices
            adjacency_list[j1].append(j2)
            adjacency_list[j2].append(j1)

    # Debug information for adjacency list
    print("\n🔍 Adjacency List:")
    for junction, neighbors in adjacency_list.items():
        print(f"Junction {junction}: Neighbors {neighbors}")

    # Function to check if the graph is a single closed loop
    def is_single_closed_loop(adjacency_list):
        visited = set()
        start_node = next(iter(adjacency_list))
        queue = deque([(start_node, None)])
        dead_ends = []

        while queue:
            node, parent = queue.popleft()
            if node in visited:
                continue
            visited.add(node)
            neighbors = adjacency_list[node]
            if len(neighbors) != 2:
                dead_ends.append(node)
            for neighbor in neighbors:
                if neighbor != parent:
                    queue.append((neighbor, node))

        # Check if all nodes are visited and the graph is connected
        is_closed_loop = len(visited) == len(adjacency_list) and all(len(neighbors) == 2 for neighbors in adjacency_list.values())
        return is_closed_loop, dead_ends

    is_closed_loop, dead_ends = is_single_closed_loop(adjacency_list)

    if is_closed_loop:
        print("✅ The lineJunctionMatrix forms a single closed loop.")
    else:
        print("❌ The lineJunctionMatrix does not form a single closed loop.")
        if dead_ends:
            print(f"🔍 Dead ends found at junctions: {dead_ends}")



# Generate matrices
planeLineMatrix, lineJunctionMatrix = generate_matrices()

# Remove "element_id" before exporting JSON
output_planes = [
    {k: v for k, v in plane.items() if k != "element_id"} for plane in planes
]

# Fix `type` to remove "ifc" prefix
for plane in output_planes:
    if "type" in plane:
        plane["type"] = plane["type"].replace("ifc", "")

# Fix `type` for semantics elements (doors/windows)
for semantic in semantics:
    if "type" in semantic:
        semantic["type"] = semantic["type"].replace("ifc", "")

# Structure final output
output_data = {
    "junctions": junctions,
    "lines": lines,
    "planes": output_planes,  # Excluding element_id
    "semantics": semantics,
    "planeLineMatrix": planeLineMatrix,
    "lineJunctionMatrix": lineJunctionMatrix,
    "cuboids": cuboids,
    "manhattan": manhattan_structures,
}

with open(output_json_path, 'w') as f:
    json.dump(output_data, f)

print(f"✅ JSON export complete: {output_json_path}")



def analyze_voids():
    print("\n🔍 Analyzing IfcRelVoidsElement (Wall Openings)...")

    wall_opening_map = {}  # {wall_id: [opening_ids]}
    for rel_void in ifc_model.by_type("IfcRelVoidsElement"):
        opening = rel_void.RelatedOpeningElement
        wall = rel_void.RelatingBuildingElement
        if opening and wall and wall.is_a("IfcWall"):
            wall_id = wall.id()
            opening_id = opening.id()
            wall_opening_map.setdefault(wall_id, []).append(opening_id)
    for wall_id, openings in wall_opening_map.items():
        if len(openings) > 1:
            print(f"🧱 Wall {wall_id} has {len(openings)} openings: {openings}")
    print("✅ Finished analyzing openings in walls.")

# Run the analysis
analyze_voids()

