"""
Read stl file, find holes and fill them
Author: Dr.-Ing. Omar ELSAYED
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import trimesh

def fill_loop(boundary_linked_list):
    filled_faces = []
    num_vertices = len(boundary_linked_list) - 1
    mid = 1 + ((num_vertices) // 2)

    for i in range(0, mid - 1):
        filled_faces.append([boundary_linked_list[i], boundary_linked_list[i + 1], boundary_linked_list[num_vertices - i]])
        filled_faces.append([boundary_linked_list[i + 1], boundary_linked_list[num_vertices - i - 1], boundary_linked_list[num_vertices - i]])

    if (num_vertices) % 2 == 0:
        filled_faces.append([boundary_linked_list[mid - 1], boundary_linked_list[mid], boundary_linked_list[mid + 1]])

    return filled_faces

# Load the original mesh from STL file
original_mesh = trimesh.load('convex_hole.stl')

edges_sorted = original_mesh.edges_sorted
original_faces = original_mesh.faces

edges_unique, edge_counts = np.unique(edges_sorted, axis=0, return_counts=True)
boundary_edges = edges_unique[edge_counts == 1]

# Step 1: Build adjacency dictionary for boundary vertices
adjacency_dict = {}
for edge in boundary_edges:
    if edge[0] not in adjacency_dict:
        adjacency_dict[edge[0]] = []
    if edge[1] not in adjacency_dict:
        adjacency_dict[edge[1]] = []
    adjacency_dict[edge[0]].append(edge[1])
    adjacency_dict[edge[1]].append(edge[0])

# Step 2: Function to extract a loop from the adjacency dictionary
def extract_loop(start_vertex, adjacency_dict, visited):
    boundary_linked_list = []
    current_vertex = start_vertex
    previous_vertex = -1
    while True:
        if current_vertex in visited:
            break
        visited.add(current_vertex)
        boundary_linked_list.append(current_vertex)
        if current_vertex in adjacency_dict:
            neighbors = adjacency_dict[current_vertex]
            if neighbors:
                next_vertex = neighbors[0] if neighbors[0] != previous_vertex else neighbors[-1]
                if next_vertex == start_vertex:
                    boundary_linked_list.append(next_vertex)
                    visited.add(next_vertex) 
                    break
                previous_vertex = current_vertex
                current_vertex = next_vertex
                if current_vertex == start_vertex:  
                    visited.add(current_vertex)  
                    break
            else:
                break
        else:
            break
    
    if boundary_linked_list and boundary_linked_list[0] == boundary_linked_list[-1]:
        return boundary_linked_list[0:-1]
    else:
        return None

# Step 3: Detect multiple loops
visited = set()
loops = []

for vertex in adjacency_dict:
    if vertex not in visited:
        loop = extract_loop(vertex, adjacency_dict, visited)
        if loop:
            loops.append(loop)
        else:
            print(f"Loop starting at vertex {vertex} is not closed and cannot be filled.")

# Step 4: Fill each loop and collect faces
all_filled_faces = original_faces.tolist()


# Step 4: Fill each loop and collect faces
all_filled_faces = original_faces.tolist()
for loop in loops:
    filled_faces = fill_loop(loop)
    all_filled_faces.extend(filled_faces)
    
# Function to get coordinates from the original mesh
def get_coordinates(mesh):
    vertices = mesh.vertices
    return vertices[:, :2]  # Consider only x and y coordinates for 2D plotting

# Animate the filling process
fig, ax = plt.subplots()
ax.set_aspect('equal')

# Get coordinates from the original mesh
coordinates = get_coordinates(original_mesh)

# Initialize the plot with original faces
for face in original_faces:
    polygon = plt.Polygon([coordinates[face[0]], coordinates[face[1]], coordinates[face[2]]], closed=True, edgecolor='black', facecolor='lightblue')
    ax.add_patch(polygon)

# Draw boundary edges
for edge in boundary_edges:
    ax.plot([coordinates[edge[0], 0], coordinates[edge[1], 0]], [coordinates[edge[0], 1], coordinates[edge[1], 1]], 'k-')

# Plot vertices
ax.scatter(coordinates[:, 0], coordinates[:, 1], color='black')

ax.set_xticks([])
ax.set_yticks([])
ax.set_aspect('equal')
ax.set_title('Filling Mesh Animation')

# Calculate frames based on new faces added in each step
frames_per_step = []
num_original_faces = len(original_faces)
for loop in loops:
    filled_faces = fill_loop(loop)
    frames_per_step.append(len(filled_faces))

# Update function for animation
def update(frame):
    ax.clear()
    ax.set_aspect('equal')
    
    # Draw original faces
    for face in original_faces:
        polygon = plt.Polygon([coordinates[face[0]], coordinates[face[1]], coordinates[face[2]]], closed=True, edgecolor='black', facecolor='lightblue')
        ax.add_patch(polygon)
    
    # Draw filled faces up to the current frame
    current_frame = 0
    for step_faces in frames_per_step:
        if frame >= current_frame + step_faces:
            current_frame += step_faces
        else:
            for i in range(current_frame, frame + 1):
                face = all_filled_faces[num_original_faces + i]
                polygon = plt.Polygon([coordinates[face[0]], coordinates[face[1]], coordinates[face[2]]], closed=True, edgecolor='black', facecolor='lightblue')
                ax.add_patch(polygon)
            break
    
    # Draw boundary edges
    for edge in boundary_edges:
        ax.plot([coordinates[edge[0], 0], coordinates[edge[1], 0]], [coordinates[edge[0], 1], coordinates[edge[1], 1]], 'k-')

    # Plot vertices
    ax.scatter(coordinates[:, 0], coordinates[:, 1], color='black')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect('equal')
    ax.set_title('Filling Mesh Animation')

    return ax

ani = animation.FuncAnimation(fig, update, frames=sum(frames_per_step), interval=2, repeat=False)
ani.save('./animation.gif', writer='imagemagick', fps=30)
