# MorphoNeuro

MorphoNeuro is a Blender add-on designed for importing, constructing, and analyzing 3D neurosurgical corridors, either using neuronavigation coordinates or 3D meshes (surface or volumes). Meshes can be built using surface scanning techniques or traditional medical imaging (CT, MRI). MorphoNeuro has been tested on surface scans.

## Requirements
- Blender 4.2 or newer (older versions might be compatible, but have not been tested)

## Installation
1. Download the add-on as ZIP
2. In Blender, go to Edit > Preferences > Add-ons
3. Click "Install..." and select the downloaded ZIP file
4. Enable the add-on by checking the box next to "Object: MorphoNeuro"
5. MorphoNeuro can now be accessed from the 3D View sidebar (press N to open) under the "MorphoNeuro" tab.

## Setup
All calculations will be in millimeters, except volumes (cubic centimeters). If you're working with 3D reconstructions, ensure you set your Blender units accordingly (Scene properties > Unit scale = 0.001, Unit length = millimeters).

## Features

### 1. Modeling
Creates a volumetric surgical corridor either importing coordinates from external neuronavigation systems, selecting vertices manually on 3D meshes, or modeling automatically based on the mesh geometrical boundaries.

#### Neuronavigation
1.  Provide the path to a coordinate file (`.json`, `.dat`, or `.txt`)
    *   `.json`: Expected format is a list of lists, where each inner list represents a group of points (e.g., `[[apex1, base1a, base1b], [apex2, base2a, base2b]]`).
    *   `.dat`/`.txt`: Stryker's Intellect annotations.dat files are supported out of the box. Others (eg. Medtronic) might need manual preparsing to JSON (see above).
2.  Hit "Create"
- **Options**
    **Grouping:** Coordinates are grouped based on the `Grouping Number`. Useful if you're exporting a single file with a lot of coordinates and you have a fixed number of coordinates per target. The first point in each group becomes the apex, the rest form the base. Set to 0 to disable. Parsed `.dat`/`.txt` data is automatically saved as a `.json` file alongside the original.

#### Manual Selection to Mesh
1.  Enter Edit Mode on a mesh
2.  Select vertices sequentially by pressing Shift + RMB. The *first* selected vertex becomes the apex, and subsequent selections form the base polygon in order.
3.  Hit "Create"
- **Options**
    *   `Create Normalized Plane`: If checked (default: True), creates an additional normalized mesh.
        *   `Projection Distance`: Specifies the distance from the apex for the normalized base (default: 100 units).
        *   `Use Spherical Geometry`: If checked, uses a spherical geometry to project the base vertices; otherwise, the centroid is used.

#### Automatic corridor
1.  Enter Edit mode on a mesh
2.  Select a single vertex (surgical target)
3.  Hit "Create"
- **Options**
    *   `Num Vertices`: Target number of vertices for the base of the final corridor after decimation (default: 1000). The internal calculation uses a higher precision (10000 vertices) before decimation.
    *   `Distance from Apex`: The distance from the apex at which the base plane is normalized to (default: 100 units).

### 2. Polishing
Tools to normalize corridors to a fixed height for comparison between approaches/subjects, and to average coordinates of repeated measurements.

#### Normalize meshes
Creates a new mesh named `"{Original Name} normalized"` using the original apex and the new projected base coordinates.

1.  Select one or more surgical corridors
2.  Hit "Run"

- **Options**
    *   `Projection Distance`: The target distance from the apex for the normalized base (default: 100 units).
    *   `Projection Method`:
        *   `SPHERICAL`: Projects base vertices onto a sphere centered at the apex with radius `Projection Distance`.
        *   `LEAST SQUARES`: Fits a plane to the original base vertices using a least squares cost function. Projects the base vertices onto a flat plane located `Projection Distance` away from the apex along the plane's normal (ensuring the normal points generally from apex towards base).
        *   `CENTROID DIRECTION`: Calculates the centroid of the original base vertices and projects the base vertices onto a flat plane perpendicular to the apex-centroid direction, located `Projection Distance` from the apex along that direction.

#### Average coordinates
Generates a single mesh representing the geometric average of two or more selected meshes, point by point. Assumes meshes have similar topology.

1.  Select two or more surgical corridors
2.  Hit "Run"

- **Options**
    *   `Align meshes`: Perform ICP registration before averaging (default: False).

### 3. Analysis
#### Volume, area, angle calculations
1.  Select one or more surgical corridors
2.  Hit one of the icons. The result will be both displayed in the 3D view and in the calculation history panel with an option to copy or delete each value. 


#### Vertex Distances
Calculates the Euclidean distance between corresponding vertices of a reference (active) mesh and one or more target (selected) meshes.

1.  Select two or more meshes
2.  Hit "Calculate"

- **Options**
    *   `Align meshes`: Perform ICP registration before distance calculation (default: True).
    *   `Export to CSV`: Save the results to a CSV file (default: True).
        *   `CSV File Path`: Location to save the CSV file (default: `//vertex_distances.csv`).

### Pairwise Distances
Compares the distances between pairs of vertices within two selected meshes (Mesh A and Mesh B).

1.  Select two meshes
2.  Hit "Calculate"

- **Options**
    *   `Export to CSV`: Save the results to a CSV file (default: True).
        *   `CSV File Path`: Location to save the CSV file (default: `//pairwise_distances.csv`).

### 4. Other tools
#### ICP Registration
Aligns one or more selected "source" meshes to an active "target" mesh using the Iterative Closest Point (ICP) algorithm.

1.  Select the surgical corridor(s) to be moved. The active object will be the target every other corridor will be registered to.
2.  Hit "Run"

- **Options**
    *   `Max Iterations`: Maximum number of ICP refinement steps (default: 50).
    *   `Convergence Threshold`: Stops if the improvement in the mean squared error between iterations is less than this value (default: 0.0001).
    *   `Use PCA Pre-alignment`: Perform initial alignment using Principal Component Analysis (default: True). Useful if working with neuronavigation meshes with similar geometry for a quick alignment using.
    *   `Use Normal Compatibility`: Consider surface normals during point matching (default: False).
    *   `Normal Weight`: If `Use Normal Compatibility` is True, this weight (0.0-1.0) balances distance vs. normal alignment. 0 ignores normals, 1 uses only normals (default: 0.3).
