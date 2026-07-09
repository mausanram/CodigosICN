import numpy as np
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from astropy.io import fits
import skimage as sk
import numpy.ma as ma

# --- 1. Load and Process FITS Data ---
path = "../images/AGO2024/proc_skp_m-009_microchip_T_170__NSAMP_324_NROW_650_NCOL_700_EXPOSURE_0_img_100.fits"
hdu_list = fits.open(path)
extension = 1

# Extract and mask active area data
active_area = hdu_list[extension - 1].data[:, :550]
active_area_mask = sk.measure.label(active_area >= np.max(active_area), connectivity=2)
data = ma.masked_array(active_area, mask=(active_area_mask > 0))

# Clean up masked data for 3D plotting
filled_data = data.filled(0)

# Downsampling prevents performance lags
downsample_factor = 1
filled_data = filled_data[::downsample_factor, ::downsample_factor]
shape_y, shape_x = filled_data.shape

# --- 2. Setup pyqtgraph 3D Window ---
app = pg.mkQApp("GLMeshItem 3D Bars Example")
w = gl.GLViewWidget()
w.show()
w.setWindowTitle('pyqtgraph example: 3D FITS Mesh Bars')
w.setCameraPosition(distance=300, elevation=100, azimuth=90)
w.pan(shape_y / 2, shape_x / 2, 0)

# --- 3. Grid Environment Setup ---
gx = gl.GLGridItem(); gx.rotate(90, 0, 1, 0); gx.translate(-10, 0, 10); w.addItem(gx)
gy = gl.GLGridItem(); gy.rotate(90, 1, 0, 0); gy.translate(0, -10, 10); w.addItem(gy)
gz = gl.GLGridItem(); w.addItem(gz)

# --- 4. Custom 3D Mesh Generation ---
num_bars = shape_y * shape_x
verts = np.empty((num_bars * 8, 3))
faces = np.empty((num_bars * 12, 3), dtype=np.uint32)
face_colors = np.empty((num_bars * 12, 4))

# Set up the colormap
cm = pg.colormap.get('plasma')
# cm.reverse()

# --- NEW: Logarithmic Transformation ---
# np.log1p handles data points safely even if filled values are 0.
# We apply a scaling multiplier afterwards to ensure the 3D bars are tall enough to see.
h_matrix = np.log1p(filled_data) * 10.0  

# # Precompute min/max bounds for proper color scaling
# h_matrix = filled_data * 0.00001

# h_min, h_max = h_matrix.min(), np.percentile(h_matrix, 99.0)
h_min, h_max = h_matrix.min(), h_matrix.max()
h_range = (h_max - h_min) if (h_max - h_min) > 0 else 1.0

# Define template cube's faces layout
cube_faces = np.array([
    [0, 1, 2], [0, 2, 3],  # Bottom
    [4, 5, 6], [4, 6, 7],  # Top
    [0, 1, 5], [0, 5, 4],  # Front
    [1, 2, 6], [1, 6, 5],  # Right
    [2, 3, 7], [2, 7, 6],  # Back
    [3, 0, 4], [3, 4, 7]   # Left
])

bar_idx = 0
bar_width = 0.9  

for y in range(shape_y):
    for x in range(shape_x):
        h = h_matrix[y, x]
        v_offset = bar_idx * 8
        f_offset = bar_idx * 12
        
        # 8 unique corner vertices layout for a single 3D bar
        verts[v_offset:v_offset+8] = np.array([
            [x,             y,             0],
            [x + bar_width, y,             0],
            [x + bar_width, y + bar_width, 0],
            [x,             y + bar_width, 0],
            [x,             y,             h],
            [x + bar_width, y,             h],
            [x + bar_width, y + bar_width, h],
            [x,             y + bar_width, h]
        ])
        
        # Map faces to the assigned structural vertices index
        faces[f_offset:f_offset+12] = cube_faces + v_offset
        
        # Calculate colormap value based on individual bar height
        norm_h = (h - h_min) / h_range
        color = cm.map(norm_h) / 255.0  
        
        # Apply the color to all 12 triangular mesh faces of this bar
        face_colors[f_offset:f_offset+12] = color
        
        bar_idx += 1

# --- 5. Build and Add Mesh Item ---
# bg_mesh = gl.GLMeshItem(vertexes=verts, faces=faces, faceColors=face_colors, shader='shaded')
bg_mesh = gl.GLMeshItem(vertexes=verts, faces=faces, faceColors=face_colors, shader=None)

w.addItem(bg_mesh)

if __name__ == '__main__':
    # print(pg.colormap.listMaps())
    pg.exec()
