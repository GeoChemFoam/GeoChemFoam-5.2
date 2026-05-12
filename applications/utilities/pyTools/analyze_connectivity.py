import numpy as np
from scipy.ndimage import label
import argparse
import os
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

try:
    import cc3d
    CC3D_AVAILABLE = True
except ImportError:
    CC3D_AVAILABLE = False
    print("Warning: cc3d not available. Install with 'pip install connected-components-3d' for faster 3D analysis.")
    print("Falling back to scipy.ndimage.label (slower).")

def load_raw_image(filepath, shape):
    """
    Load a raw image file into a numpy array.
    
    Args:
        filepath (str): Path to the raw image file
        shape (tuple): Shape of the image (height, width) for 2D or (depth, height, width) for 3D
    
    Returns:
        numpy.ndarray: The loaded image data
    """
    return np.fromfile(filepath, dtype=np.uint8).reshape(shape)

def create_phase_mask(image, phase_ranges):
    """
    Create a binary mask for a specific phase based on value ranges.
    
    Args:
        image (numpy.ndarray): Input image data
        phase_ranges (list): List of (min, max) tuples defining the phase ranges
    
    Returns:
        numpy.ndarray: Binary mask where True indicates the specified phase
    """
    mask = np.zeros_like(image, dtype=bool)
    for min_val, max_val in phase_ranges:
        mask |= (image >= min_val) & (image <= max_val)
    return mask

def calculate_porosity(image, pore_ranges):
    """
    Calculate the porosity of the image based on pore ranges.
    
    Args:
        image (numpy.ndarray): Input image data
        pore_ranges (list): List of (min, max) tuples defining the pore ranges
    
    Returns:
        float: Porosity as a percentage
    """
    pore_mask = create_phase_mask(image, pore_ranges)
    return (np.sum(pore_mask) / pore_mask.size) * 100

def get_connected_pores_mask_2d(image, pore_ranges, direction):
    """
    Get a mask of connected pores in the specified direction for 2D data.
    
    Args:
        image (numpy.ndarray): Input 2D image data
        pore_ranges (list): List of (min, max) tuples defining the pore ranges
        direction (str): Direction to analyze ('x', 'y')
    
    Returns:
        numpy.ndarray: Binary mask where True indicates connected pores
    """
    # Create binary mask for pores
    pore_mask = create_phase_mask(image, pore_ranges)
    
    # Label all connected components
    labeled_array, num_features = label(pore_mask)
    
    if num_features == 0:
        return np.zeros_like(pore_mask, dtype=bool)
    
    # For 2D micromodels, we want to find pores that span across the image
    # in the specified direction
    
    if direction == 'x':
        # Check which components span from left to right
        spanning_components = []
        for i in range(1, num_features + 1):
            component_mask = (labeled_array == i)
            # Get column indices where this component exists
            cols = np.where(np.any(component_mask, axis=0))[0]
            if len(cols) > 0:
                # Check if component spans a significant portion of the width
                span_ratio = (cols.max() - cols.min() + 1) / image.shape[1]
                # Also check if it touches both edges
                touches_left = cols.min() == 0
                touches_right = cols.max() == image.shape[1] - 1
                
                if touches_left and touches_right:
                    spanning_components.append(i)
                elif span_ratio > 0.5:  # Spans more than 50% of width
                    spanning_components.append(i)
    
    elif direction == 'y':
        # Check which components span from top to bottom
        spanning_components = []
        for i in range(1, num_features + 1):
            component_mask = (labeled_array == i)
            # Get row indices where this component exists
            rows = np.where(np.any(component_mask, axis=1))[0]
            if len(rows) > 0:
                # Check if component spans a significant portion of the height
                span_ratio = (rows.max() - rows.min() + 1) / image.shape[0]
                # Also check if it touches both edges
                touches_top = rows.min() == 0
                touches_bottom = rows.max() == image.shape[0] - 1
                
                if touches_top and touches_bottom:
                    spanning_components.append(i)
                elif span_ratio > 0.5:  # Spans more than 50% of height
                    spanning_components.append(i)
    
    # If no spanning components found, use the largest component
    if not spanning_components:
        print(f"Warning: No components span across the image in {direction}-direction.")
        print("Using the largest connected component instead.")
        component_sizes = np.bincount(labeled_array.ravel())[1:]  # Skip background (0)
        largest_component = np.argmax(component_sizes) + 1
        spanning_components = [largest_component]
    
    # Create mask for spanning components
    connected_pores_mask = np.zeros_like(pore_mask, dtype=bool)
    for comp_id in spanning_components:
        connected_pores_mask |= (labeled_array == comp_id)
    
    return connected_pores_mask

def get_connected_pores_mask_3d(image, pore_ranges, direction):
    """
    Get a mask of connected pores in the specified direction for 3D data.
    Uses efficient 3D connectivity analysis with cc3d when available.
    
    Args:
        image (numpy.ndarray): Input 3D image data
        pore_ranges (list): List of (min, max) tuples defining the pore ranges
        direction (str): Direction to analyze ('x', 'y', 'z')
    
    Returns:
        numpy.ndarray: Binary mask where True indicates connected pores
    """
    # Create binary mask for pores
    pore_mask = create_phase_mask(image, pore_ranges)
    
    # Map direction to numpy axis
    direction_map = {'x': 2, 'y': 1, 'z': 0}
    axis = direction_map[direction]
    
    print(f"Performing 3D connected component analysis...")
    
    # Use cc3d if available (much faster), otherwise fall back to scipy
    if CC3D_AVAILABLE:
        # cc3d is much faster for 3D volumes - use connectivity=6 for face-connected only
        labeled_array = cc3d.connected_components(pore_mask.astype(np.uint8), connectivity=6)
        stats = cc3d.statistics(labeled_array)
        num_features = len(stats['voxel_counts']) - 1  # Exclude background
        
        if num_features == 0:
            return np.zeros_like(pore_mask, dtype=bool)
        
        print(f"Found {num_features:,} connected components using cc3d")
        
        # Get bounding boxes efficiently using cc3d
        # cc3d bounding boxes format: (slice(z1,z2), slice(y1,y2), slice(x1,x2))
        cc3d_bboxes = stats['bounding_boxes'][1:]  # Skip background
        voxel_counts = stats['voxel_counts'][1:]  # Skip background
        
        # Convert cc3d slice format to [z_min, y_min, x_min, z_max, y_max, x_max]
        bboxes = []
        for bbox_slices in cc3d_bboxes:
            if len(bbox_slices) == 3:
                z_slice, y_slice, x_slice = bbox_slices
                z_min, z_max = z_slice.start, z_slice.stop - 1
                y_min, y_max = y_slice.start, y_slice.stop - 1  
                x_min, x_max = x_slice.start, x_slice.stop - 1
                bboxes.append([z_min, y_min, x_min, z_max, y_max, x_max])
            else:
                # Handle unexpected format
                print(f"Warning: Unexpected bbox format: {bbox_slices}")
                bboxes.append([0, 0, 0, 1, 1, 1])
        
    else:
        # Fallback to scipy (slower) - use structure for 6-connectivity (face-connected only)
        # Create 3D structure for 6-connectivity (face neighbors only)
        structure = np.zeros((3, 3, 3), dtype=bool)
        structure[1, 1, 0] = True  # front
        structure[1, 1, 2] = True  # back
        structure[1, 0, 1] = True  # top
        structure[1, 2, 1] = True  # bottom
        structure[0, 1, 1] = True  # left
        structure[2, 1, 1] = True  # right
        labeled_array, num_features = label(pore_mask, structure=structure)
        
        if num_features == 0:
            return np.zeros_like(pore_mask, dtype=bool)
        
        print(f"Found {num_features:,} connected components using scipy (consider installing cc3d for speed)")
        
        # Calculate bounding boxes manually (slower)
        bboxes = []
        voxel_counts = []
        for i in range(1, num_features + 1):
            coords = np.where(labeled_array == i)
            if len(coords[0]) > 0:
                bbox = [coords[j].min() for j in range(3)] + [coords[j].max() for j in range(3)]
                bboxes.append(bbox)
                voxel_counts.append(len(coords[0]))
            else:
                bboxes.append([0, 0, 0, 0, 0, 0])
                voxel_counts.append(0)
    
    # Find spanning components efficiently using vectorized operations
    spanning_components = []
    
    # Vectorized spanning detection
    for i, (bbox, voxel_count) in enumerate(zip(bboxes, voxel_counts), 1):
        if voxel_count == 0:
            continue
            
        z_min, y_min, x_min, z_max, y_max, x_max = bbox
        
        # Check spanning in the specified direction
        if axis == 0:  # z-direction
            touches_front = z_min == 0
            touches_back = z_max == image.shape[0] - 1
            span_ratio = (z_max - z_min + 1) / image.shape[0]
        elif axis == 1:  # y-direction
            touches_front = y_min == 0
            touches_back = y_max == image.shape[1] - 1
            span_ratio = (y_max - y_min + 1) / image.shape[1]
        else:  # x-direction
            touches_front = x_min == 0
            touches_back = x_max == image.shape[2] - 1
            span_ratio = (x_max - x_min + 1) / image.shape[2]
        
        # Component spans if it touches both boundaries or spans >70% of dimension
        if touches_front and touches_back:
            spanning_components.append(i)
        elif span_ratio > 0.7:
            spanning_components.append(i)
    
    # If no spanning components, use percolation analysis on largest components only
    if not spanning_components:
        print(f"No components span across the volume in {direction}-direction.")
        print("Analyzing percolation using largest components...")
        
        # Only analyze the largest components (top 10% by size or max 100 components)
        total_pores = np.sum(pore_mask)
        size_threshold = max(total_pores * 0.001, 100)  # At least 0.1% of total pores or 100 voxels
        
        # Get indices of large components, sorted by size
        large_component_indices = []
        for i, voxel_count in enumerate(voxel_counts, 1):
            if voxel_count >= size_threshold:
                large_component_indices.append((i, voxel_count))
        
        # Sort by size and take top 100
        large_component_indices.sort(key=lambda x: x[1], reverse=True)
        large_component_indices = large_component_indices[:100]
        
        print(f"Analyzing {len(large_component_indices)} largest components for percolation...")
        
        # Vectorized boundary touching analysis for large components
        for comp_id, _ in large_component_indices:
            bbox = bboxes[comp_id - 1]
            z_min, y_min, x_min, z_max, y_max, x_max = bbox
            
            # Count boundary touches using bbox
            boundaries_touched = 0
            if z_min == 0: boundaries_touched += 1
            if z_max == image.shape[0] - 1: boundaries_touched += 1
            if y_min == 0: boundaries_touched += 1
            if y_max == image.shape[1] - 1: boundaries_touched += 1
            if x_min == 0: boundaries_touched += 1
            if x_max == image.shape[2] - 1: boundaries_touched += 1
            
            # Include if touches multiple boundaries
            if boundaries_touched >= 3:
                spanning_components.append(comp_id)
        
        # If still no components, use the largest
        if not spanning_components:
            print("Using the largest connected component as fallback.")
            largest_idx = np.argmax(voxel_counts)
            spanning_components = [largest_idx + 1]
        else:
            print(f"Found {len(spanning_components)} percolating components.")
    else:
        print(f"Found {len(spanning_components)} spanning components.")
    
    # Create connected pores mask efficiently
    connected_pores_mask = np.zeros_like(pore_mask, dtype=bool)
    
    if CC3D_AVAILABLE and len(spanning_components) > 1:
        # Use cc3d for efficient multi-component extraction
        component_mask = np.isin(labeled_array, spanning_components)
        connected_pores_mask = component_mask & pore_mask
    else:
        # Standard approach for single component or when cc3d not available
        for comp_id in spanning_components:
            connected_pores_mask |= (labeled_array == comp_id)
    
    return connected_pores_mask

def get_connected_pores_mask(image, pore_ranges, direction):
    """
    Get a mask of connected pores in the specified direction.
    Automatically detects 2D vs 3D data and uses appropriate method.
    
    Args:
        image (numpy.ndarray): Input image data
        pore_ranges (list): List of (min, max) tuples defining the pore ranges
        direction (str): Direction to analyze ('x', 'y', 'z')
    
    Returns:
        numpy.ndarray: Binary mask where True indicates connected pores
    """
    if len(image.shape) == 2:
        # 2D data
        if direction == 'z':
            print("Warning: z-direction specified for 2D data. Using x-direction instead.")
            direction = 'x'
        return get_connected_pores_mask_2d(image, pore_ranges, direction)
    elif len(image.shape) == 3:
        # 3D data
        return get_connected_pores_mask_3d(image, pore_ranges, direction)
    else:
        raise ValueError(f"Unsupported image dimensions: {image.shape}")

def create_modified_image(image, pore_ranges, solid_ranges, direction):
    """
    Create a new image with only connected pores, replacing unconnected pores with median solid value.
    
    Args:
        image (numpy.ndarray): Input image data
        pore_ranges (list): List of (min, max) tuples defining the pore ranges
        solid_ranges (list): List of (min, max) tuples defining the solid ranges
        direction (str): Direction to analyze ('x', 'y', 'z')
    
    Returns:
        tuple: (modified_image, unconnected_pores_mask)
    """
    # Create a copy of the original image
    modified_image = image.copy()
    
    # Get connected pores mask
    connected_pores_mask = get_connected_pores_mask(image, pore_ranges, direction)
    
    # Get solid mask
    solid_mask = create_phase_mask(image, solid_ranges)
    
    # Calculate median solid value
    median_solid = np.median(image[solid_mask])
    
    # Replace unconnected pores with median solid value
    pore_mask = create_phase_mask(image, pore_ranges)
    unconnected_pores = pore_mask & ~connected_pores_mask
    modified_image[unconnected_pores] = int(median_solid)
    
    return modified_image, unconnected_pores

def save_raw_image(image, filepath):
    """
    Save a numpy array as a raw image file.
    
    Args:
        image (numpy.ndarray): Image data to save
        filepath (str): Path to save the raw image file
    """
    image.astype(np.uint8).tofile(filepath)

def visualize_slice_2d(image, modified_image, unconnected_pores, pore_ranges, output_path):
    """
    Create a visualization for 2D data.
    """
    # Create a custom colormap: black for background, white for pores, red for unconnected pores
    colors = ['black', 'white', 'red']
    cmap = ListedColormap(colors)
    
    # Create a visualization array
    vis_array = np.zeros_like(image, dtype=np.uint8)
    pore_mask = create_phase_mask(image, pore_ranges)  # Use actual pore ranges
    vis_array[pore_mask] = 1  # Pores in white
    vis_array[unconnected_pores] = 2  # Unconnected pores in red
    
    num_unconnected = np.sum(unconnected_pores)
    
    # Create the figure - add more subplots if there are unconnected pores to show
    if num_unconnected > 0:
        fig = plt.figure(figsize=(18, 6))
        gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1])
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
        ax3 = fig.add_subplot(gs[2])
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
    
    # Plot original with unconnected pores highlighted
    ax1.imshow(vis_array, cmap=cmap, interpolation='nearest')
    ax1.set_title(f'Original Micromodel\n(Red = {num_unconnected} Unconnected Pores)')
    ax1.axis('off')
    
    # Plot modified
    ax2.imshow(modified_image, cmap='gray')
    ax2.set_title('Modified Micromodel\n(Unconnected Pores Removed)')
    ax2.axis('off')
    
    # If there are unconnected pores, create a zoomed view
    if num_unconnected > 0:
        # Find bounding box of unconnected pores
        rows, cols = np.where(unconnected_pores)
        min_row, max_row = rows.min(), rows.max()
        min_col, max_col = cols.min(), cols.max()
        
        # Add padding to show context
        padding = 20
        min_row = max(0, min_row - padding)
        max_row = min(image.shape[0], max_row + padding)
        min_col = max(0, min_col - padding)
        max_col = min(image.shape[1], max_col + padding)
        
        # Create zoomed view
        zoomed_vis = vis_array[min_row:max_row, min_col:max_col]
        ax3.imshow(zoomed_vis, cmap=cmap, interpolation='nearest')
        ax3.set_title(f'Zoomed View\n(Region: {min_row}:{max_row}, {min_col}:{max_col})')
        ax3.axis('off')
        
        # Add text showing exact pixel count
        fig.suptitle(f'Connectivity Analysis Results\nTotal pores: {np.sum(pore_mask):,}, Unconnected: {num_unconnected}, Connected: {np.sum(pore_mask) - num_unconnected:,}', 
                    fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

def visualize_slice_3d(image, modified_image, unconnected_pores, pore_ranges, slice_idx, direction, output_path):
    """
    Create a comprehensive 3D visualization showing multiple perspectives and projections.
    """
    # Create a custom colormap: black for background, white for pores, red for unconnected pores
    colors = ['black', 'white', 'red']
    cmap = ListedColormap(colors)
    
    # Create visualization array
    vis_array = np.zeros_like(image, dtype=np.uint8)
    pore_mask = create_phase_mask(image, pore_ranges)
    vis_array[pore_mask] = 1  # Pores in white
    vis_array[unconnected_pores] = 2  # Unconnected pores in red
    
    # Calculate statistics
    num_unconnected = np.sum(unconnected_pores)
    total_pores = np.sum(pore_mask)
    
    # Create comprehensive figure
    fig = plt.figure(figsize=(20, 16))
    
    # Create a grid layout for multiple subplots
    gs = fig.add_gridspec(4, 5, height_ratios=[1, 1, 1, 0.8], width_ratios=[1, 1, 1, 1, 1])
    
    # === ROW 1: Maximum Intensity Projections (MIP) ===
    
    # MIP along Z-axis (view from top)
    ax1 = fig.add_subplot(gs[0, 0])
    mip_z = np.max(vis_array, axis=0)
    ax1.imshow(mip_z, cmap=cmap, interpolation='nearest')
    ax1.set_title('MIP: Z-axis\n(Top View)')
    ax1.axis('off')
    
    # MIP along Y-axis (view from front)
    ax2 = fig.add_subplot(gs[0, 1])
    mip_y = np.max(vis_array, axis=1)
    ax2.imshow(mip_y, cmap=cmap, interpolation='nearest')
    ax2.set_title('MIP: Y-axis\n(Front View)')
    ax2.axis('off')
    
    # MIP along X-axis (view from side)
    ax3 = fig.add_subplot(gs[0, 2])
    mip_x = np.max(vis_array, axis=2)
    ax3.imshow(mip_x, cmap=cmap, interpolation='nearest')
    ax3.set_title('MIP: X-axis\n(Side View)')
    ax3.axis('off')
    
    # Connectivity statistics per slice
    ax4 = fig.add_subplot(gs[0, 3])
    direction_map = {'x': 2, 'y': 1, 'z': 0}
    primary_axis = direction_map[direction]
    
    # Calculate connected pores per slice along primary direction
    slice_connectivity = []
    for i in range(image.shape[primary_axis]):
        if primary_axis == 0:  # z-slices
            slice_pores = np.sum(pore_mask[i, :, :])
            slice_unconnected = np.sum(unconnected_pores[i, :, :])
        elif primary_axis == 1:  # y-slices
            slice_pores = np.sum(pore_mask[:, i, :])
            slice_unconnected = np.sum(unconnected_pores[:, i, :])
        else:  # x-slices
            slice_pores = np.sum(pore_mask[:, :, i])
            slice_unconnected = np.sum(unconnected_pores[:, :, i])
        
        connectivity_ratio = (slice_pores - slice_unconnected) / slice_pores if slice_pores > 0 else 0
        slice_connectivity.append(connectivity_ratio)
    
    ax4.plot(slice_connectivity, 'b-', linewidth=2)
    ax4.set_title(f'Connectivity Ratio\nper {direction.upper()}-slice')
    ax4.set_xlabel(f'{direction.upper()}-slice index')
    ax4.set_ylabel('Connected Ratio')
    ax4.grid(True, alpha=0.3)
    ax4.set_ylim(0, 1)
    
    # 3D scatter plot of unconnected pores
    ax5 = fig.add_subplot(gs[0, 4], projection='3d')
    if num_unconnected > 0:
        unc_coords = np.where(unconnected_pores)
        # Subsample if too many points (for performance)
        if len(unc_coords[0]) > 1000:
            indices = np.random.choice(len(unc_coords[0]), 1000, replace=False)
            z_coords = unc_coords[0][indices]
            y_coords = unc_coords[1][indices]
            x_coords = unc_coords[2][indices]
        else:
            z_coords, y_coords, x_coords = unc_coords
        
        ax5.scatter(x_coords, y_coords, z_coords, c='red', s=1, alpha=0.6)
        ax5.set_title(f'3D Distribution\nof Unconnected Pores')
    else:
        ax5.text(0.5, 0.5, 0.5, 'No Unconnected\nPores Found', 
                transform=ax5.transAxes, ha='center', va='center')
        ax5.set_title('3D Distribution\nof Unconnected Pores')
    
    ax5.set_xlabel('X')
    ax5.set_ylabel('Y')
    ax5.set_zlabel('Z')
    
    # === ROW 2: Central slices in each direction ===
    
    # Central Z-slice
    ax6 = fig.add_subplot(gs[1, 0])
    z_center = image.shape[0] // 2
    ax6.imshow(vis_array[z_center, :, :], cmap=cmap, interpolation='nearest')
    ax6.set_title(f'Central Z-slice\n(slice {z_center})')
    ax6.axis('off')
    
    # Central Y-slice
    ax7 = fig.add_subplot(gs[1, 1])
    y_center = image.shape[1] // 2
    ax7.imshow(vis_array[:, y_center, :], cmap=cmap, interpolation='nearest')
    ax7.set_title(f'Central Y-slice\n(slice {y_center})')
    ax7.axis('off')
    
    # Central X-slice
    ax8 = fig.add_subplot(gs[1, 2])
    x_center = image.shape[2] // 2
    ax8.imshow(vis_array[:, :, x_center], cmap=cmap, interpolation='nearest')
    ax8.set_title(f'Central X-slice\n(slice {x_center})')
    ax8.axis('off')
    
    # Porosity distribution along primary direction
    ax9 = fig.add_subplot(gs[1, 3])
    porosity_per_slice = []
    for i in range(image.shape[primary_axis]):
        if primary_axis == 0:  # z-slices
            total_pixels = image.shape[1] * image.shape[2]
            pore_pixels = np.sum(pore_mask[i, :, :])
        elif primary_axis == 1:  # y-slices
            total_pixels = image.shape[0] * image.shape[2]
            pore_pixels = np.sum(pore_mask[:, i, :])
        else:  # x-slices
            total_pixels = image.shape[0] * image.shape[1]
            pore_pixels = np.sum(pore_mask[:, :, i])
        
        porosity = (pore_pixels / total_pixels) * 100
        porosity_per_slice.append(porosity)
    
    ax9.plot(porosity_per_slice, 'g-', linewidth=2)
    ax9.set_title(f'Porosity per\n{direction.upper()}-slice (%)')
    ax9.set_xlabel(f'{direction.upper()}-slice index')
    ax9.set_ylabel('Porosity (%)')
    ax9.grid(True, alpha=0.3)
    
    # Volume rendering preview (simplified)
    ax10 = fig.add_subplot(gs[1, 4], projection='3d')
    # Create a simplified 3D representation by showing connected components
    # Use 6-connectivity structure for consistency with main analysis
    structure = np.zeros((3, 3, 3), dtype=bool)
    structure[1, 1, 0] = True  # front
    structure[1, 1, 2] = True  # back
    structure[1, 0, 1] = True  # top
    structure[1, 2, 1] = True  # bottom
    structure[0, 1, 1] = True  # left
    structure[2, 1, 1] = True  # right
    labeled_pores, num_components = label(pore_mask, structure=structure)
    if num_components > 0:
        # Show largest connected component
        component_sizes = np.bincount(labeled_pores.ravel())[1:]
        largest_comp_id = np.argmax(component_sizes) + 1
        largest_comp_coords = np.where(labeled_pores == largest_comp_id)
        
        # Subsample for visualization
        if len(largest_comp_coords[0]) > 2000:
            indices = np.random.choice(len(largest_comp_coords[0]), 2000, replace=False)
            z_coords = largest_comp_coords[0][indices]
            y_coords = largest_comp_coords[1][indices]
            x_coords = largest_comp_coords[2][indices]
        else:
            z_coords, y_coords, x_coords = largest_comp_coords
        
        ax10.scatter(x_coords, y_coords, z_coords, c='blue', s=0.5, alpha=0.1)
        ax10.set_title('3D Connected\nPore Network')
    else:
        ax10.text(0.5, 0.5, 0.5, 'No Connected\nPores Found', 
                 transform=ax10.transAxes, ha='center', va='center')
        ax10.set_title('3D Connected\nPore Network')
    
    ax10.set_xlabel('X')
    ax10.set_ylabel('Y')
    ax10.set_zlabel('Z')
    
    # === ROW 3: Modified image views ===
    
    # Modified central slices
    ax11 = fig.add_subplot(gs[2, 0])
    ax11.imshow(modified_image[z_center, :, :], cmap='gray')
    ax11.set_title(f'Modified Z-slice\n(slice {z_center})')
    ax11.axis('off')
    
    ax12 = fig.add_subplot(gs[2, 1])
    ax12.imshow(modified_image[:, y_center, :], cmap='gray')
    ax12.set_title(f'Modified Y-slice\n(slice {y_center})')
    ax12.axis('off')
    
    ax13 = fig.add_subplot(gs[2, 2])
    ax13.imshow(modified_image[:, :, x_center], cmap='gray')
    ax13.set_title(f'Modified X-slice\n(slice {x_center})')
    ax13.axis('off')
    
    # Before/after comparison histogram
    ax14 = fig.add_subplot(gs[2, 3])
    ax14.hist(image.ravel(), bins=50, alpha=0.7, label='Original', color='blue')
    ax14.hist(modified_image.ravel(), bins=50, alpha=0.7, label='Modified', color='red')
    ax14.set_title('Intensity Distribution')
    ax14.set_xlabel('Pixel Value')
    ax14.set_ylabel('Frequency')
    ax14.legend()
    ax14.grid(True, alpha=0.3)
    
    # Connected components size distribution
    ax15 = fig.add_subplot(gs[2, 4])
    if num_components > 0:
        component_sizes = np.bincount(labeled_pores.ravel())[1:]
        ax15.semilogy(sorted(component_sizes, reverse=True), 'o-', markersize=4)
        ax15.set_title('Connected Component\nSizes (log scale)')
        ax15.set_xlabel('Component Rank')
        ax15.set_ylabel('Component Size')
        ax15.grid(True, alpha=0.3)
    else:
        ax15.text(0.5, 0.5, 'No Connected\nComponents Found', 
                 transform=ax15.transAxes, ha='center', va='center')
        ax15.set_title('Connected Component\nSizes')
    
    # === ROW 4: Summary statistics ===
    
    # Create a text summary
    ax16 = fig.add_subplot(gs[3, :])
    ax16.axis('off')
    
    summary_text = f"""
    3D CONNECTIVITY ANALYSIS SUMMARY
    
    Volume Shape: {image.shape} (Z × Y × X)    |    Analysis Direction: {direction.upper()}
    
    POROSITY STATISTICS:
    • Total Volume: {np.prod(image.shape):,} voxels
    • Total Pores: {total_pores:,} voxels ({total_pores/np.prod(image.shape)*100:.2f}%)
    • Connected Pores: {total_pores - num_unconnected:,} voxels ({(total_pores - num_unconnected)/np.prod(image.shape)*100:.2f}%)
    • Unconnected Pores: {num_unconnected:,} voxels ({num_unconnected/np.prod(image.shape)*100:.6f}%)
    
    CONNECTIVITY STATISTICS:
    • Total Connected Components: {num_components}
    • Largest Component: {max(component_sizes) if num_components > 0 else 0:,} voxels
    • Connectivity Efficiency: {(total_pores - num_unconnected)/total_pores*100:.4f}%
    • Average Slice Connectivity: {np.mean(slice_connectivity)*100:.2f}%
    """
    
    ax16.text(0.05, 0.95, summary_text, transform=ax16.transAxes, 
             fontsize=10, verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))
    
    # Add color legend
    legend_elements = [
        plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='black', markersize=10, label='Background'),
        plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='white', markersize=10, label='Connected Pores'),
        plt.Line2D([0], [0], marker='s', color='w', markerfacecolor='red', markersize=10, label='Unconnected Pores')
    ]
    ax16.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(0.98, 0.95))
    
    plt.suptitle(f'3D Pore Connectivity Analysis: {direction.upper()}-direction', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

def visualize_slice(image, modified_image, unconnected_pores, pore_ranges, slice_idx, direction, output_path):
    """
    Create a visualization before and after processing.
    Automatically detects 2D vs 3D and uses appropriate method.
    """
    if len(image.shape) == 2:
        visualize_slice_2d(image, modified_image, unconnected_pores, pore_ranges, output_path)
    else:
        visualize_slice_3d(image, modified_image, unconnected_pores, pore_ranges, slice_idx, direction, output_path)

def validate_phase_coverage(image, phase_ranges_list):
    """
    Validate that all pixels in the image are covered by at least one phase.
    
    Args:
        image (numpy.ndarray): Input image data
        phase_ranges_list (list): List of phase ranges for each phase
    
    Raises:
        ValueError: If there are unclassified pixels
    """
    # Create a mask for all classified pixels
    classified_mask = np.zeros_like(image, dtype=bool)
    for phase_ranges in phase_ranges_list:
        classified_mask |= create_phase_mask(image, phase_ranges)
    
    # Check for unclassified pixels
    unclassified = ~classified_mask
    if np.any(unclassified):
        unclassified_values = np.unique(image[unclassified])
        raise ValueError(f"Found unclassified pixels with values: {unclassified_values}")

def parse_phase_ranges(phase_str):
    """
    Parse phase ranges from string format 'min1-max1,min2-max2,...'
    """
    ranges = []
    for range_str in phase_str.split(','):
        min_val, max_val = map(int, range_str.split('-'))
        ranges.append((min_val, max_val))
    return ranges

def main():
    parser = argparse.ArgumentParser(description='Analyze connectivity of phases in a raw image.')
    parser.add_argument('filepath', help='Path to the raw image file')
    parser.add_argument('--shape', nargs='+', type=int, required=True,
                      help='Shape of the image (e.g., "256 256" for 2D or "64 256 256" for 3D)')
    parser.add_argument('--pores', type=str, required=True,
                      help='Pore phase ranges (e.g., "0-50,200-255")')
    parser.add_argument('--solids', type=str, required=True,
                      help='Solid phase ranges (e.g., "51-199")')
    parser.add_argument('--micropores', type=str,
                      help='Micropore phase ranges (e.g., "100-150")')
    parser.add_argument('--direction', type=str, choices=['x', 'y', 'z'], required=True,
                      help='Direction to analyze connectivity')
    parser.add_argument('--output', type=str, required=True,
                      help='Path to save the modified image')
    parser.add_argument('--slice', type=int,
                      help='Index of the slice to visualize (optional)')
    parser.add_argument('--visualize', action='store_true',
                      help='Create visualization of before/after slices')
    
    args = parser.parse_args()
    
    # Convert shape argument to tuple
    shape = tuple(args.shape)
    
    try:
        # Load the image
        image = load_raw_image(args.filepath, shape)
        
        # Parse phase ranges
        pore_ranges = parse_phase_ranges(args.pores)
        solid_ranges = parse_phase_ranges(args.solids)
        
        # Collect all phase ranges for validation
        phase_ranges_list = [pore_ranges, solid_ranges]
        
        # Handle optional micropores
        if args.micropores:
            micropore_ranges = parse_phase_ranges(args.micropores)
            phase_ranges_list.append(micropore_ranges)
        
        # Validate that all pixels are classified
        validate_phase_coverage(image, phase_ranges_list)
        
        # Calculate initial porosity
        initial_porosity = calculate_porosity(image, pore_ranges)
        
        # Create modified image with only connected pores
        modified_image, unconnected_pores = create_modified_image(image, pore_ranges, solid_ranges, args.direction)
        
        # Calculate final porosity
        final_porosity = calculate_porosity(modified_image, pore_ranges)
        
        # Save the modified image
        save_raw_image(modified_image, args.output)
        
        # Create visualization if requested
        if args.visualize:
            vis_output = os.path.splitext(args.output)[0] + '_visualization.png'
            visualize_slice(image, modified_image, unconnected_pores, pore_ranges, args.slice, args.direction, vis_output)
            print(f"\nVisualization saved to: {vis_output}")
        
        # Print results
        print(f"\nImage Analysis:")
        print(f"Image dimensions: {len(shape)}D")
        print(f"Image shape: {shape}")
        
        # Print unique values and their distribution
        unique_values, counts = np.unique(image, return_counts=True)
        print(f"Unique pixel values: {unique_values}")
        print(f"Pixel value distribution:")
        for val, count in zip(unique_values, counts):
            percentage = (count / image.size) * 100
            print(f"  Value {val}: {count:,} pixels ({percentage:.2f}%)")
        
        # Print connectivity information
        pore_mask = create_phase_mask(image, pore_ranges)
        if len(shape) == 3:
            # Use 6-connectivity for 3D analysis consistency
            structure = np.zeros((3, 3, 3), dtype=bool)
            structure[1, 1, 0] = True  # front
            structure[1, 1, 2] = True  # back
            structure[1, 0, 1] = True  # top
            structure[1, 2, 1] = True  # bottom
            structure[0, 1, 1] = True  # left
            structure[2, 1, 1] = True  # right
            labeled_pores, num_components = label(pore_mask, structure=structure)
        else:
            # 2D uses default 8-connectivity (4-connectivity is too restrictive for 2D micromodels)
            labeled_pores, num_components = label(pore_mask)
        if num_components > 0:
            component_sizes = np.bincount(labeled_pores.ravel())[1:]  # Skip background
            print(f"Total connected pore components: {num_components}")
            print(f"Largest component size: {max(component_sizes)} pixels")
            print(f"Smallest component size: {min(component_sizes)} pixels")
        
        print(f"\nPorosity Analysis:")
        print(f"Initial porosity: {initial_porosity:.12f}%")
        print(f"Final porosity (after removing unconnected pores): {final_porosity:.12f}%")
        print(f"Porosity reduction: {initial_porosity - final_porosity:.12f}%")
        print(f"\nModified image saved to: {args.output}")
        
        if len(shape) == 2:
            print(f"2D connectivity analysis was performed in the {args.direction}-direction")
            print("Note: For 2D micromodels, 'x' means left-right connectivity, 'y' means top-bottom connectivity")
        else:
            print(f"3D connectivity analysis was performed in the {args.direction}-direction")
            
    except Exception as e:
        print(f"Error: {str(e)}")

if __name__ == "__main__":
    main() 