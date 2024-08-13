# Polygon Smoothing Algorithm
# by: Ulrike Holfeld & Nicolas Martinez Heredia


# Test data including their parameters
    # inFCName:        test  |  test2   |   test3  |   shape (large dataset)
    # offset_dist:     0.1   |  30000   |   0.05   |   0.01
    # reduction:       0.8   |  0.8     |   0.8    |   0.8
    # niterations:     3     |  3       |   3      |   2

    # The data varies in terms of how the polygon starts and ends (convex, concave) and tests how the algorithm closes the polygon
    # test.shp:     last vertex triple: convex,   first vertex triple: convex
    # test2.shp:    last vertex triple: concave,  first vertex triple: convex
    # test3.shp:    last vertex triple: concave,  first vertex triple: concave

import os
import math
import shutil
import time
import tkinter as tk
from tkinter import filedialog
from shapely.geometry import LineString, Polygon, Point
import shapefile
import matplotlib.pyplot as plt
    

# Source: N. Prechtel, 2024, shiftVec.py..................................................... # 
def createSubdir(workspace, subdirList):
    '''produce subdirs from a list of sudir names in a workspace'''
    for subdir in subdirList:
        if not os.path.isdir(workspace + "/" + subdir):
            os.mkdir(workspace + "/" + subdir)


# Source: N. Prechtel, 2024, shiftVec.py..................................................... # 
def controlExtension(inName, ext):
    '''enforce user-defined extension to the input file name'''
    if ext[0] != '.':
        ext = '.' + ext
    if inName.rfind(".") > 0:
        inName = inName[:inName.rfind(".")] + ext
    else:
        inName = inName + ext
    return inName


# Source: N. Prechtel, 2024, shiftVec.py..................................................... # 
def completePath(workspace, subdir, nameList):
    '''merge project dir, subdirectory and file names into complete path'''
    for ix in range(len(nameList)):
        nameList[ix] = workspace + "/" + subdir + "/" + nameList[ix]
    return nameList


# Source: N. Prechtel, 2024, shiftVec.py..................................................... #  
def checkExistence(pathList):
    '''do specified datasets exist under the entered paths ?'''
    check = True
    for data in pathList:
        if not os.path.isfile(data):
            check = False
            break
    return check


def iterativeSmoothing(inFC, offset_dist, niterations, reduction):
    '''Iteratively smoothing polygon outline using reduction of offset distance'''
    # Read Shapefile
    sf = shapefile.Reader(inFC[:inFC.rfind('.')])
    geoType = sf.shapeType
    if not (geoType == 5 or geoType == 15):     #5, 15 are polygon
        print('! input has incompatible geometry')
        return
    if len(sf.shapes()) != 1:  # Check if there is only one shape
        print('! input contains more than one polygon')
        return
    
    final_smoothed_polygon = None

    # Get the single polygon shape
    shape = sf.shapes()[0]
    
    fig, axs = plt.subplots(1, niterations, figsize=(5 * niterations, 5))

    if niterations == 1:
        axs = [axs]  # Ensure axs is iterable for single iteration

    for i, ax in enumerate(axs, start=1):
        offset_lines = []
        polygon = Polygon(shape.points) if final_smoothed_polygon is None else final_smoothed_polygon
        
        # Plot original polygon
        _plotPolygon(ax, polygon, 'blue', offset_lines)

        # Smooth the polygon outline
        smooth_polygon_points = []
        smooth_polygon_points, offset_lines = _smoothPolygon(polygon, offset_dist)

        # Plot smoothed polygon
        _plotPolygon(ax, Polygon(smooth_polygon_points), 'red', offset_lines)

        ax.set_title(f"Iteration {i}")

        # Reduce the offset distance for the next iteration
        offset_dist *= reduction

        # Store the final smoothed polygon after the last iteration
        final_smoothed_polygon = Polygon(smooth_polygon_points)

    plt.tight_layout()
    plt.show()
    
    return final_smoothed_polygon


def _plotPolygon(ax, polygon, color, offset_lines=None):
    '''Plot the polygon and its offset lines'''
    x, y = polygon.exterior.xy
    ax.plot(x, y, color=color)
    if offset_lines:
        for line in offset_lines:
            x, y = line.xy
            ax.plot(x, y, color='green')


def _smoothPolygon(polygon, offset_dist):
    '''Smooth the polygon outline using offset distance and return the smoothed points and offset lines'''
    smooth_polygon_points = [] # List to store the smoothed points
    offset_lines = [] # List to store the offset lines
    skip_next_iteration = False  # Initialize the flag
    prev_offset_line = []  # Keep track of offset lines from consecutive iterations   
    section_type = {} # Dictionary to store the type of each section (concave or not concave)
    len_coords = len(polygon.exterior.coords)  # Number of vertices in the polygon

    # Iterate through the vertices of the polygon
    for i in range(len(polygon.exterior.coords)-2):
        # Check if the flag is set to skip the current iteration
        if skip_next_iteration:
            # Add i and "not concave" to the section_type dictionary
            section_type[i] = "not concave"
            # Reset the flag
            skip_next_iteration = False  
            continue  # Skip this iteration

        # Get the three consecutive points
        p1, p2, p3 = _consecutivePoints(polygon, i, len_coords)
        #print(p1, p2, p3) 

        # Add original point p1 to smooth_polygon_points if the previous section was not concave
        if section_type.get(i-1) != "concave":
            smooth_polygon_points.append(p1)
            
        # Calculate the cross-product of the vectors formed by the points
        cross_product = _cross(p1,p2,p3)
        
        # Determine convex or concave based on the sign of the cross-product
        if cross_product < 0:
            # Convex section: Calculate offset point perpendicular to the edge (p1, p2)
            offset_point, offset_line = _convexOffsetPoint(p1,p2, offset_dist)
            smooth_polygon_points.append(offset_point)
            offset_lines.append(offset_line)

            prev_offset_line = [] # Reset the prev_offset_line

            section_type[i] = "not concave" 

        elif cross_product > 0:
            # Concave section: Calculate offset point along the bisector direction
            offset_point, offset_line = _concaveOffsetPoint(p1, p2, p3, offset_dist)   
            smooth_polygon_points.append(offset_point)
            offset_lines.append(offset_line)

            # If offset-line and previous offset-line intersect, replace the last 2 points with the intersection point
            if section_type.get(i-1) == "concave":
                _handleIntersection(prev_offset_line, offset_line, smooth_polygon_points)

            # Update prev_offset_line
            prev_offset_line = [offset_line]

            # Check if it's not the last iteration
            if i != len(polygon.exterior.coords) - 3:
                # Check if the next iteration results in a convex section
                next_cross_product = _cross(p2, p3, Point(polygon.exterior.coords[(i + 3) % len(polygon.exterior.coords)]))
                if next_cross_product < 0:
                    # Set the flag to skip the next iteration
                    skip_next_iteration = True
                    # Append p3
                    smooth_polygon_points.append(p3)
            
            section_type[i] = "concave"

        else:
            # Collinear points: Cross-product is zero
            prev_offset_line = [] # Reset the prev_offset_line
            section_type[i] = "not concave"

    # Handle the last triplet separately
    smooth_polygon_points, offset_lines = _closePolygon(polygon, smooth_polygon_points, offset_dist, offset_lines, section_type, prev_offset_line)

    # Print the dictionary
    #for key, value in section_type.items():
    #    print(f"Iteration {key}: {value}")

    return smooth_polygon_points, offset_lines


def _consecutivePoints(polygon, index, len_coords):
    '''Get the three consecutive points from the polygon's exterior coordinates'''
    # Using "%"" to wrap the index to the start of the list when it reaches the end
    p1 = Point(polygon.exterior.coords[index % len_coords])
    p2 = Point(polygon.exterior.coords[(index + 1) % len_coords])
    p3 = Point(polygon.exterior.coords[(index + 2) % len_coords])
    return p1, p2, p3


def _cross(p1, p2, p3):
    '''Calculate the cross-product of the vectors formed by the points'''
    v1 = (p2.x - p1.x, p2.y - p1.y)
    v2 = (p3.x - p2.x, p3.y - p2.y)
    cross = (v1[0] * v2[1]) - (v1[1] * v2[0])
    return cross


def _convexOffsetPoint(p1, p2, offset_dist):
    # Adapted code based on: N. Prechtel, 2024, shiftVec.py..................................................... # 
    
    '''Calculate the offset point perpendicular to the line (p1, p2) and return the offset point and offset line'''  
    # Calculate the midpoint (center) of the line segment
    cent = Point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2)
    
    # Calculate the angle of the line segment
    bear = math.atan2(p2.y - p1.y, p2.x - p1.x)

    # Calculate offset point
    dx = offset_dist * math.cos(bear + math.pi/2.)
    dy = offset_dist * math.sin(bear + math.pi/2.)
    xn = cent.x + dx
    yn = cent.y + dy
    offset_point = Point(xn, yn)

    # Create offset line
    offset_line = LineString([cent, offset_point])

    return offset_point, offset_line


def _concaveOffsetPoint(p1, p2, p3, offset_dist):
    '''Calculate the offset point along the bisector direction and return the offset point and offset line'''
    # Calculate the bisector angle
    angle1 = math.atan2(p2.y - p1.y, p2.x - p1.x)
    angle2 = math.atan2(p3.y - p2.y, p3.x - p2.x)
    bisector_angle = (angle1 + angle2) / 2
    bisector_angle -= math.pi / 2 # Rotate 90 degrees to point in the correct direction
   
    # Adjust the bisector angle if it's pointing inward
    if angle2 > angle1:
        bisector_angle += math.pi  # Add 180 degrees (pi radians)

    # Calculate the coordinates of the new endpoint
    xn = p2.x + offset_dist * math.cos(bisector_angle)
    yn = p2.y + offset_dist * math.sin(bisector_angle)
    
    offset_point = Point(xn, yn)

    # Create offset line
    offset_line = LineString([p2, offset_point])

    return offset_point, offset_line


def _handleIntersection(prev_offset_line, offset_line, smooth_polygon_points):
    '''Replace the last 2 points with the intersection point if the offset lines intersect'''
    # Check if prev_offset_line exists and intersects with current offset_line
    if prev_offset_line:
        for i in prev_offset_line:
            if offset_line.intersects(i):
                # Get the intersection point
                intersection_point = offset_line.intersection(i)
                # Replace the 2 offset points with 1 intersection point
                smooth_polygon_points[-2:] = [intersection_point]
                break


def _closePolygon(polygon, smooth_polygon_points, offset_dist, offset_lines, section_type, prev_offset_line):
    '''Handle the last triplet separately and close the polygon'''
    p1 = Point(polygon.exterior.coords[-2])  # Second-to-last point
    p2 = Point(polygon.exterior.coords[-1])  # Last point (is also first point of the polygon)
    p3 = Point(polygon.exterior.coords[1])   # Second point (not [0] because the last and first points are the same)

    #print(p1, p2, p3)

    # Add original point p1 to smooth_polygon_points if the highest key (= previous section) was not concave
    if section_type[max(section_type.keys())] != "concave":
        smooth_polygon_points.append(p1)
    
    # Calculate the cross-product of the vectors formed by the points
    cross_product = _cross(p1,p2,p3)
    
    # Determine convex/concave/collinear based on the sign of the cross-product
    if cross_product < 0:        
        # Convex section when previous section was not concave: Calculate offset point perpendicular to the edge (p1, p2)
        if section_type[max(section_type.keys())] != "concave":
            offset_point, offset_line = _convexOffsetPoint(p1,p2, offset_dist)
            smooth_polygon_points.append(offset_point)
            offset_lines.append(offset_line)

        # Add the last point of the polygon (to close the polygon)
        smooth_polygon_points.append(smooth_polygon_points[0])

        section_type[len(section_type)] = "not concave"

    elif cross_product > 0:
        # Concave section: Calculate offset point along the bisector direction
        offset_point, offset_line = _concaveOffsetPoint(p1, p2, p3, offset_dist)   
        smooth_polygon_points.append(offset_point)
        offset_lines.append(offset_line)

        # If offset-line and previous offset-line intersect, replace the last 2 points with the intersection point
        _handleIntersection(prev_offset_line, offset_line, smooth_polygon_points)

        # Handle different cases to close the polygon (depending on whether there was a previous intersection and what the first section type is):

        # There was an intersection with the last offset line
        if smooth_polygon_points[-1] != offset_point:   
            # First section is a concave section too
            if section_type[sorted(section_type.keys())[0]] == "concave":
                # Check for intersection 
                if offset_line.intersects(offset_lines[0]):
                    intersection_point = offset_line.intersection(offset_lines[0]) # Get the intersection point
                    smooth_polygon_points.pop(0)  # Remove the first element
                    smooth_polygon_points[-1] = intersection_point  # Replace the last element with the intersection point
                    smooth_polygon_points.pop(1)  # Remove the second point
                # No intersection
                else:
                    # Replace the first point with the last one
                    smooth_polygon_points[0] = smooth_polygon_points[-1]
            
            # First section is not concave (either convex or collinear)
            else:
                # Replace the first point with the last one
                smooth_polygon_points[0] = smooth_polygon_points[-1]
                # Only for convex sections
                if  _cross(first_cross_product) != 0:
                    # Delete the second point and the first offset line
                    smooth_polygon_points.pop(1)  # Remove the second point
                    offset_lines.pop(0)  # Remove the first offset line
        
        # There was no intersection with the last offset line
        else:
            smooth_polygon_points[0] = smooth_polygon_points[-1] # Replace the first point with the last one
            # Calculate the cross-product of the first section
            first_cross_product = _cross(Point(polygon.exterior.coords[(0)]),Point(polygon.exterior.coords[(1)]),Point(polygon.exterior.coords[(2)]))  
            # Check if the first section was convex
            if section_type[sorted(section_type.keys())[0]] != "concave" and first_cross_product != 0:
                smooth_polygon_points.pop(1)  # Remove the second point
                offset_lines.pop(0)  # Remove the first offset line
            
        section_type[len(section_type)] = "concave"
    

    else:
        # Collinear points
        smooth_polygon_points.append(smooth_polygon_points[0]) #Equals p2
        section_type[len(section_type)] = "not concave"
    
    return smooth_polygon_points, offset_lines


def copy_projection(inFC, outFC):
    '''Copy the projection file (.prj) from input to output'''
    shutil.copyfile(inFC[:inFC.rfind('.')] + '.prj',
                    outFC[:outFC.rfind('.')] + '.prj')


def writeShape(outFC, sf):
    '''Write the smoothed polygon to a shapefile'''
    w = shapefile.Writer(outFC[:outFC.rfind('.')], shapeType=shapefile.POLYGON)
    w.field('ID', 'N')  # Define a numeric attribute field
    w.poly([sf.exterior.coords]) # Write the polygon shape
    w.record(0)  # Assign ID 0
    w.close() # Close the shapefile writer
    
    print('Wrote shapefile:', outFC)


def runProgram():
    ''' Initiates the code with user input from the tkinter GUI'''
    worksp = workspace_entry.get()
    inFCName = inFCName_entry.get()
    outFCName = outFCName_entry.get()
    offset_dist = float(offset_dist_entry.get())
    reduction = float(reduction_entry.get())
    niterations = int(niterations_entry.get())

    ''' Code execution workflow'''
    stime = time.time()
    createSubdir(worksp, ['Shape', 'Output']) # Creates these 2 folders if they dont exist yet
    inFCName = controlExtension(inFCName, '.shp') # Output: shape.shp 
    inFC = completePath(worksp, 'Shape', [inFCName])[0]
    outFCName = controlExtension(outFCName, '.shp')
    outFC = completePath(worksp, 'Output', [outFCName])[0]

    if not checkExistence([inFC]):
        print('! input data incomplete or in wrong folder')
        exit(1)

    final_smoothed_polygon = iterativeSmoothing(inFC, offset_dist, niterations, reduction)
    copy_projection(inFC, outFC)
    writeShape(outFC, final_smoothed_polygon)
    extime = time.time()-stime
    print('Exectution time [s]: '+str(round(extime,3)))


def browseDirectory():
    ''' Opens a file dialog to browse for a directory and inserts the path into the workspace entry field'''
    directory = filedialog.askdirectory()
    workspace_entry.delete(0, tk.END)
    workspace_entry.insert(0, directory)


''' GUI '''

root = tk.Tk()

workspace_label = tk.Label(root, text="Workspace")
workspace_label.pack()
workspace_entry = tk.Entry(root)
workspace_entry.pack()
browse_button = tk.Button(root, text="Browse", command=browseDirectory)
browse_button.pack()

inFCName_label = tk.Label(root, text="Input Shapefile Name")
inFCName_label.pack()
inFCName_entry = tk.Entry(root)
inFCName_entry.pack()

outFCName_label = tk.Label(root, text="Output Shapefile Name")
outFCName_label.pack()
outFCName_entry = tk.Entry(root)
outFCName_entry.pack()

offset_dist_label = tk.Label(root, text="Offset Distance")
offset_dist_label.pack()
offset_dist_entry = tk.Entry(root)
offset_dist_entry.pack()

reduction_label = tk.Label(root, text="Reduction")
reduction_label.pack()
reduction_entry = tk.Entry(root)
reduction_entry.pack()

niterations_label = tk.Label(root, text="Number of Iterations")
niterations_label.pack()
niterations_entry = tk.Entry(root)
niterations_entry.pack()

run_button = tk.Button(root, text="Run Program", command=runProgram)
run_button.pack()

root.mainloop()