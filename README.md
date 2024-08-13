# Project: Implementation of efficient polygon smoothing without area shrinkage
Class: Subject-Specific GIS Applications, Technical University of Dresden\
Group: Ulrike Holfeld, Nicolas Martinez Heredia\
Date: 01.03.2024

# Objective
Implementation of an algorithm which efficiently smooths a polygon around an arbitrary cell/point distribution without ‘loosing’ cells after application. 

# Approach
The polygon is assumed to be uniformly digitised clockwise. New “smoothing” vertices are added as follows: In a convex section, one new point which is offset from the centre of the segment in normal direction by a distance # d is entered. In concave sections, the original vertex is replaced by a new one offset to the left in bisector direction. If two such bisector offsets (all shown in red in the figure) will intersect within distance d, their # intersection point is used instead of two offset vertices. The smoothing can be done a couple of times while steadily reducing the distances d iteration by iteration. 

# Test data
To use the test data (in the folder "Shape"), please specify the variables in the GUI. The parameters that can be used for each dataset are added to the code as a comment.

