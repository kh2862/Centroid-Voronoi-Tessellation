# Centroid-Voronoi-Tessellation
Computing centroid Voronoi tessellation by Fortune's algorithm

Naive version

+ Assume points pairwise different, should be modify
+ Assume points do not lie on boundary, should be modify
+ how to determine converge? the code falls on this issue
+ circle events happens one after another, do not know the treatment is okay or not(eg. insert site degenerate case)
+ correct: geomerty, fill boundary, event class
+ important: the final doubly connected edge is not correct, some 0-length edges exist. but make_polygon process eliminate this.
+ todo: overload point==!=; comments; bounding box,precision put where?; event destroy; return 1.0; this pointer, do_voro_te more initialize
+ implement rope and use it in the algorithm(by splay tree)
