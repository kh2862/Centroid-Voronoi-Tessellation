#include "geometry.h"
#include "Voronoi_Tessellation.h"
#include "Centroid_Voronoi_Tessellation.h"

int main()
{
    Model model;
    model.input_generator("in.txt");
    model.Lloyd();
    model.output("out.txt");
    return 0;
}
