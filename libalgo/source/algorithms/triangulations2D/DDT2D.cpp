#include "DDT2D.h"

//Initialize pointer to swapping criteria (static)
double ( *DDT2D::pcriterion ) ( const Point3DCartesian <double> *, const Point3DCartesian <double> *, const Point3DCartesian <double> *, const Point3DCartesian <double> *, const Point3DCartesian <double> *, const Point3DCartesian <double> * )  =  &SwappingCriteria::getAbn <double>;


