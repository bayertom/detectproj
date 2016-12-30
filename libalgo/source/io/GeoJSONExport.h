#ifndef GEOJSONEXPORT_H
#define	GEOJSONEXPORT_H

#include <iostream>
#include "libalgo/source/algorithms/graticule/Graticule.h"

class GeoJSONExport {
    
public:
    template <typename Point>
    static void exportGraticule ( std::ostream& out, const typename TMeridiansListF <typename Point::Type> ::Type & meridians, const typename TParallelsListF <typename Point::Type> ::Type & parallels,
                                       const Container <Point *> &points, const typename Point::Type font_height, const typename Point::Type lat_step, const typename Point::Type lon_step );
    
private:
    template <typename GraticulePart, typename Point>
    static void processGraticuleElements ( std::ostream& out, GraticulePart & part, const Container <Point*> & points, const char * layer_name_edges, const char * layer_name_generators,
                                           const typename Point::Type font_height, const typename Point::Type step, bool isLast );
    
    static void createHeader( std::ostream& out);
    
    static void endHeader( std::ostream& out);
    
    template <typename T>
    static void createLine ( std::ostream& out, const char * layer_name, const T x1, const T y1, const T x2, const T y2 );

    template <typename T>
    static void createPoint ( std::ostream& out, const char * layer_name, const T x, const T y );

    template <typename T>
    static void createText ( std::ostream& out, const char * layer_name, const char * text, const T x, const T y, const T rotation, const T height );
    
};

#include "GeoJSONExport.hpp"

#endif	/* GEOJSONEXPORT_H */

