#ifndef GEOJSONEXPORT_HPP
#define	GEOJSONEXPORT_HPP

template <typename Point>
void GeoJSONExport::exportGraticule ( std::ostream& out, const typename TMeridiansListF <typename Point::Type> ::Type & meridians,
                                       const typename TParallelsListF <typename Point::Type> ::Type & parallels, const Container <Point *> &points,
                                       const typename Point::Type font_height, const typename Point::Type lat_step, const typename Point::Type lon_step )
{
    //Create header section
    createHeader ( out );

    //Process all meridians
    typename TMeridiansListF <typename Point::Type> ::Type::const_iterator i_meridians = meridians.begin();

    for ( i_meridians = meridians.begin(); i_meridians != meridians.end(); ++i_meridians )
    {
        bool isLast = false;
        if (parallels.size() == 0 && i_meridians == --meridians.end()) {
            isLast = true;
        }
        processGraticuleElements ( out, *i_meridians, points, "Meridians", "Meridian_labels", font_height, lon_step, isLast );
    }

    //Process all parallels
    typename TParallelsListF <typename Point::Type> ::Type::const_iterator i_parallels = parallels.begin();

    for ( i_parallels = parallels.begin(); i_parallels != parallels.end(); ++i_parallels )
    {
        bool isLast = false;
        if (i_parallels == --parallels.end()) {
            isLast = true;
        }
        processGraticuleElements ( out, *i_parallels, points, "Parallels", "Parallel_labels", font_height, lat_step, isLast );
    }

    //End header section
    endHeader ( out );
}

template <typename GraticulePart, typename Point>
void GeoJSONExport::processGraticuleElements ( std::ostream& out, GraticulePart & part, const Container <Point*> & points, const char * layer_graticule_name, const char * layer_labels_name,
                const typename Point::Type font_height, const typename Point::Type step, bool isLast )
{
        //Export meridian or parallel (defined as a template parameter GraticulePart)
        TIndexList ind_points = part.getPointsIndices();

        //Process all points
        const unsigned int n = ind_points.size();

        for ( unsigned int i = 1; i < n; i++ )
        {
                //Get actual point
                Point * p = points [ind_points[i]];

                //Get previous point
                Point * p_previous = points [ind_points[i - 1]];

                //Create a line> part of the meridian or parallel
                createLine ( out, layer_graticule_name, p_previous->getX(), p_previous->getY(), p->getX(), p->getY() );

                out << ",";
        }

        //Create label
        char point_id_text [255];

        if ( n > 1 )
        {
                //Set accuracy depending on a step
                if ( step > 1.0 ) sprintf ( point_id_text, "%3.1f", part.getCoord () );
                else if ( step > 0.1 ) sprintf ( point_id_text, "%3.2f", part.getCoord () );
                else sprintf ( point_id_text, "%3.3f", part.getCoord () );

                //Compute bearing
                const typename Point::Type bearing = Bearing::getBearing ( points [ind_points[0.5 * n - 1 ]], points [ind_points[0.5 * n]] );

                //Create label for meridian/parallel
                createText ( out, layer_labels_name, point_id_text, points[ind_points[0.5 * n - 1 ]]->getX() + 0.5 * font_height * cos ( bearing * M_PI / 180 ), points [ind_points[0.5 * n - 1]]->getY() + 0.5 * font_height * sin ( bearing * M_PI / 180 ), bearing, font_height );
                if (!isLast) {
                    out << ",";
                }
        }
}

inline void GeoJSONExport::createHeader( std::ostream& out) {
    out << "{";
    out << "\"type\":\"FeatureCollection\",";
    out << "\"features\":[";
}

inline void GeoJSONExport::endHeader( std::ostream& out) {
    out << "]";
    out << "}";
}

template <typename T>
void GeoJSONExport::createLine ( std::ostream& out, const char * layer_name, const T x1, const T y1, const T x2, const T y2 )
{
    out << "{";
    out << "\"geometry\":{";
    out << "\"type\":\"LineString\",";
    out << "\"coordinates\":[";
    out << "[" << x1 << "," << y1 << "],";
    out << "[" << x2 << "," << y2 << "]";
    out << "]";
    out << "},";
    out << "\"type\":\"Feature\",";
    out << "\"properties\":{";
    out << "\"Layer\":\"" << layer_name << "\"";
    out << "}";
    out << "}";
}


template <typename T>
void GeoJSONExport::createPoint ( std::ostream& out, const char * layer_name, const T x, const T y )
{
    out << "{";
    out << "\"geometry\":{";
    out << "\"type\":\"Point\",";
    out << "\"coordinates\":[" << x << "," << y << "]";
    out << "},";
    out << "\"type\":\"Feature\",";
    out << "\"properties\":{";
    out << "\"Layer\":\"" << layer_name << "\"";
    out << "}";
    out << "}";
}


template <typename T>
void GeoJSONExport::createText ( std::ostream& out, const char * layer_name, const char * text, const T x, const T y, const T rotation, const T height )
{
    out << "{";
    out << "\"geometry\":{";
    out << "\"type\":\"Point\",";
    out << "\"coordinates\":[" << x << "," << y << "]";
    out << "},";
    out << "\"type\":\"Feature\",";
    out << "\"properties\":{";
    out << "\"Layer\":\"" << layer_name << "\",";
    out << "\"label\":\"" << text << "\"";
    out << "}";
    out << "}";
}

#endif	/* GEOJSONEXPORT_HPP */
