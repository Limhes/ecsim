#include "simulation.h"
#include "electrodes.h"

extern const double CONST_PI; // in simulation.h
extern const string NO_GEOM; // in electrodes.h

void Electrode::update()
{
    if (m_type == 0) // disk
    {
        electrodeGeometryFactor = 0.0;
        electrodeArea = CONST_PI*m_geom1*m_geom1;
    }
    else if (m_type == 1) // square
    {
        electrodeGeometryFactor = 0.0;
        electrodeArea = m_geom1*m_geom1;
    }
    else if (m_type == 2) // rectangle
    {
        electrodeGeometryFactor = 0.0;
        electrodeArea = m_geom1*m_geom2;
    }
    else if (m_type == 3) // cylinder
    {
        electrodeGeometryFactor = 1.0;
        electrodeArea = 2*CONST_PI*m_geom1*m_geom2;
    }
    else if (m_type == 4) // sphere
    {
        electrodeGeometryFactor = 2.0;
        electrodeArea = 4*CONST_PI*m_geom1*m_geom1;
    }
    else if (m_type == 5) // hemisphere
    {
        electrodeGeometryFactor = 2.0;
        electrodeArea = 2*CONST_PI*m_geom1*m_geom1;
    }
    epsilon = sqrt(electrodeArea);
}

const vector<string> Electrode::electrodeTypes({"Disk", "Square", "Rectangle", "Cylinder", "Sphere", "Hemisphere"});
const vector<string> Electrode::electrodeGeom1({"Radius", "Width", "Width", "Radius", "Radius", "Radius"});
const vector<string> Electrode::electrodeGeom2({NO_GEOM, NO_GEOM, "Length", "Length", NO_GEOM, NO_GEOM});
