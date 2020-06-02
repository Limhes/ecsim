#ifndef ELECTRODES_H
#define ELECTRODES_H

#include <iostream>
#include <cmath>
#include <vector>

const std::string NO_GEOM = "(n/a)";

class Electrode
{
public:
    static const std::vector<std::string> electrodeTypes;
    static const std::vector<std::string> electrodeGeom1;
    static const std::vector<std::string> electrodeGeom2;

    // electrode parameters:
    double electrodeArea, electrodeGeometryFactor; // electrode area [m2] and geometry (0.0 for planar, 1.0 for cylindrical, 2.0 for spherical)
    double epsilon;

    double m_geom1, m_geom2;
    std::size_t m_type;

    Electrode() { electrodeArea = 0.0; electrodeGeometryFactor = 0.0; epsilon = 0.0; m_geom1 = 0.0; m_geom2 = 0.0; m_type = 0; }
    ~Electrode() {}

    void setType(std::size_t idx) { m_type = idx; update(); }

    void setGeom1(double d) { m_geom1 = d; update(); }
    void setGeom2(double d) { m_geom2 = d; update(); }
private:
    void update();
};

#endif // ELECTRODES_H
