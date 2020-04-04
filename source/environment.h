#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <iostream>

using namespace std;

class Environment
{
public:
    Environment() : temperature(298.15) {}
    ~Environment() {}
    double temperature;

    void setTemperature(double _t) { temperature = _t; }
};

#endif // ENVIRONMENT_H
