#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <QObject>
#include <iostream>

using namespace std;

class Environment : public QObject
{
    Q_OBJECT

public:
    Environment(double _t) { temperature = _t; }
    ~Environment() {}
    double temperature;

public slots:
    void setTemperature(double _t) { temperature = _t; }
};

#endif // ENVIRONMENT_H
