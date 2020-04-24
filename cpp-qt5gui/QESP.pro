QT += core gui
greaterThan(QT_MAJOR_VERSION, 4): QT += widgets printsupport

TARGET = QESP
DESTDIR = build
TEMPLATE = app
RESOURCES += resources.qrc

CONFIG += debug
#CONFIG += release
CONFIG += c++11

DEFINES += QT_DEPRECATED_WARNINGS

SOURCES += \
    coefs_alpha_beta.cpp \
    electrodes.cpp \
    experiment.cpp \
    main.cpp \
    mainwindow.cpp \
    model.cpp \
    qcustomplot.cpp \
    simulation.cpp \
    system.cpp

HEADERS += \
    coefs_alpha_beta.h \
    delegates.h \
    electrodes.h \
    environment.h \
    experiment.h \
    mainwindow.h \
    model.h \
    qcustomplot.h \
    simulation.h \
    system.h

FORMS += \
    mainwindow.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
