TEMPLATE = app
CONFIG += console
CONFIG -= qt

HEADERS += \
    ../../includes/cholesky.h \
    ../../includes/fleximatrix.h \
    ../../includes/vector.h

SOURCES += \
    test4.cpp \
    ../matrixmaker.cpp \
    ../ircmatrix.cpp \
    ../cholesky.cpp \
    ../fleximatrix.cpp \
    ../vector.cpp



INCLUDEPATH = ../../includes
DEFINES = DEBUG
