TEMPLATE = app
CONFIG += console
CONFIG -= qt

HEADERS += \
    ../../includes/cholesky.h

SOURCES += \
    test4.cpp \
    ../matrixmaker.cpp \
    ../ircmatrix.cpp \
    ../cholesky.cpp



INCLUDEPATH = ../../includes
DEFINES = DEBUG
