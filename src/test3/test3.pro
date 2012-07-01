TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    test3.cpp \
    ../matrixmaker.cpp \
    ../ircmatrix.cpp \
    ../cholincpreconditioner.cpp

DEFINES = DEBUG
INCLUDEPATH += ../../includes

HEADERS += \
    ../../includes/cholincpreconditioner.h
