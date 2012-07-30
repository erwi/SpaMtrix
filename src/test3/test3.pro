TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += \
    test3.cpp \
    ../matrixmaker.cpp \
    ../ircmatrix.cpp \
    ../cholincpreconditioner.cpp \
    ../cholesky.cpp \
    ../fleximatrix.cpp \
    ../vector.cpp

DEFINES = DEBUG
INCLUDEPATH += ../../include

HEADERS += \
    ../../include/cholincpreconditioner.h \
    ../../include/cholesky.h \
    ../../include/ircmatrix.h
