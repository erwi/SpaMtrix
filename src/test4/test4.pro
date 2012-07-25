TEMPLATE = app
CONFIG += console
CONFIG -= qt

HEADERS += \
    ../../includes/cholesky.h \
    ../../includes/fleximatrix.h \
    ../../includes/vector.h \
    ../../includes/spamtrix_blas.h

SOURCES += \
    test4.cpp \
    ../matrixmaker.cpp \
    ../ircmatrix.cpp \
    ../cholesky.cpp \
    ../fleximatrix.cpp \
    ../vector.cpp \
    ../spamtrix_blas.cpp



INCLUDEPATH = ../../includes
DEFINES = DEBUG
