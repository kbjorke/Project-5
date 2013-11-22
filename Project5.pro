TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    Function.cpp \
    Integral.cpp \
    lib.cpp \
    hermite.cpp \
    gaussiandeviate.cpp \
    Hamiltonian.cpp \
    UnixTime.cpp

HEADERS += \
    Function.h \
    Integral.h \
    lib.h \
    hermite.h \
    gaussiandeviate.h \
    Hamiltonian.h \
    UnixTime.h

LIBS += -lrt
