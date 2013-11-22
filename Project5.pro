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
    UnixTime.cpp \
    problem_definitions.cpp \
    output_functions.cpp

HEADERS += \
    Function.h \
    Integral.h \
    lib.h \
    hermite.h \
    gaussiandeviate.h \
    Hamiltonian.h \
    UnixTime.h \
    problem_definitions.h \
    output_functions.h

LIBS += -lrt
