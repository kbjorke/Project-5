# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

# QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK

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
    output_functions.cpp \
    VariationalMC.cpp

HEADERS += \
    Function.h \
    Integral.h \
    lib.h \
    hermite.h \
    gaussiandeviate.h \
    Hamiltonian.h \
    UnixTime.h \
    problem_definitions.h \
    output_functions.h \
    VariationalMC.h

LIBS += -lrt
