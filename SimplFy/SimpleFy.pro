TARGET = simpleFy
DEPENDPATH +=
INCLUDEPATH += libs/vcglib/ #Specify the vcglib diretory location here!
CONFIG += console stl debug_and_release
QMAKE_CXXFLAGS += -fpermissive -std=c++11
TEMPLATE = app
HEADERS += \
    src/simplifier/base.h \
    src/simplifier/clustering.h \
    src/simplifier/mesh_properties.h \
    src/simplifier/quadricdecimator.h \
    src/file_handler.h
SOURCES += \
    libs/vcglib/wrap/ply/plylib.cpp \
    src/main.cpp

# Mac specific Config required to avoid to make application bundles
CONFIG -= app_bundle
