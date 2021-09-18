#-------------------------------------------------
#
# Project created by QtCreator 2016-10-17T23:19:28
#
#-------------------------------------------------
QMAKE_LFLAGS += -pthread    # без неё не работают потоки STL

QT += core gui
#QT += datavisualization
QT += charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Solid.run
TEMPLATE = app


# для intel mkl
INCLUDEPATH += "/usr/include/mkl"
QMAKE_LFLAGS += -m64
QMAKE_CXXFLAGS += -m64
QMAKE_LIBS += -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -lm -ldl



#debug {
#    # Конфигурация debug? Это qmake не волнует! Он всё равно добавит флаги из release! Ату его!
#    CONFIG -= release # Эта строчка уберёт ключ -O2
#
#    QMAKE_CXXFLAGS += -O0 #Или QMAKE_CXXFLAGS_DEBUG
#    QMAKE_CFLAGS += -O0   #Или QMAKE_CFLAGS_DEBUG
#}

#release {
#    # Эти строчки уберут ключ -O2 -------------------------
#    QMAKE_CXXFLAGS_RELEASE -= -O2
#    QMAKE_CFLAGS_RELEASE -= -O2
#
#    QMAKE_CXXFLAGS += -O3  #Или QMAKE_CXXFLAGS_RELEASE
#    QMAKE_CFLAGS += -O3    #Или QMAKE_CFLAGS_RELEASE
#}
#QMAKE_CFLAGS_RELEASE -= -O2
#QMAKE_CFLAGS_RELEASE += -O3

   QMAKE_CXXFLAGS_RELEASE -= -O2
   QMAKE_CFLAGS_RELEASE -= -O2
   #QMAKE_CXXFLAGS -= -O2
   #QMAKE_CFLAGS -= -O2
   QMAKE_CXXFLAGS += -O3  #Или QMAKE_CXXFLAGS_RELEASE
   QMAKE_CFLAGS += -O3    #Или QMAKE_CFLAGS_RELEASE


QMAKE_CXXFLAGS += -march=native -mtune=native
QMAKE_CFLAGS += -march=native -mtune=native
#QMAKE_CXXFLAGS += -march=native -mtune=native -funroll-loops
#QMAKE_CFLAGS += -march=native -mtune=native -funroll-loops
#QMAKE_CXXFLAGS += -Q -march=native -mtune=native
#QMAKE_CFLAGS += -Q -march=native -mtune=native
#QMAKE_CXXFLAGS += -march=native -funroll-loops
#QMAKE_CFLAGS += -march=native -funroll-loops
#-funroll-loops


#DESTDIR = D:/Projects/SolidGUI/release

#QMAKE_CFLAGS -= -nologo
#QMAKE_CXXFLAGS -= -nologo
#CONFIG += warn_off
#QMAKE_CFLAGS_WARN_ON -= -Wall
#QMAKE_CXXFLAGS_WARN_ON -= -Wall

#QMAKE_CFLAGS_WARN_ON -= -W1
#QMAKE_CFLAGS_WARN_ON -= -W2
#QMAKE_CFLAGS_WARN_ON -= -W3
#QMAKE_CFLAGS_WARN_ON += -W2

#QMAKE_CXXFLAGS += /wd4996
#QMAKE_CXXFLAGS += -implicit conversion changes signedness

#QMAKE_CXXFLAGS += -Wsign-conversion

INCLUDEPATH += $$PWD/elementary

INCLUDEPATH += $$PWD/os
INCLUDEPATH += $$PWD/gui
INCLUDEPATH += $$PWD/sle

INCLUDEPATH += $$PWD/interpolation
INCLUDEPATH += $$PWD/grid

INCLUDEPATH += $$PWD/solving
INCLUDEPATH += $$PWD/solving/heat
INCLUDEPATH += $$PWD/solving/solid


SOURCES += \
    elementary/elementary.cpp \
    elementary/fem.cpp \
    elementary/integration.cpp \
    grid/grid.cpp \
    grid/regionIndexation.cpp \
    gui/visualization2d.cpp \
    interpolation/funParser.cpp \
    interpolation/interpolation.cpp \
    mainwindow.cpp \
    os/console.cpp \
    os/timeinterval.cpp \
    sle/slau3x3.cpp \
    sle/slausolving.cpp \
    solving/heat/heat.cpp \
    solving/logger.cpp \
    solving/solid/post.cpp \
    solving/solid/solid.cpp \
    solving/solid/solidContact.cpp \
    solving/solid/solidPlastic.cpp \
    solving/solid/solid_base.cpp \
    solving/solving.cpp \
    solving/solvingThread.cpp \
    tests.cpp \
    main.cpp \
    tests/test_Kirsch.cpp \
    tests/test_contactCylinder.cpp \
    tests/test_contactSphere.cpp \
    tests/test_crack_plate_gen.cpp \
    tests/test_creep.cpp \
    tests/test_formovka.cpp \
    tests/test_sphere_T.cpp \
    tests/test_sphere_Te.cpp \
    tests/test_sphere_creep.cpp \
    tests/test_sphere_ep.cpp

HEADERS  += \
    __old/_trash.h \
    __old/__todo.h \
    elementary/elementary.h \
    elementary/fem.h \
    elementary/integration.h \
    grid/grid.h \
    grid/grid_base.h \
    grid/regionIndexation.h \
    gui/visualization2d.h \
    interpolation/funParser.h \
    interpolation/interpolation.h \
    mainwindow.h \
    os/console.h \
    os/timeinterval.h \
    sle/slau3x3.h \
    sle/slausolving.h \
    solving/heat/heat.h \
    solving/logger.h \
    solving/logger_base.h \
    solving/solid/post.h \
    solving/solid/solid.h \
    solving/solid/solidContact.h \
    solving/solid/solidPlastic.h \
    solving/solid/solid_base.h \
    solving/solving.h \
    solving/solvingThread.h \
    tests.h \
    tests/test_Kirsch.h \
    tests/test_contactCylinder.h \
    tests/test_contactSphere.h \
    tests/test_crack_plate_gen.h \
    tests/test_creep.h \
    tests/test_formovka.h \
    tests/test_sphere_T.h \
    tests/test_sphere_Te.h \
    tests/test_sphere_creep.h \
    tests/test_sphere_ep.h

FORMS    += mainwindow.ui

#CONFIG(debug, debug|release){
#    DESTDIR=$$shadowed($$ROOT_PWD)/debug
#}else{
#    DESTDIR=$$shadowed($$ROOT_PWD)/release
#}

#STATECHARTS +=
