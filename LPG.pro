TEMPLATE = lib
DEFINES += LPG_LIBRARY
CONFIG += c++11
QT += core
# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

INCLUDEPATH += /usr/include

SOURCES += \
    src/lpg.cpp \
    src/utils.cpp

HEADERS += \
    include/LPG_global.h \
    include/lpg.h \
    include/utils.h

# Default rules for deployment.
unix {
    target.path = /usr/lib
}
!isEmpty(target.path): INSTALLS += target
