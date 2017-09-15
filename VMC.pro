TEMPLATE = subdirs
CONFIG += ordered
SUBDIRS +=
    src
    app
    tests

tests.depends = src
app.depends = src

OTHER_FILES += \
    defaults.pri
