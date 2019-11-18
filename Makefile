# Makefile for TimeDomainTHz
# gcc 4.9.2
# GNU Make version 4.0

SHELL = /bin/sh

SRCDIR = ./src
OBJDIR = ./obj
LIBDIR = ./lib

INCPATH      = -I $(SRCDIR)

CC           = gcc
CXX          = g++
LINK         = g++

CFLAGS       = -O2 -g -Wall -std=c99
CXXFLAGS     = -O2 -g -Wall
LFLAGS       =
LIBS         =

####### Output directory

OBJECTS_DIR   = ./obj

####### Files

HEADER = 	$(SRCDIR)/fields.h \
	$(SRCDIR)/vector.h \
	$(SRCDIR)/screen.h \
	$(SRCDIR)/global.h

SRC = 	$(SRCDIR)/fields.cpp \
	$(SRCDIR)/vector.cpp \
	$(SRCDIR)/screen.cpp \
	$(SRCDIR)/GaussianWavePacket.cpp

OBJ = 	$(OBJDIR)/fields.o \
	$(OBJDIR)/screen.o \
	$(OBJDIR)/vector.o

TARGETOBJ = $(OBJDIR)/GaussianWavePacket.o

TARGET = GaussianWavePacket

INCPATH = -I $(LIBDIR)

####### Implicit rules

.SUFFIXES: .o .c .cpp

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o $@ $<

$(OBJDIR)/%.o: $(LIBDIR)/%.c
	$(CC) -c $(CFLAGS) $(INCPATH) -o $@ $<

####### Build rules

.PHONY : first all clean

first: $(TARGET)

all: $(TARGET)

$(TARGET):  $(OBJ)  $(TARGETOBJ)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJ) $(TARGETOBJ) $(LIBS)

$(OBJ): $(HEADER)

clean:
	-rm $(OBJ)
	-rm $(TARGETOBJ)
	-rm $(TARGET)

