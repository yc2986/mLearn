CC = g++
TARGET = ../mLearn
SOURCE = mLearn.cpp
DEPS = ../lib/libmodel.so
OBJDIR = obj
OBJ = $(addprefix $(OBJDIR)/, mLearn.o)

CPPFLAGS = -w -std=c++11 -I../include -I../ -Ofast -DNDEBUG
LINK = -Wl,--no-as-needed -L../lib -lmodel

$(TARGET): $(OBJ)
	$(CC) $(OBJ) $(LINK) -o $(TARGET)

$(OBJ): $(SOURCE)
	$(CC) $(CPPFLAGS) -c $(SOURCE) -o $(OBJ)

.PHONY: clean
clean:
	rm $(TARGET) $(OBJ)