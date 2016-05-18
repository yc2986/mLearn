GARBAGE_OBJ = $(shell find . -name '*.o')
GARBAGE_SO = $(shell find . -name '*.so')
EXEC = mLearn

all:
	+$(MAKE) -C source
	+$(MAKE) -C test_app

.PHONY: clean
clean:
	rm $(GARBAGE_OBJ) $(GARBAGE_SO) $(EXEC)