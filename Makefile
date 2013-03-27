CC = g++
CFLAGS = -w
LIBS = -lm

.SUFFIXES: .cpp .o

#Dependency Format:
#  Target: Source
#  [TAB]Command

#Pre Defined Macros:
#  $(CC) $(CFLAGS) ...
#  $@ :full name of the current target.
#  $< :source file of the current target. 
#  $? :list of updated dependency files.

OBJS = IOFile.o LOGFile.o DOS.o main.o
defalt: plotDOS.exe

plotDOS.exe: $(OBJS)
	$(CC) $(OBJS) -o $@
.cpp.o:
	$(CC) $(CFLAGS) -c $< 

#just to exe command
.PHONY: clear clean
cl clear clean:
	rm -vf *.o *.exe

