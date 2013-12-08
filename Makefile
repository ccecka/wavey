#
# 'make'        build executable file
# 'make clean'  removes all .o and executable files
#

DEPSDIR := $(shell mkdir -p .deps/; echo .deps)
# Get the shell name to determine the OS
UNAME := $(shell uname)

# define the C compiler to use
CC := gcc
ifeq ($(UNAME), Linux)
CXX := g++ -std=gnu++0x
endif
ifeq ($(UNAME), Darwin)
CXX := $(shell for i in 4.7 4.6 4.5; do if g++-mp-$$i -v >/dev/null 2>&1; then echo g++-mp-$$i; exit; fi; done; echo false) -std=gnu++0x
OBJC := gcc
endif
LINK := $(CXX)

# define any compile-time flags
CFLAGS := -O3 -W -Wall -pedantic #-Wextra
ifeq ($(PROFILE),1)
CFLAGS += -g -pg
endif
DEPCFLAGS = -MD -MF $(DEPSDIR)/$*.d -MP

# define any directories containing header files other than /usr/include
#   include directories like -Ipath/to/files
INCLUDES = -I. -I/usr/local/include

# define any libraries to link into executable
#   To link in libraries (libXXX.so or libXXX.a) use -lXXX options
LIBS += -lfftw3 -lgsl -lgslcblas -lm
#LIBS += -lfftw3 -llapack -lblas -lgsl -lgslcblas -lm

##################
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
##################

# suffix replacement rule for building .o's from .cpp's
#   $<: the name of the prereq of the rule (a .cpp file)
#   $@: the name of the target of the rule (a .o file)
.cpp.o:
	$(CXX) $(CFLAGS) $(DEPCFLAGS) $(DEFS) $(INCLUDES) -c -o $@ $<

# 'make' - default rule
all: mlfmm

# make sure object files are up-to-date, then compile MAIN
mlfmm: mlfmm.o
	$(LINK) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

time_interp: time_interp.o
	$(LINK) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

time_uniform: time_uniform.o
	$(LINK) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

time_sphere: time_sphere.o
	$(LINK) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

kappa_error: kappa_error.o
	$(LINK) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LIBS)

# 'make clean' - deletes all .o and temp files, exec, and dependency file
clean:
	-$(RM) *.o *~
	-$(RM) mlfmm time_interp time_uniform time_sphere kappa_error
	$(RM) -r $(DEPSDIR)

always:
	@:

DEPFILES := $(wildcard $(DEPSDIR)/*.d) $(wildcard $(DEPSDIR)/*/*.d)
ifneq ($(DEPFILES),)
include $(DEPFILES)
endif

# define rules that do not actually generate the corresponding file
.PHONY: clean all always
