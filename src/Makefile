LIBNAME = GeFiCa
LIBRARY = lib$(LIBNAME).so
ROOTMAP = $(LIBRARY:.so=.rootmap)
DIC_SRC = $(LIBNAME).cc

LINKDEF = LinkDef.h
SOURCES = $(filter-out $(DIC_SRC), $(wildcard *.cc))
HEADERS = $(filter-out $(LINKDEF), $(wildcard *.h))
OBJECTS = $(SOURCES:.cc=.o)
DEPENDS = $(SOURCES:.cc=.d)

RC = root-config
CXX = $(shell $(RC) --cxx)
LIBS = $(shell $(RC) --libs)
CXXFLAGS = $(shell $(RC) --cflags) -O2 -Wall -fPIC
# ------------------------------------------------------------------------------
# first target is the default target, before "make all", make will include all
# other makefiles specified by the include command
all: $(LIBRARY)
	@echo
	@echo "Please add current folder to (DY)LD_LIBARY_PATH."
	@echo 

# $(LIBRARY) depends on all *.o files.
#  The flag "-shared" is used to create shared libs
#  $@ represents the target, that is, $(LIBRARY)
#  $^ represents all the prerequisites, i.e., all *.o files
$(LIBRARY): $(DIC_SRC:.cc=.o) $(OBJECTS)
	@echo
	@echo "creating shared library:"
	$(CXX) $^ $(CXXFLAGS) -shared $(LIBS) -o $@

# An xxx.o file depends on xxx.cc. It is created with the command:
# 	$(CXX) -c xxx.cc -o xxx.o
# Since this is obvious, "make" automatically does it. 
# There is no need to explicitly write down the rules to do it.

# create dictionary and rootmap file for ROOT interactive session
$(DIC_SRC): $(HEADERS) $(LINKDEF)
	@echo 
	@echo "creating ROOT dictionary and rootmap files:" 
	rootcling -f $@ $(CXXFLAGS) -s lib$(LIBNAME) -rml lib$(LIBNAME) -rmf $(ROOTMAP) $^
	@echo 
	@echo "creating object files:" 

# include *.d files, which are makefiles defining dependencies between files,
# this will replace the default target if it is placed before "all"
ifeq ($(filter info clean tags, $(MAKECMDGOALS)),)
  -include $(DEPENDS) # - tells make to continue if dependence files do not exist
endif
# rules to create *.d files
%.d:%.cc
	@echo creating $@
	@set -e; rm -f $@; \
	  $(CXX) -MM $(CXXFLAGS) $< > $@.$$$$; \
	  sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	  rm -f $@.$$$$ 
clean:
	rm -f `cat .gitignore`
tags:
	ctags --c-kinds=+p $(HEADERS) $(SOURCES)
.PHONY: all clean tags
