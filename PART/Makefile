# Collect all libraries
LIBS = -larmadillo -lopenblas -larpack -lsuperlu
ARMADILLO_LIB = $(mkArmadilloLib)
ARPACK_LIB = $(mkArpackLib)
OPENBLAS_LIB = $(mkOpenblasLib)
SUPERLU_LIB = $(mkSuperluLib)
ALL_LIB_DIR =  $(ARMADILLO_LIB) $(ARPACK_LIB) $(OPENBLAS_LIB) $(SUPERLU_LIB)
libraries := $(foreach lib,$(ALL_LIB_DIR),-L$(lib)) 
libraries_rpath := $(foreach lib,$(ALL_LIB_DIR),-Wl,-rpath=$(lib)) 
LINK_FLAGS = $(libraries) $(LIBS) $(libraries_rpath)

# Collect all includes
ARMADILLO_INC = $(mkArmadilloInc)
OPENBLAS_INC = $(mkOpenblasInc)
SUPERLU_INC = $(mkSuperluInc)
ALL_INC_DIR = $(ARMADILLO_INC) $(OPENBLAS_INC) $(SUPERLU_INC)
includes := $(foreach inc,$(ALL_INC_DIR),-I$(inc))
INCLS  = -I./include -I./include/sigpack $(includes)

# C++ compiler and flags
CXX = g++
STDFLAGS = -std=c++11
WARNFLAGS = -Wall
OPTFLAGS_REL = -O2
OPTFLAGS_DEB = -g -O0
OPTFLAGS = $(OPTFLAGS_REL)
UNIXFLAG = -Dunix
CXXFLAGS = $(STDFLAGS) $(WARNFLAGS) $(UNIXFLAG)

# Location of Doxygen configuration file
DOC_DIR = ../doc

# Source directory
SRCS_DIR = ./src

# Create object .o from source
OBJS=$(subst .cpp,.o,$(wildcard $(SRCS_DIR)/*.cpp))

# Create exe from object
exe_obj=$(filter $(SRCS_DIR)/main%.o,$(OBJS))
EXEC=$(exe_obj:$(SRCS_DIR)/%.o=%)

# Phony targets
.PHONY = check-env all example real doxy clean distclean debug prof

# "make" does "make all" by deafault
.DEFAULT_GOAL = all

# check if environmental variables have been defined, otherwise error
ENV_VARS := mkArmadilloLib mkArpackLib mkOpenblasLib mkSuperluLib mkArmadilloInc mkOpenblasInc mkSuperluInc
fatal_if_undefined = $(if $(findstring undefined,$(origin $1)),$(error Error: variable [$1] is undefined))
$(foreach 1,$(ENV_VARS),$(fatal_if_undefined))

# make example and real
all: CXXFLAGS += $(OPTFLAGS)
all: $(EXEC)

# debug mode
debug: CXXFLAGS += $(OPTFLAGS_DEB)
debug: $(EXEC)

# profile mode
prof: CXXFLAGS += -pg
prof: all

# create documentation with doxygen
# @ >/dev/null 2>&1 suppress command output and errors (there are too many lines)
doxy:
	@doxygen $(DOC_DIR)/Doxyfile
	@$(MAKE) -C $(DOC_DIR)/latex
	@cp $(DOC_DIR)/latex/refman.pdf $(DOC_DIR)
	@mv $(DOC_DIR)/refman.pdf $(DOC_DIR)/DocumentationPART_Raciti_Riva.pdf
	
# remove all .o and exe
clean:
	$(RM) -f $(EXEC) $(OBJS)

# remove all files generated
distclean:
	$(MAKE) clean
	$(RM) -r $(DOC_DIR)/latex $(DOC_DIR)/html $(DOC_DIR)/doxygen_sqlite3.db
	$(RM) *.out *.bak *~

# Compile
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLS) $(OBJS) -o $@ $(LINK_FLAGS)

# Link
$(SRCS_DIR)/%.o: $(SRCS_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLS) $< -c -o $@