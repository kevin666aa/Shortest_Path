AR=ar
AS=as
CC=gcc
CPP=cpp
CXX=g++
LD=ld
OBJCOPY = objcopy
OBJDUMP = objdump
STRIP = strip

INC_DIR = ./include
SRC_DIR = ./src
OBJ_DIR = ./obj
BIN_DIR = ./bin
TESTBIN_DIR = ./testbin

#Uncomment the next line to enter debug
#DEBUG_MODE=TRUE

PKGS = expat

ifdef DEBUG_MODE
DEFINES  += -DDEBUG
CFLAGS   += -g -ggdb -D_GLIBCXX_DEBUG
endif



INCLUDE  += -I $(INC_DIR)
CFLAGS   +=  -Wall `pkg-config --cflags $(PKGS)`
LDFLAGS  += `pkg-config --libs $(PKGS)`
TESTLDFLAGS  +=  -lgtest -lgtest_main  `pkg-config --libs $(PKGS)` -lpthread
CPPFLAGS += -std=c++14 
PROJ_NAME = proj4
TEST_NAME = test

PROJ_OBJS = $(OBJ_DIR)/XMLParser.o      \
            $(OBJ_DIR)/MapRouter.o      
MAIN_OBJS = $(OBJ_DIR)/main.o
TEST_OBJS = $(OBJ_DIR)/tests.o
            
all: directories test $(BIN_DIR)/$(PROJ_NAME)
test: $(TESTBIN_DIR)/$(TEST_NAME)
	$(TESTBIN_DIR)/$(TEST_NAME)

$(BIN_DIR)/$(PROJ_NAME): directories $(PROJ_OBJS) $(MAIN_OBJS)
	$(CXX) $(PROJ_OBJS) $(MAIN_OBJS) -o $(BIN_DIR)/$(PROJ_NAME) $(CFLAGS) $(CPPFLAGS) $(DEFINES) $(LDFLAGS) 
		
$(TESTBIN_DIR)/$(TEST_NAME): directories $(PROJ_OBJS) $(TEST_OBJS)
	$(CXX) $(PROJ_OBJS) $(TEST_OBJS) -o $(TESTBIN_DIR)/$(TEST_NAME) $(CFLAGS) $(CPPFLAGS) $(DEFINES) $(TESTLDFLAGS) 

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(DEFINES) $(INCLUDE) -c $< -o $@
	
.PHONY: directories
directories:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(BIN_DIR)
	mkdir -p $(TESTBIN_DIR)
	
.PHONY: clean_directories
clean_directories:
	rm -rf $(OBJ_DIR)
	rm -rf $(BIN_DIR)
	rm -rf $(TESTBIN_DIR)
	
.PHONY: clean
clean: clean_directories directories
	
