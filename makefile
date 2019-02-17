CC = g++ -std=c++14 -Wall -pedantic
TARGET:=siberneic
TEST_TARGET:=stest
RM := rm -rf
TEST := test
SRC_DIR := src
INC_DIR := inc 
BUILD_DIR := bin
TEST_BIN_DIR := $(TEST)/bin
SRC_EXT := cpp
BINARY_DIR = $(BUILD_DIR)/obj


# OCL_INC  = -I/opt/AMDAPPSDK-3.0/include/
# OCL_INC  = -I/opt/intel/opencl/SDK/include
# OCL_LIB  = -L/usr/lib64/OpenCL/vendors/amd/

LIBS := -lOpenCL
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S), Darwin)
	LIBS := -framework OpenCL
endif


SRC := $(wildcard $(SRC_DIR)/*.$(SRC_EXT))
SRC += $(wildcard $(SRC_DIR)/utils/*.$(SRC_EXT))
TEST_SRC := $(wildcard $(TEST)/*.$(SRC_EXT))

OBJ := $(patsubst $(SRC_DIR)/%,$(BINARY_DIR)/%,$(SRC:.$(SRC_EXT)=.o))
TEST_OBJ := $(patsubst $(TEST)/%,$(TEST_BIN_DIR)/%,$(TEST_SRC:.$(SRC_EXT)=.o))

CXXFLAGS = $(CC)

all: CXXFLAGS += -O3
all: $(TARGET)

debug: CXXFLAGS += -ggdb -O0
debug: $(TARGET)

test: SRC := $(shell ls $(SRC_DIR)| grep [^sibernetic]\.cpp)
test: SRC := $(wildcard $(SRC_DIR)/$(SRC))
test: OBJ := $(patsubst $(SRC_DIR)/%,$(BINARY_DIR)/%,$(SRC:.$(SRC_EXT)=.o))
test: $(TEST_TARGET)


$(TARGET): $(OBJ)
	@echo 'Building target: $@'
	$(CXXFLAGS) $(OCL_LIB) -o $(BUILD_DIR)/$(TARGET) $(OBJ) $(LIBS)
	@echo 'Finished building target: $@'

$(BINARY_DIR)/%.o: $(SRC_DIR)/%.$(SRC_EXT)
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BINARY_DIR)
	@mkdir -p $(BINARY_DIR)/utils
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXXFLAGS) $(OCL_INC) -I$(INC_DIR) -c -o "$@" "$<"
	@echo 'Finished building: $<'

$(TEST_BIN_DIR)/%.o: $(TEST)/%.$(SRC_EXT)
	@mkdir -p $(TEST_BIN_DIR)
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXXFLAGS) -I$(INC_DIR) -c -o "$@" "$<"
	@echo 'Finished building: $<'


$(TEST_TARGET): $(OBJ) $(TEST_OBJ)
	@echo ' $(SRC)'
	@echo ' $(OBJ)'
	@echo 'Building target: $@'
	$(CXXFLAGS) $(OCL_LIB) -o $(TEST_BIN_DIR)/$(TEST_TARGET) $(TEST_OBJ) $(OBJ) $(LIBS) -lgtest
	@echo 'Finished building target: $@'

clean:
	@echo "Cleaning...";
	$(RM) $(BUILD_DIR)
	$(RM) $(TEST_BIN_DIR)

radio: 
	@echo 'Building files: $(CPP_FILES)'

.PHONY: all clean debug mac_os test
