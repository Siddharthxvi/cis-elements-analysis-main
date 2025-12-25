# Compiler and compilation flags
CXX = g++
CXXFLAGS = --std=c++17 -Wall -Wextra -O2

# Source files
SOURCES = Main.cpp Analysis.cpp FindPromoter.cpp FindSpacers.cpp

# Object files (derived from source files)
OBJECTS = $(SOURCES:.cpp=.o)

# Target executable name
TARGET = Main

# Default rule: Build the executable
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)

# Phony targets to avoid filename conflicts
.PHONY: all clean
