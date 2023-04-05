CXXSTD   := c++17

CXX      := g++
CXXFLAGS  = -std=$(CXXSTD) -fno-rtti -Wno-unused-variable -Wno-unused-but-set-variable -Wno-parentheses -Wall -O3 -I../
LDFLAGS   =
LDLIBS    = $(LDFLAGS) -lstdc++ -lpthread -latomic -lm

TEST_SRC := test.cpp test_parts.cpp test_pid.cpp test_constexpr.cpp test_bounded_pid.cpp
TEST_OBJ := $(TEST_SRC:.cpp=.o)
TESTS    := $(TEST_OBJ:.o=)

all: $(TESTS)
build: all

$(TEST_OBJ): $(TEST_SRC)

$(TEST_OBJ): Makefile

test: test.o
test_parts: test_parts.o
test_pid: test_pid.o
test_constexpr: test_constexpr.o
test_bounded_pid: test_bounded_pid.o

clean:
	rm -f $(TEST_OBJ) $(TESTS)

rebuild: clean all

run: all
	@./test

.PHONY: all clean rebuild run
