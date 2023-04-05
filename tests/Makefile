CXXSTD   := c++17

CXX      := g++
CXXFLAGS  = -std=$(CXXSTD) -fno-rtti \
			-pedantic \
			-g -Og \
			-Wno-unused-variable \
			-Wno-unused-but-set-variable \
			-Wno-parentheses \
			-Wall \
			-MMD -MP \
			-I ../
LDFLAGS   =
LDLIBS    = $(LDFLAGS) -lstdc++ -lm

TEST_SRC := test.cpp test_parts.cpp test_pid.cpp test_constexpr.cpp test_bounded_pid.cpp test_disc_diff.cpp
TEST_OBJ := $(TEST_SRC:.cpp=.o)
TESTS    := $(TEST_OBJ:.o=)
DEPENDS  := $(TEST_SRC:.cpp=.d)

all: $(TESTS)
build: all

$(TEST_OBJ): $(TEST_SRC)

$(TEST_OBJ): Makefile

test: test.o
test_parts: test_parts.o
test_pid: test_pid.o
test_constexpr: test_constexpr.o
test_bounded_pid: test_bounded_pid.o
test_disc_diff: test_disc_diff.o

-include $(DEPENDS)
%.o: %.c Makefile
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP -c $< -o $@

clean:
	rm -f $(TEST_OBJ) $(DEPENDS) $(TESTS)

rebuild: clean all

run: all
	@./test

.PHONY: all clean rebuild run
