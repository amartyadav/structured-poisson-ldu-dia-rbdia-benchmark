CC = g++
CFLAGS_BASE = -Wall -Wextra -g -std=c++17 -I./include -fopenmp
LDFLAGS = -fopenmp -lm

SRCS = src/dia_solver.cpp src/rb_dia_solver.cpp src/ldu_solver.cpp src/main.cpp src/problem.cpp src/utils.cpp src/math_sanity_check.cpp

OBJS = $(SRCS:.cpp=.o)

TARGET = poisson_benchmark

# Profiling object files for each optimisation level
OBJS_O1 = $(SRCS:.cpp=.O1.o)
OBJS_O2 = $(SRCS:.cpp=.O2.o)
OBJS_O3 = $(SRCS:.cpp=.O3.o)

TARGET_O1 = poisson_benchmark_O1
TARGET_O2 = poisson_benchmark_O2
TARGET_O3 = poisson_benchmark_O3

.PHONY: all profile clean

all: $(TARGET)

profile: $(TARGET_O1) $(TARGET_O2) $(TARGET_O3)

$(TARGET): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

$(TARGET_O1): $(OBJS_O1)
	$(CC) $(OBJS_O1) -o $@ $(LDFLAGS) -llikwid

$(TARGET_O2): $(OBJS_O2)
	$(CC) $(OBJS_O2) -o $@ $(LDFLAGS) -llikwid

$(TARGET_O3): $(OBJS_O3)
	$(CC) $(OBJS_O3) -o $@ $(LDFLAGS) -llikwid

%.o: %.cpp
	$(CC) $(CFLAGS_BASE) -O2 -c $< -o $@

%.O1.o: %.cpp
	$(CC) $(CFLAGS_BASE) -O1 -DLIKWID_PERFMON -c $< -o $@

%.O2.o: %.cpp
	$(CC) $(CFLAGS_BASE) -O2 -DLIKWID_PERFMON -c $< -o $@

%.O3.o: %.cpp
	$(CC) $(CFLAGS_BASE) -O3 -march=native -DLIKWID_PERFMON -c $< -o $@

clean:
	rm -rf $(TARGET) $(TARGET_O1) $(TARGET_O2) $(TARGET_O3) \
		$(OBJS) $(OBJS_O1) $(OBJS_O2) $(OBJS_O3)
