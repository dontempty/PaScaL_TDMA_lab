# =======================================
#            Makefile for MPI project
# =======================================

# 컴파일러 및 플래그
CXX        = g++
CXXFLAGS   = -O2 -std=c++17 -static-libgcc -static-libstdc++
INCLUDES   = -I"C:/Program Files (x86)/Microsoft SDKs/MPI/Include"
LDFLAGS    = -L"C:/Program Files (x86)/Microsoft SDKs/MPI/Lib/x64" -lmsmpi

# 타겟 및 소스
TARGET     = main_lab.exe
SRCS       = main_lab.cpp global.cpp mpi_topology.cpp mpi_subdomain.cpp
OBJS       = $(SRCS:.cpp=.o)

# MPI 실행 옵션
MPIEXEC    = mpiexec
NPROCS     = 4
INPUT_FILE = ../run/PARA_INPUT.txt

.PHONY: all run clean

# 기본 빌드
all: $(TARGET)

# 링크
$(TARGET): $(OBJS)
	@echo "[Link] $@"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

# 컴파일 규칙 (헤더 의존성 포함)
main_lab.o: main_lab.cpp global.hpp mpi_topology.hpp mpi_subdomain.hpp
	@echo "[Compile] $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

global.o: global.cpp global.hpp
	@echo "[Compile] $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

mpi_topology.o: mpi_topology.cpp mpi_topology.hpp
	@echo "[Compile] $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

mpi_subdomain.o: mpi_subdomain.cpp mpi_subdomain.hpp global.hpp mpi_topology.hpp
	@echo "[Compile] $<"
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# MPI 실행
run: all
	@echo "[Run] $(MPIEXEC) -n $(NPROCS) ./$(TARGET) $(INPUT_FILE)"
	$(MPIEXEC) -n $(NPROCS) ./$(TARGET) $(INPUT_FILE)

# 오브젝트 및 실행파일 삭제
clean:
	@echo "[Clean] removing objects and executable"
	rm -f $(OBJS) $(TARGET)
