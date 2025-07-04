# ── 컴파일러 및 옵션 ─────────────────────────────────────
CXX       := g++
CXXFLAGS  := -O2 -std=c++17 -static-libgcc -static-libstdc++
# MPI 경로: 공백 없는 8.3 이름 사용
MPI_INC   := -IC:/PROGRA~2/Microsoft\ SDKs/MPI/Include
MPI_LIB   := -LC:/PROGRA~2/Microsoft\ SDKs/MPI/Lib/x64 -lmsmpi

# ── 소스, 오브젝트, 디펜던시 ───────────────────────────────
SRC    := main_lab.cpp \
          global.cpp \
          mpi_topology.cpp \
          mpi_subdomain.cpp \
          solve_theta.cpp \
          ../scr/pascal_tdma.cpp \
		  ../scr/tdmas.cpp \
		  ../scr/para_range.cpp

OBJ    := $(SRC:.cpp=.o)
DEP    := $(SRC:.cpp=.d)
TARGET := main_lab.exe

# ── MPI 실행 설정 ─────────────────────────────────────────
MPIEXEC    := mpiexec
NPROCS     := 4
INPUT_FILE := ../run/PARA_INPUT.txt

.PHONY: all clean run

# ── 기본 빌드 ─────────────────────────────────────────────
all: $(TARGET)

$(TARGET): $(OBJ)
	@echo "[Link] $@"
	$(CXX) $(CXXFLAGS) -o $@ $^ $(MPI_LIB)

# ── 패턴 규칙 + 자동 의존성 생성 ───────────────────────────
%.o: %.cpp
	@echo "[Compile] $<"
	$(CXX) $(CXXFLAGS) $(MPI_INC) -MMD -MP -c $< -o $@

-include $(DEP)

# ── 실행 ───────────────────────────────────────────────────
run: $(TARGET)
	@echo "[Run] $(MPIEXEC) -n $(NPROCS) ./$(TARGET) $(INPUT_FILE)"
	$(MPIEXEC) -n $(NPROCS) ./$(TARGET) $(INPUT_FILE)

# ── 정리 ───────────────────────────────────────────────────
clean:
	@echo "[Clean] Removing objects and dependencies"
	rm -f $(TARGET) $(OBJ) $(DEP)
