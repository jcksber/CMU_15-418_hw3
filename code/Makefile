APP_NAME=wireroute

OBJS=wireroute.o

default: $(APP_NAME)

# Compile for Xeon Phi
$(APP_NAME): CXX = icc -m64 -std=c++11
$(APP_NAME): CXXFLAGS = -I. -O3 -Wall -openmp -offload-attribute-target=mic -DRUN_MIC

# Compile for CPU
cpu: CXX = g++ -m64 -std=c++11
cpu: CXXFLAGS = -I. -O3 -Wall -fopenmp -Wno-unknown-pragmas

# Compilation Rules
$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

cpu: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(APP_NAME) $(OBJS)

%.o: %.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

submit:
	cd jobs && ./batch_generate.sh && cd ../latedays && ./submit.sh
clean:
	/bin/rm -rf *~ *.o $(APP_NAME) jobs/$(USER)_*.job latedays/$(USER)_*

# For a given rule:
# $< = first prerequisite
# $@ = target
# $^ = all prerequisite

