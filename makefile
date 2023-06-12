CXX = g++
CXXFLAGS = -fopenmp
OBJS = \
		main.o\
		Random.o\
		class.o\
		ShapeSampling.o\
		restructure.o\
		force.o\
		ODE_solver.o\

all: $(OBJS)
	$(CXX) -g $(CXXFLAGS)  $(OBJS) 

main.o: main.cpp Random.hpp vector3.hpp parameters.hpp class.hpp ShapeSampling.hpp	
Random.o:Random.cpp Random.hpp
class.o:class.cpp class.hpp vector3.hpp  force.hpp
ShapeSampling.o:ShapeSampling.cpp ShapeSampling.hpp Random.hpp ODE_solver.hpp restructure.hpp 
restructure.o:restructure.cpp restructure.hpp Random.hpp class.hpp ODE_solver.hpp
force.o:force.cpp force.hpp  class.hpp vector3.hpp
ODE_solver.o:ODE_solver.cpp ODE_solver.hpp force.hpp class.hpp vector3.hpp

clean:
	rm -f *.o a.out *.vtk d*.off d*.dat e* nohup.out
