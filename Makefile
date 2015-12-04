all: raytrace.cpp
	g++ -O3 -o raytrace raytrace.cpp
clean: 
	rm raytrace