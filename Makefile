# Compile the main executable
ex: cycle_length.cpp
	g++ -Wall -Werror -pedantic -g --std=c++11 cycle_length.cpp -o cycl.exe
fn: functional_graph.cpp
	g++ -Wall -Werror -pedantic -g --std=c++11 functional_graph.cpp -o func.exe

# Remove automatically generated files
clean :
	rm -rvf *.exe *~ *.out *.dSYM *.stackdump
