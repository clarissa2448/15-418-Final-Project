APP_NAME=luby_mpi
OBJS += luby_mpi.o
OBJS += luby_mpi_blocked_assignment.o
OBJS += luby_mpi_blocked_every_pair.o
OBJS += util.o
OBJS += luby_mpi_blocked_assignment_async.o
OBJS += main.o

CXX = mpic++ -std=c++11
CXXFLAGS = -I. -O3 #-Wall -Wextra

default: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp
	$(CXX) $< $(CXXFLAGS) -c -o $@

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class