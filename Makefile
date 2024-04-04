CC= g++
CFLAGS= -O3 -fPIC -I/Users/bronsonjeong/Documents/c++libraries/boost_1_78_0 -lrt -Wall -Wextra -w -pthread -std=c++11
#CFLAGS= -O3 -fPIC -I/Users/bronsonjeong/Documents/c++libraries/boost_1_78_0 -Wall -Wextra -w -std=c++17 -lpthread

all: cadmix_mt

cadmix_mt: admix_mt.cpp
		$(CC) $(CFLAGS) -o cadmix_mt admix_mt.cpp

.PHONY: cadmix_mt
