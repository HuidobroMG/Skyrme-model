all: minim6_cos.exe


COPTS = -O5 -Wall -W -Wshadow -Wcast-qual -Wcast-align -Wconversion \
	-Winline -Wold-style-cast -Woverloaded-virtual -Wno-unused
	
minim6_cos.exe: minim6_cos.cc  fields_6_cos.c
	g++ ${COPTS} -o minim6_cos.exe minim6_cos.cc  fields_6_cos.c \
	-L/lib64/ -lgsl -lgslcblas


.compress: 
	touch .compress; ./save.com

.PHONY: clean 
clean: 
	rm -i *.exe *~ *.*~ .compress

