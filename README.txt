ECE6280
Crypto project


instructions 

open a terminal

go to the folder where you unzipped the file

start bash by typing "bash" 

compile with "gcc index_calculus.c -o index_calculus -lgmp -L=. "

run the program with "LD_LIBRARY_PATH=. ./index_calculus  your_number"

The -L and LD_PRIMARY_PATH are some tricks to use libraries and files not in their normal location. In this case, they are in the current folder, which is unusual. 
gmp.h contains the function definitions, libgmp.a the compilation instructions for the library and libgm.so.10.0 is the library itself.




