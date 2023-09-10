Bethe Ansatz for 1D Hubbard model
=================================

Authors: Gerald Knizia
         Chong Sun ([email](sunchong137@gmail.comn))

Originally written by Gerald Knizia, see `notes_from_knizia.txt`.
This version removed the dependency on boost, and rewrote the code to be
compatible for Python 3.

# Running the code
In this directory, type

`make`

to compile `bethe_ansatz.cpp` file.

Then you can either run the code by

`./bethe_ansatz U Q`

where `U` is the Hubbard U value, and `Q` is the initial guess of the filling. `Q = 1` corresponds to half-filling.

Or you can run the python interface:

`python bethe_ansatz.py`
