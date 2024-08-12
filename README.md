# 14.5 Cosmological Distances

## Description
Files contain Python code that analyse the questions for project 14.5 that require programming. The files used are detailed below:

### lookback.py
Computes the lookback time for the 4 different universe cases for a range of redshifts. Outputs a graph and values at specific points of redshift that have beeb tabulated in the report.

### distances.py
Computes the dimensionless luminosity and angular diameter distances as functions of redshift over the range 0 to 7. Depending on the dimensionless parameters from the Friedmann equation, different angular diameter distances are used. Ouptuts a graph and values for tabulation.

### comoving.py
Computes the average value of the V/V_max test. First the program reads and stores the quasar.dat file as arrays. Uses numerical integration of the luminosity distance used to find the comoving volume for z. By looping through values of possible z_max, the program finds the solution using f/f0 and luminosity distances to find D_Lmax, z_max and therefore V_max.