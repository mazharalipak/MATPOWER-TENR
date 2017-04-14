# Matpower-TENR
TENR 
TENR Algorithm
There are three different folders for different choices of Transversality.....
LU,QR and SVD are in one folder (check Readme1 for instruction), (The algorithm is robust and ready to use).
The recommendation is to try LU,QR and SVD as the code is fully debugged and results are better and faster comparing to right_eig and determinant....
Right_eig vector Transversality in separate folder (check Readme2. for instruction)
Determinant Transversality in separate folder (check Readme3. for instruction)
Kindly, also check the files with values of alpha for different choices of Transversality......
Some test case can be initiated with a slightly different "initial state" i.e. not with Lambda^{0}=0 but Lambda^{0}>0.
IEEE 300 case requires a bit tunning for QR and LU transversality except SVD.
