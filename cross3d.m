function C = cross3d(A, B)
 
C = [A(2)*B(3) - A(3)*B(2); ...
     A(3)*B(1) - A(1)*B(3); ...
     A(1)*B(2) - A(2)*B(1)];
 
end