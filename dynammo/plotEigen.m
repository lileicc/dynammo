function plotEigen(A, varargin)

eigenvalue = eig(A);
polar(angle(eigenvalue), abs(eigenvalue), varargin{:});