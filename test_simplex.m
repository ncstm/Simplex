function test_simplex = test_simplex()

SubjectToEquation = [1 2 3 3 4 -5 2 4 5 -4];

Coefficients = [
    1 4 0 6 4 3 2 2 6 5;
    2 4 6 -4 7 6 -3 -1 4 6;
    8 2 0 -5 -4 9 7 4 -4 9;
    1 0 2 3 4 5 -3 -5 0 6;
    -2 -3 -5 -6 10 0 0 0 -4 -9;
    3 1 2 5 2 4 0 5 7 9;
    3 5 6 6 -4 3 5 0 3 -4;
    4 5 1 2 3 5 -5 -3 -5 -2;
    1 3 3 4 9 0 5 6 7 9;
    2 6 4 2 2 0 -1 4 -5 -6
];

RHS = [25; 15; 16; -16; -3; 5; 14; 15; 4; 10];

Equivalencies = [-1; 0; -1; 1; -1; 1; -1; -1; 1; 0];

Min = true;
constrains = [0 SubjectToEquation; RHS Coefficients];
tic
result = simplex(standard_form(constrains,Equivalencies, Min));
toc
disp(result)
