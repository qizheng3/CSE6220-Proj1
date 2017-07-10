% random number generator for input
clc;
clear;

n2 = 5;
m = 100;
n = 2^n2;
n = 5000000;
constants = 1.0*(rand(n,1)-0.5);
values = 1+0.0001*(rand(m,1)-0.5);

file1 = fopen('sample-constants.txt','w');
fprintf(file1,'%d\n',n);
fprintf(file1,'%8.9f\n',constants);
fclose(file1);

file2 = fopen('sample-values.txt','w');
fprintf(file2,'%d\n',m);
fprintf(file2,'%8.9f\n',values);
fclose(file2);
