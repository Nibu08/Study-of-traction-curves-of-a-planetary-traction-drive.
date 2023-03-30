function [rho]=findrfk(Rax,Rbx,Ray,Rby)
rho=1/Rax+1/Rbx+1/Ray+1/Rby;
F=(1/Rax+1/Rbx-(1/Ray+1/Rby))/rho;
disp(F);
