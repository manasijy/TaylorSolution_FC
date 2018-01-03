function [DCM] = DC(g)

DCM(1,1) = g(1);
DCM(1,2) = g(2);
DCM(1,3) = g(3);
DCM(2,1) = g(4);
DCM(2,2) = g(5);
DCM(2,3) = g(6);
DCM(3,1) = g(7);
DCM(3,2) = g(8);
DCM(3,3) = g(9);
