function [ u ] = pipeFlowProfileOptimizer( r,rPipe,u0,a )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
u = 2*u0*(1-(r/rPipe).^a);

end

