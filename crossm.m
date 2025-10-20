%% Cross Matrix Function
% Robert Halverson
% 1/224/2021

function M = crossm(u)

M = [0 -u(3) u(2);u(3) 0 -u(1);-u(2) u(1) 0];