function [x] = MDP_ObjectiveFunction(w2,std,S)
f = -(w2'*std)*(w2'*S*w2)^(-.5);
end