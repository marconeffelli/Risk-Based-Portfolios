function [x] = ERC_ObjectiveFunction(w1,S)
x = 0;
R = S*w1;
for i=1:size(w1)
    for j=1:size(w1)
        x = x + (w1(i)*R(i)-w1(j)*R(j))^2;
    end
end
x = x/(w1'*R);
end