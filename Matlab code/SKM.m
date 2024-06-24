function [t]= SKM(x,y,z)
t = [0 -z y;
     z  0 -x;
     -y x  0];
end
