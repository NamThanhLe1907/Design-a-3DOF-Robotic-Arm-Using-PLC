function [t]= SKMN(B)
x = B(1,1);y = B(2,1);z = B(3,1);
B = [x;y;z];
t = [0 -z y;
     z  0 -x;
     -y x  0];
end