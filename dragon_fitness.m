function [f,acceptance_probablity] = dragon_fitness(a,c,s)

%DRAGON_FITNESS Summary of this function goes here
%   Detailed explanation goes here
%a stands for allignment
%s stands for seperation 
%c stands for cohesion 
f=0;
atc_diff=0; % allignment to cohesion 
atc_diff=(abs(a-c)/c)*100;
ats_diff=(abs(a-s)/s)*100;
if atc_diff<30
   if ats_diff>40
      f=1;
      
   end
end
acceptance_probablity=1-atc_diff/100;

end

