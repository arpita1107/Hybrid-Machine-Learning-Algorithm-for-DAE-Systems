function [H] = electrode_nn_Jacob(h, Z0) 

% <---------------------------------------------------------------------------->
% This function computes numerical gradient matrix for any function vector  
% ------------------------------------------------------------------------------

Z0;

electrode_global;

nZ0 = length(Z0);

for i = 1 : nZ0 
    
    eps = Z0(i) / 1000;
    
    if (eps == 0)
        
        eps = 1 / 1000;
        
    end     
    
    Zp = Z0;    
    Zp(i) = Zp(i) + eps ;
    fp = feval(h, Zp);
    Zn = Z0;
    Zn(i) = Zn(i) - eps ;
    fn = feval(h, Zn)  ; 
    H(:,i) = (fp - fn) / (2 * eps);
    
end     

end