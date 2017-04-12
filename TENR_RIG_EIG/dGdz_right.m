function [ dGdz ] = dGdz_right( J, z )
%dGdz_right Derivative wrt z
%   Detailed explanation goes here
    n = size(J,1);

    if issparse(J)
        zrow = sparse(ones(n,1),1:n,z);
    else
        zrow = z';
    end
    
    dGdz = [J; 2*zrow];
end

