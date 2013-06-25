function [objctv_val, objctv_grad] = objctv1_decomp(...
    x, x1x2, P_dmnt_egnvctrs, P_dmnt_egnvals, biases)

grads = zeros(size(P_dmnt_egnvctrs{1},1),2);
dstncs_sqrd = zeros(1,2);

for k=1:2
    
    [dstncs_sqrd(k), grads(:,k)] = dstnc_objctv_decomp(...
        x-x1x2{k}(:), P_dmnt_egnvctrs{k}, P_dmnt_egnvals{k}, biases(k));
    
end

%objective = (dstnc1 - dstnc2)^2
%grad = 2*(dstnc1 - dstnc2)*(grad1 - grad2);

dstnc_diff = (sqrt(dstncs_sqrd(1)) - sqrt(dstncs_sqrd(2)));
objctv_val = dstnc_diff^2;
objctv_grad = 2*dstnc_diff*(grads(:,1) - grads(:,2));

% dx1 = x - x1x2{1};
% dx2 = x - x1x2{2};
% 
% 
% P1dx1 = P1\dx1;
% P2dx2 = P2\dx2;
% dx1P1dx1 = dx1'*P1dx1;
% dx2P2dx2 = dx2'*P2dx2;
% 
% %objective = (dstnc1 - dstnc2)^2
% objctv_val  = dx1P1dx1 - dx2P2dx2;
% objctv_grad = (4*(objctv_val))*(P1dx1 - P2dx2);
% objctv_val  = objctv_val*objctv_val;

end
