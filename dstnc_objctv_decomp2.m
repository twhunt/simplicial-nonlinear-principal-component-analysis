function [mtrc, mtrc_grad] = dstnc_objctv_decomp2(...
    x, Q_small_egnvctrs, Q_small_egnvals, Q_large_egnval)



x_crndt_vctr = Q_small_egnvctrs'*x(:);
Q_small_egnvctrs_x_crndt_vctr_prdct = Q_small_egnvals(:).*x_crndt_vctr;

%the square of the metric
mtrc = x_crndt_vctr(:)'*Q_small_egnvctrs_x_crndt_vctr_prdct;

%the gradient before scaling by the inverse of the metric
mtrc_grad = Q_small_egnvctrs*Q_small_egnvctrs_x_crndt_vctr_prdct;

if ~isempty(Q_large_egnval)
    
    x_perp = x(:) - Q_small_egnvctrs*x_crndt_vctr;
    
    mtrc = mtrc + Q_large_egnval*sum(x_perp.^2);
    mtrc_grad = mtrc_grad + Q_large_egnval*x_perp;
    
end

mtrc = sqrt(mtrc);
mtrc_grad = (1/mtrc)*mtrc_grad;


% [QQ RR] = qr(P_dmnt_egnvctrs);
% %yy = QQ(:, (numel(P_dmnt_egnvals)+1):end)'*x;
% %objctv_grad = 2*[P_dmnt_egn_vctr_prdct(:).' yy(:).'];
% 
% DD = zeros(5);
% for k=1:3
%     DD(k,k) = P_dmnt_egnvals(k) + bias;
% end
% 
% for k=4:5
%    DD(k,k) = bias; 
% end
% 
% PP = ...
%     [P_dmnt_egnvctrs QQ(:, (numel(P_dmnt_egnvals)+1):end)] ...
%     *DD ...
%     *[P_dmnt_egnvctrs QQ(:, (numel(P_dmnt_egnvals)+1):end)]'
% 
% Q = inv(PP);

end
