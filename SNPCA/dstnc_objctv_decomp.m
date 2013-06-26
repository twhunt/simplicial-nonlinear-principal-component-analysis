function [objctv_val, objctv_grad] = dstnc_objctv_decomp(...
    x, P_dmnt_egnvctrs, P_dmnt_egnvals, bias)

%compute induced metric inv(P + epsilon*I) where P has only a few
%significant eigenvalues

%P = V*D*V'
%D = [D_d 0; 0 0] V=[Q1 Q2]
%inv(P + epsilon*I) = V*[inv(D_d + epsilon*I) 0; 0 (1/epsilon)*I]
%x'*inv(P+epsilon*I)*x = (Q1'*x)'*inv(D_d + epsilon*I)*(Q1*x)
%                        + (1/epsilon)*(Q2'*x)'*(Q2*x)

x_norm_sqrd = x(:).'*x(:);

bias_inv = 1/bias;
%P_dmnt_egnvals + bias is vector plus scalar
%(P_dmnt_egnvctrs + bias)) is diagonal of inverse of diagonal matrix

x_dmnt_egnvctrs_crndt_vctr = P_dmnt_egnvctrs'*x;
P_dmnt_egn_vctr_prdct = ...
    (1./(P_dmnt_egnvals + bias)).*x_dmnt_egnvctrs_crndt_vctr;

x_dmnt_egnvctrs_crndt_vctr_norm_sqrd = ...
    x_dmnt_egnvctrs_crndt_vctr'*x_dmnt_egnvctrs_crndt_vctr;

term1 = x_dmnt_egnvctrs_crndt_vctr'*P_dmnt_egn_vctr_prdct;
term2 = bias_inv*(x_norm_sqrd - x_dmnt_egnvctrs_crndt_vctr_norm_sqrd);

objctv_val  = (term1 + term2);

objctv_grad = P_dmnt_egnvctrs*P_dmnt_egn_vctr_prdct;
objctv_grad = objctv_grad ...
    + bias_inv*(x - P_dmnt_egnvctrs*x_dmnt_egnvctrs_crndt_vctr);
objctv_grad = 2*objctv_grad';


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
