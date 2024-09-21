function singvals = optshrink_impl(singvals,beta,loss)

%y = @(x)( (1+sqrt(beta)).*(x<=beta^0.25) + sqrt((x+1./x) ...
%     .* (x+beta./x)).*(x>(beta^0.25)) );
%     assert(sigma>0)
%     assert(prod(size(sigma))==1)

x = @(y)( sqrt(0.5*((y.^2-beta-1 )+sqrt((y.^2-beta-1).^2 - 4*beta) ))...
    .* (y>=1+sqrt(beta)));

% this is found to be not exactly right for all y's < 1+sqrt(beta).
%opt_fro_shrink = @(y)( sqrt(max(((y.^2-beta-1).^2 - 4*beta),0) ) ./ y);
opt_fro_shrink = @(y)( sqrt((y.^2-beta-1).^2 - 4*beta) ./ y);

opt_op_shrink = @(y)(max(x(y),0));
opt_nuc_shrink = @(y)(max(0, (x(y).^4 - sqrt(beta)*x(y).*y - beta)) ...
    ./((x(y).^2) .* y));

switch loss
    case 'fro'
        %singvals = sigma * opt_fro_shrink(singvals/sigma);
        % to fix the bug:
        singvals1 = opt_fro_shrink(singvals);
        singvals1(singvals<1+sqrt(beta))=0;
        singvals= singvals1;
    case 'nuc'
        %         y = singvals/sigma;
        %         singvals = sigma * opt_nuc_shrink(y);
        singvals = opt_nuc_shrink(y);
        singvals((x(y).^4 - sqrt(beta)*x(y).*y - beta)<=0)=0;
    case 'op'
        %         singvals = sigma * opt_op_shrink(singvals/sigma);
        singvals = opt_op_shrink(singvals);
    otherwise
        error('loss unknown')
end

end
