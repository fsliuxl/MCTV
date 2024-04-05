function x= myPCG_rpca(x,rhs,sizeD)
  temp =rhs(:);
  [x, ~] = pcg(@(x) Fun(x), temp, 1e-4,1000,[],[],x);   
    function y = Fun(x)
        y = diff_x(diff_x(x,sizeD),sizeD) + diff_y(diff_y(x,sizeD),sizeD) + x;
    end
end