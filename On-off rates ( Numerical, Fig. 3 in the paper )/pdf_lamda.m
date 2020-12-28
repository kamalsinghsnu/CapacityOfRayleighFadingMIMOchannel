function pdf = pdf_lamda(y)



global m n;

% m = min(t,r)
% n = max(t,r) 



temp1 = (  exp(-y) .* (y.^(n-m))  ) ./ m; %outer term
        
%outer summation loop
        
Sum1 = 0;
        
for k=0:1:m-1,
           
            
            temp2 = factorial(k) ./ factorial(k + n-m);
      
            Sum2 = 0;
            
            for j = 0:1:k,
                
               
                temp3 = ((-y).^j) .*  nchoosek(k + n - m, k - j) ./ factorial(j);
                
                
                %temp3 =  factorial(n-m+k) ./ (factorial(k-j) .* factorial(n-m+j) .* factorial(j) )     ;
                
                Sum2 = Sum2  +  temp3;
                
                
            end
            
            Sum1 = Sum1 + (Sum2.^2) .* temp2;
        
end

pdf = temp1 .* Sum1;










