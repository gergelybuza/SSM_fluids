function y = gausslobatto(obj,y_st,y_end)
% grid according to gauss lobatto points

grid = zeros(1,obj.N);
for i = 1:obj.N  
    grid(i) = cos((i-1)*pi/(obj.N-1));
end   
% map grid linearly to [y_st,y_end]
y = (y_end - y_st)*grid/2 + (y_end + y_st)/2;
y = y';
end