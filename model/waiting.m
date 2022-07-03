function waiting(i,j)
    
    if i == 1
    fprintf('( >)>')
    elseif i == j
    fprintf('\b\b\b\b\b\n')
    elseif ~rem(i,2)==0
    fprintf('\b\b\b\b\b( >)>')
    elseif rem(i,2)==0
    fprintf('\b\b\b\b\b(< )>')
    end

end
