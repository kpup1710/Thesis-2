function isSmaller = IsSmallerThan( a, b)
format long;
% a < b <=> a(i) < b(i) \forall i = 1, p
% number of dimensions
p = length(a);

for i = 1 : p
    if a(i) - b(i) > 0
       isSmaller = false;
       return;
    end
end

isSmaller = true;

end

