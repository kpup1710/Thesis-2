function check = IsInsideBox(v, a, b)
   % check if v in [a, b]
    check = false;
    if IsSmallerThan(a, v) && IsSmallerThan(v, b)
        check = true;
    end
end
