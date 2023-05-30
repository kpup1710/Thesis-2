function arr_new_vk = GetNewProperElements(vk, wk, p, m, M)
format long
    e = eye(p);
    arr_new_vk = [];
    for i = 1 : p
        v_k_i = vk-e(:,i).*(vk - wk);
        if ~IsInsideBox(v_k_i, m, M) %any(i == improperIndices) ||
            continue;
        end
        arr_new_vk = [arr_new_vk ;v_k_i'];
    end
    
end

