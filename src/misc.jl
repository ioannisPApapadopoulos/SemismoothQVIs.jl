function latex_table(its::Matrix)
    nr, nc = size(its)
    text = ""
    for i = 1:nr
        for j = 1:nc
            val = j > 1 ? Int(its[i,j]) : its[i,j]
            text = text * "$val"
            # val = its[i,j]
            if j == nc 
                text = text * "\\\\\n"
            else
                text = text * " &"
            end
        end
    end
    return text
end
