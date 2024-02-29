function latex_table(its::Matrix; cap=Inf, brackets=[])
    nr, nc = size(its)
    text = ""
    for i = 1:nr
        for j = 1:nc
            val = j > 1 ? Int(its[i,j]) : its[i,j]
            if j>1 && val==cap
                text = text * "-"
            else
                text = text * "$val"
                if !isempty(brackets) && j>1
                    text = text * " ($(brackets[i,j-1]))"
                end
            end
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
