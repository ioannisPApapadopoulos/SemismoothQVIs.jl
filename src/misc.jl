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

function test2_latex_table(its, eoc, converge)
    nr, nc = size(its)
    text = ""

    text = text * "{" * repeat("c |", nc) *"c} \n\n"
    text = text * "R"
    for j in 1:nc
        val = its[1,j]
        text = text * "& $val"
    end
    text = text * "\\\\ \n"

    for i = 2:nr
        for j = 1:nc
            text * ""
            if !converge[i-1,j]
                text = text * "& \\begin{tabular}{@{}c@{}} -  \\end{tabular}   \n"
            else
                text = text * "& \\begin{tabular}{@{}c@{}} $(Int(its[i,j])) \\\\ ($(eoc[i-1,j]))  \\end{tabular}   \n"
            end
            if j == nc
                text = text *"\n\n"
            end
        end
    end
    return text
end

function test4_latex_table(its::Matrix, its2::Matrix)
    nr, nc = size(its)
    text = ""
    for i = 1:nr
        text = text * "$(its[i,1])  &"
        for j = 2:nc
            val = Int(its[i,j])
            val2 = Int(its2[i,j])
            text = text * "$val ($val2)"
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