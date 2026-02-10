function z = getZstep(infotable,filename)
    rowIdx = strcmp({infotable.name}, filename);
    idx = find(rowIdx==1);
    z = infotable(idx).Var2;
end