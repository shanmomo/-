## FRESCO计算结果

f = open("FRESCO/He_Zr/201.txt","r" )
f_lines = readlines(f)
θf = [parse(Float64,l[1:10]) for l in f_lines]
σf = [parse(Float64,l[11:end]) for l in f_lines]
close(f)
