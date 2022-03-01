## FRESCO计算结果
f = open("FRESCO/He_Zr/201.txt","r" )
n = countlines( f )
seekstart( f )
θf,σf = zeros(n),zeros(n)
for i = 1:n
    x,y = split( readline( f ), "," )
    θf[i] = parse(Float64, x )
    σf[i] = parse(Float64, y )
end
close( f )
