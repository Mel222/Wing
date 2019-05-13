function rhsvec = build_rhs(xs,ys,alpha)

np = length(xs) - 1;
rhsvec = zeros(np+1,1);

psiFS = transpose(ys)*cos(alpha)-transpose(xs)*sin(alpha);

rhsvec(np) = 0;
rhsvec(np+1) = 0;

for i=1:np-1
    
    rhsvec(i) = psiFS(i)-psiFS(i+1);
    
end


