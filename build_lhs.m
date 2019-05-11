function lhsmat = build_lhs(xs,ys)

np = length(xs) - 1;
psip = zeros(np,np+1);
lhsmat = zeros(np+1,np+1);


for i=1:np
    [infa, infb] = panelinf(xs(1),ys(1),xs(1+1),ys(1+1),xs(i),ys(i)); 
    psip(i,1)=infa;
    [infa, infb] = panelinf(xs(np),ys(np),xs(np+1),ys(np+1),xs(i),ys(i));
    psip(i,np+1)=infb;
     
     for j=2:np
        
        [infaP0, infbP0] = panelinf(xs(j),ys(j),xs(j+1),ys(j+1),xs(i),ys(i));
        [infaP1, infbP1] = panelinf(xs(j-1),ys(j-1),xs(j),ys(j),xs(i),ys(i));
        psip(i,j) = infaP0 - infbP1;
        
     end
end

lhsmat(1,:)=0;
lhsmat(1,1)=1;

lhsmat(np+1,:)=0;
lhsmat(np+1,np+1)=1;

for i=2:np
    for j=1:np+1
        if i == np
            lhsmat(i,j) = psip(1,j)-psip(i,j);
        else
  
            lhsmat(i,j) = psip(i+1,j)-psip(i,j);
        end
    end
end
