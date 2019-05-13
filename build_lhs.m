function lhsmat = build_lhs2(xs,ys)

np = length(xs) - 1;
psip = zeros(np,np+1);
lhsmat = zeros(np+1,np+1);


for i=1:np
    for j=1:np+1
        if j == 1
            
            [infa, infb] = panelinf(xs(j),ys(j),xs(j+1),ys(j+1),xs(i),ys(i)); 
            psip(i,j)=infa;
            
        elseif j == np+1
            
            [infa, infb] = panelinf(xs(j-1),ys(j-1),xs(j),ys(j),xs(i),ys(i));
            psip(i,j)=infb;
            
        else     
        
            [infaP0, infbP0] = panelinf(xs(j),ys(j),xs(j+1),ys(j+1),xs(i),ys(i));
            [infaP1, infbP1] = panelinf(xs(j-1),ys(j-1),xs(j),ys(j),xs(i),ys(i));
            psip(i,j) = infaP0 + infbP1;
        end 
     end
end


lhsmat(np+1,:)=0;
lhsmat(np+1,np+1)=1;

lhsmat(np,:)=0;
lhsmat(np,1)=1;

for i=1:np-1
    for j=1:np+1
  
            lhsmat(i,j) = psip(i+1,j)-psip(i,j);
    end
end
