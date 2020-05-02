function feature=main_chaincode(im)
s=strel('disk',4,0);
g=imclose(~im, s);
B = bwboundaries(g,8);
hB=size(B,1);
feature=zeros(1,8);
for sim=1:hB
    b=B{sim}; 
    c = fchcode2(b,8);
    for i = 1:length(c.fcc)
        feature(1,1+c.fcc(i))=feature(1,1+c.fcc(i))+1;        
    end
end
end