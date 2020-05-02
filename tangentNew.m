function[vect] = tangentNew(img)
% function tangentNew calculates the tangent angles formed by traversing the entire contour

img=imgaussfilt(rgb2gray(img));
image=imcomplement(imbinarize(img));
image=bwskel(image);
%preprocessing steps
[row,col]=size(image);
[boundaries,label,nlabel,A]=bwboundaries(image,8,'noholes');
countPixels=zeros(1,nlabel);
%[L,nL]=bwlabel(image);
for i=1:row
    for j=1:col
        if label(i,j)>0
            countPixels(1,label(i,j))=countPixels(1,label(i,j))+1;
        end
    end
end

k=1;
m=sum(countPixels);
tangent=zeros(1,m);
for i=1:nlabel
    p=1;
    x=boundaries{i}(:,2);
    y=boundaries{i}(:,1);
    window=3; %the window size chosen for calculating the tangent angles

    while((p+window)<=numel(x))
        tangent(1,k)=atan((y(p+window,1)-y(p,1))/(x(p+window,1)-x(p,1)));
        k=k+1;
        p=p+1;
    end
end
%tangent is the required feature vector
y=sort(tangent);
q1=median(y(find(y<median(y))));
q2=median(y);
q3=median(y(find(y>median(y))));
vect=[max(tangent) min(tangent) mean(tangent) std(tangent) q1 q2 q3] %vect= vector containing the statistical features
end