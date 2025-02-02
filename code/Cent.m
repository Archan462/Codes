function [vect] = Cent(im)
%function 'cent' computes the centoid distance vector of a given contour
im=imcomplement(imbinarize(rgb2gray(im)));
count=0;


g=bwskel(im);
for i=1:size(g,1)
    for j=1:size(g,2)
        if g(i,j)==1
            count=count+1;
        end
    end
end
sx=0;
sy=0;
for i=1:size(g,1)
    for j=1:size(g,2)
        if g(i,j)==1
          sx=sx+i;
          sy=sy+j;
        end
    end
end
xcent=sx/(count+1); %xcent= x-coordinate of centroid
ycent=sy/(count+1); %ycent= y-coordinate of centroid 

CentDist=zeros(count,1);
count=0;
for i=1:size(g,1)
    for j=1:size(g,2)
        if g(i,j)==1
            d=sqrt((i-xcent)^2+(j-ycent)^2);
            CentDist(count+1)=d;
            count=count+1;
        end
    end
end
A=sum(CentDist);
CentDist=CentDist/A; %Centdist= Feature vector of centroid distance
y=sort(CentDist);
q1=median(y(find(y<median(y))));
q2=median(y);
q3=median(y(find(y>median(y))));
vect=[max(CentDist) min(CentDist) mean(CentDist) std(CentDist) q1 q2 q3] %vect= vector containing the 7 statistical parameters of the feature vector

end



