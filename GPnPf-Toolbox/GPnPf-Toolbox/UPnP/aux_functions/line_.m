function line_(point1,point2,col,w)

if size(point1,1)*size(point1,2)==3
    line([point1(1),point2(1)],[point1(2),point2(2)],[point1(3),point2(3)],'color',col,'linewidth',w);
else
    line([point1(1),point2(1)],[point1(2),point2(2)],'color',col,'linewidth',w);
end

