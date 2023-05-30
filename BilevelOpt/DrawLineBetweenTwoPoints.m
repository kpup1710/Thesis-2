function DrawLineBetweenTwoPoints( point1, point2 )
p = length(point1);

if p == 2
    plot([point1(1),point2(1)],[point1(2),point2(2)], 'black');
    hold on
end

if p == 3
    plot3([point1(1),point2(1)],[point1(2),point2(2)],[point1(3),point2(3)],'black-');
%     plot3(point2(1),point2(2),point2(3),'m*','markersize',3);
    hold on
end

end

