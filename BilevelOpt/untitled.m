d1 = 600;
a1 = 450;
a2 = 350;

for phi1=-11*pi/12:0.2:11*pi/12
    for phi2= -49*pi/60:0.2:49*pi/60
        for d3 = 0:2:400
            x = a1*cos(phi1)+a2*cos(phi1+phi2);
            y = a1*sin(phi1)+a2*sin(phi1+phi2);
            z = d1-d3;
            plot3(x,y,z,'*')
            
            title('Khong gian lam viec cua Robot 3D');
            hold on
         end
     end
end