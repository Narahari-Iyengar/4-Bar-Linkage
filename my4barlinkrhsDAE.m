function zdot=my4barlinkrhsDAE(t,z,p)
x1=z(1); y1=z(2); th1=z(3);
x2=z(4); y2=z(5); th2=z(6);
x3=z(7); y3=z(8); th3=z(9);
x4=z(10); y4=z(11); th4=z(12);
x1d=z(13); y1d=z(14); th1d=z(15);
x2d=z(16); y2d=z(17); th2d=z(18);
x3d=z(19); y3d=z(20); th3d=z(21);
x4d=z(22); y4d=z(23); th4d=z(24);

names = fieldnames(p);
for i=1:length(names)
    eval([names{i} '= p.' names{i} ';' ]); % e.g.  a = p.a;
end

dots=[x1d; y1d; th1d; x2d; y2d; th2d; x3d; y3d; th3d; x4d; y4d; th4d];

M = mymassmatrixDAE(I1,I2,I3,I4,d1,d2,d3,d4,l1,l2,l3,l4,m1,m2,m3,m4,th1,th2,th3,th4);
b = myrhsstuffDAE(d1,d2,d3,d4,g,l1,l2,l3,l4,m1,m2,m3,m4,th1,th2,th3,th4,th1d,th2d,th3d,th4d);

a=M\b;
n=[a(1); a(2); a(3); a(4); a(5); a(6); a(7); a(8); a(9); a(10); a(11); a(12)];

zdot=[dots; n];
end