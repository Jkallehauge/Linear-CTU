function Ainv=Analytical_inversion3X3(ATA)

% Analytically inverts a 3X3 matrix

A=ATA(2,2)*ATA(3,3)-ATA(3,2)*ATA(2,3);
B=-(ATA(2,1)*ATA(3,3)-ATA(2,3)*ATA(3,1));
C=ATA(2,1)*ATA(3,2)-ATA(2,2)*ATA(3,1);

D=-(ATA(1,2)*ATA(3,3)-ATA(1,3)*ATA(3,2));
E=ATA(1,1)*ATA(3,3)-ATA(1,3)*ATA(3,1);
F=-(ATA(1,1)*ATA(3,2)-ATA(1,2)*ATA(3,1));

G=ATA(1,2)*ATA(2,3)-ATA(1,3)*ATA(2,2);
H=-(ATA(1,1)*ATA(2,3)-ATA(1,3)*ATA(2,1));
I=ATA(1,1)*ATA(2,2)-ATA(1,2)*ATA(2,1);

Ainv=(1/(ATA(1,1)*A+ATA(1,2)*B+ATA(1,3)*C))*[A D G; B E H; C F I];
