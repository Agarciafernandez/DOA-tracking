function Y=Rotation_Mu(X,Mu)
% ROTATION_MU rotates the vector of points so that it coincides with the
%             North Pole (a.k.a. the standard z-axis)

%% rotation (orient mean vector (mu) with the North Pole)
Np=[0,0,1]; % z-axis (North Pole)
Mu=Mu/norm(Mu); % should be unit vector
    
if Mu(3)~= 1
    Ux=cross(Np,Mu);%axis of rotation
    Ux= Ux/norm(Ux);
    thetaX=acos(Mu(3));
    Rg= rotationVectorToMatrix(Ux*thetaX);
    Y=X*Rg;
end
