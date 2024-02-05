Matrix = [Parameters_list, Manual_Label']
[u,z] = pca(Parameters_list,'NumComponents',3); 
% Ploting PCA 3D
figure
scatter3(z(:,1),z(:,2),z(:,3),40,Matrix(:,end),'filled'); 



clear  NewMatrix Matrix
Matrix = [Parameters_list, Manual_Label']
Size = size(Matrix,1); 
Vector = 1:Size; 
radomVector =  Vector(randperm(length(Vector)));
suma = 1

for i = radomVector 
    NewMatrix(suma,:) = Matrix(i,:)
    suma = suma + 1

end 

imagesc(NewMatrix)

[u,z] = pca(NewMatrix,'NumComponents',3); 

% Ploting PCA 3D
figure
scatter3(z(:,1),z(:,2),z(:,3),40,NewMatrix(:,end),'filled'); 
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');