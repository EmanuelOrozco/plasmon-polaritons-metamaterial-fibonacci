%Programa para calcular la estructura de bandas de un cristal
%Fot?nico 2D por medio del m?todo de ondas planas PWE
%Referencia: Photonics Crystals, Physics and practical modeling
%Igor A. Sukhoivanov and Igor V. Guryev
%Por: Hern?n Alejandro G?mez Urrea


%Input parameters: PhC period, radius of an
%element, permittivities
%Output data: The band structure of the 2D PhC.
clear all;

%Variables l1 and l2 contain thicknesses of each
%layer within unit cell

global numG
global kxx
global kyy
global Vs
global in1
global G

%The variable a defines the period of the structure.
%It influences on results in case of the frequency
%normalization.
a=1e-6; % in meters
%The variable r contains elements radius. Here it is
%defined as a part of the period.
r=0.4*a; % in meters
%The variable eps1 contains the information about
%the relative background material permittivity.
eps1=9;
%The variable eps2 contains the information about
%permittivity of the elements composing PhC. In our
%case permittivities are set to form the structure
%perforated membrane
eps2=1;
%The variable precis defines the number of k-vector
%points between high symmetry points
precis=5;%5
%The variable nG defines the number of plane waves.
%Since the number of plane waves cannot be
%arbitrary value and should be zero-symmetric, the
%total number of plane waves may be determined
%as (nG*2-1)^2
nG=4; %4
%The variable precisStruct defines contains the
%number of nodes of discretization mesh inside in
%unit cell discretization mesh elements.
precisStruct=150;
%The following loop carries out the definition
%of the unit cell. The definition is being
%made in discreet form by setting values of inversed

%dielectric function to mesh nodes.
nx=1;
for countX=-a/2:a/precisStruct:a/2
ny=1;
    for countY=-a/2:a/precisStruct:a/2
%The following condition allows to define the circle
%with of radius r
        if(sqrt(countX^2+countY^2)<r)
%Setting the value of the inversed dielectric
%function to the mesh node
            struct(nx,ny)=1/eps2;
%Saving the node coordinate
            xSet(nx)=countX;
            ySet(ny)=countY;
        else
            struct(nx,ny)=1/eps1;
            xSet(nx)=countX;
            ySet(ny)=countY;
        end
        ny=ny+1; 
    end
    nx=nx+1; 
end

imagesc(struct); %Para visualizar la estructura

%The computation of the area of the mesh cell. It is
%necessary for further computation of the Fourier
%expansion coefficients.
dS=(a/precisStruct)^2;%Lo estamos tomando sobre la divisi?n mas
%peque?a
%Forming 2D arrays of nods coordinates
xMesh=meshgrid(xSet(1:length(xSet)-1));%note que este es un vector fila
yMesh=meshgrid(ySet(1:length(ySet)-1))';%note que este es un vector columna
%En Matlab la coma indica esto
%Transforming values of inversed dielectric function
%for the convenience of further computation
structMesh=struct(1:length(xSet)-1,...
       1:length(ySet)-1)*dS/(max(xSet)-min(xSet))^2;%Esto para la 
   %integracion
%Defining the k-path within the Brilluoin zone. The
%k-path includes all high symmetry points and precis
%points between them
%Se parametriza teniendo en cuenta precis, est? me dar? el n?mero 
%de puntos a considerar
%Note que se construye kx y ky de tal manera que contengan todo el conjunto

%GX En este ky=0 para todo kx entre 0 y pi/a

kx(1:precis+1)=0:pi/a/precis:pi/a;
ky(1:precis+1)=zeros(1,precis+1);
%XM
%En este camino kx=pi/a para todo  ky entre 0 y pi/a
kx(precis+2:precis+precis+1)=pi/a;
ky(precis+2:precis+precis+1)=...
                       pi/a/precis:pi/a/precis:pi/a;
%MG 
%En este tenemos una l?nea recta que se recorre desde (pi/a,pi/a) hasta 
% (0,0), por eso el menos en el contador, para este caso kx=ky
kx(precis+2+precis:precis+precis+1+precis)=...
                    pi/a-pi/a/precis:-pi/a/precis:0;
ky(precis+2+precis:precis+precis+1+precis)=...
                    pi/a-pi/a/precis:-pi/a/precis:0;
%After the following loop, the variable numG will
%contain real number of plane waves used in the
%Fourier expansion
numG=1;
%The following loop forms the set of reciprocal
%lattice vectors.
for Gx=-nG*2*pi/a:2*pi/a:nG*2*pi/a
    for Gy=-nG*2*pi/a:2*pi/a:nG*2*pi/a
        G(numG,1)=Gx;
        G(numG,2)=Gy;
        numG=numG+1;
    end
end

%The next loop computes the Fourier expansion
%coefficients which will be used for matrix
%differential operator computation.
for countG=1:numG-1
  for countG1=1:numG-1
     CN2D_N(countG,countG1)=sum(sum(structMesh.*...
     exp(1i*((G(countG,1)-G(countG1,1))*...
     xMesh+(G(countG,2)-G(countG1,2))*yMesh))));
  end
end
%The next loop computes matrix differential operator
%in case of TE mode. The computation
%is carried out for each of earlier defined
%wave vectors.
for countG=1:numG-1
 for countG1=1:numG-1
  for countK=1:length(kx)%este cuenta sobre los vectores de en la ZIB
   M(countK,countG,countG1)=...
CN2D_N(countG,countG1)*((kx(countK)+G(countG,1))*...
                       (kx(countK)+G(countG1,1))+...
(ky(countK)+G(countG,2))*(ky(countK)+G(countG1,2)));

  end
 end
end
%The computation of eigen-states is also carried
%out for all wave vectors in the k-path.
for countK=1:length(kx)
%Taking the matrix differential operator for current
%wave vector.
  MM(:,:)=M(countK,:,:);
%Computing the eigen-vectors and eigen-states of
%the matrix
  [D V]=eig(MM);
%Transforming matrix eigen-states to the form of
%normalized frequency.
  dispe(:,countK)=sqrt(V*ones(length(V),1))*a/2/pi;
end


%Plotting the band structure
%Creating the output field
figure(3);
%Creating axes for the band structure output.
ax1=axes;
%Setting the option
%drawing without cleaning the plot
hold on;
%Plotting the first 8 bands
for u=1:8
  plot(abs(dispe(u,:)),'r','LineWidth',2);
%If there are the PBG, mark it with blue rectangle
  if(min(dispe(u+1,:))>max(dispe(u,:)))
    rectangle('Position',[1,max(dispe(u,:)),...
    length(kx)-1,min(dispe(u+1,:))-...
               max(dispe(u,:))],'FaceColor','b',...
                                'EdgeColor','b');
  end
end
%Signing Labeling the axes
set(ax1,'xtick',...
              [1 precis+1 2*precis+1 3*precis+1]);
set(ax1,'xticklabel',['G';'X';'M';'G']);
ylabel('\omegaa/2\pic','FontSize',16);
%xlabel('Wavevector','FontSize',16);
xlim([1 16])
set(ax1,'XGrid','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Hasta aqu? relaci?n de dispersi?n
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
   
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ahora vamos a considerar el perfil del campo para alg?n valor de k,
%en la ZIB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kxx=0*pi/(a);%cambiar este dependiendo del valor de k que quiera 
%considerar
%Falta definir la frecuencia, este valor me va a dar todos los
%autovalores de la matriz para el respectivo k, decidir cual voy a
%considerar
kyy=0*pi/(a)
 
%The next loop computes matrix differential operator
%in case of TE mode. The computation
%is carried out for each of earlier defined
%wave vectors.
%Calculo la matriz para el valor de k=(kxx,kyy) escogido
for countG=1:numG-1
 for countG1=1:numG-1
   MM1(countG,countG1)=...
CN2D_N(countG,countG1)*((kxx+G(countG,1))*...
                       (kxx+G(countG1,1))+...
(kyy+G(countG,2))*(kyy+G(countG1,2)));

  end
end

%Computing the eigen-vectors and eigen-states of
%the matrix
[V,D]=eig(MM1);
%Transforming matrix eigen-states to the form of
%normalized frequency.

%Here we sort eigen-values obtained and normalize
%them
[d,ind] = sort(diag(D));%Da los autovalores ordenados
Ds = D(ind,ind);%Matriz diagonal con los autovalores ordenados
Vs = V(:,ind);%Eigenvectores ordenados
%ind representa el numero de la banda a la que le quiero calcular la
%frecuencia para el respectivo k

in1=4;%n?mero de banda
freq=sqrt(abs(Ds(in1,in1)))*a/2/pi;

figure(2);
%Disabling window cleaning before plot of next graph
%hold on;

%Norma al cuadrado de Hz (Hz*Hz^*)
FieldHz=Field(xMesh,yMesh).*conj(Field(xMesh,yMesh));

imagesc(FieldHz); %Para visualizar el campo
%Funciones
%Calcular la dependencia del campo con respecto a x.
%Se toman los autovalores para cada banda y se desarrolla en la base de
%ondas planas



function H = Field(x,y)
global numG
global kxx
global kyy
global Vs
global in1
global G

H=0;
countt=0;
    for countG=1:numG-1
     countt=countt+1;
     H = H+Vs(countt,in1)*exp(1i*((kxx+G(countG,1)*x)+ ...
         (kyy+G(countG,2)*y)));
     end
end



   
   