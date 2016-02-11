clear all;
close all;

%% Initialization
[points, seg, tri, edge] = importfilemesh('mesh50.msh');
nNodes = size(points,2); nElems = size(tri,2); nnNodes = 2*nNodes;

% Steel material properties
E = 200E09; % Young's modulus
nu = 0.30; % Poisson's ratio
rho = 7750; % [kg/m^3] density
t = 1; % [m] thickness

g = -1000; 
P = 1E04;

Penalty = 10E15;

%% Boundary conditions

ibcd = edge(edge(:,3)==1,1:2)';

isBdNode = false(nNodes,1);
isBdNode(ibcd) = true;
ibcd = find(isBdNode);

ibcdneum = edge(edge(:,3)~=1,1:2)';

nnb=2*length(ibcd);
ibc=zeros(nnb,1); ibc(1:2:end)=2*ibcd-1; ibc(2:2:end)=2*ibcd;

%% Two-dimensional linear elasticity assembly of the stiffness matrix 
% That is , we will solve K = assemble(points,tri,E,nu)

% Constant matrices 
D = (E/((1+nu)*(1-(2*nu))))*[1-nu,nu,0;nu,1-nu,0;0,0,0.5*(1-(2*nu))];
Zh = [zeros(1,2);eye(2)];

% Loop over all elements
BData = zeros(3,6,nElems);
K = sparse(nnNodes,nnNodes);

for currElem = 1:nElems
	currNodes = tri(1:3,currElem); 
	currIndex = 2*tri([1,1,2,2,3,3],currElem)-[1;0;1;0;1;0]; 
	
	xy = zeros(3,2); 
	xy(:,1) = points(1,currNodes)'; 
	xy(:,2) = points(2,currNodes)'; 
	
	% T is the triangle matrix 
	% where det(T) = Area(triangle)
	T = [1,1,1;xy'];
	dN = T\Zh;
	
	B = zeros(3,6);
	B([1,3],[1,3,5]) = dN'; 
	B([3,2],[2,4,6]) = dN'; 
	
	BData(:,:,currElem) = B;
	
	K(currIndex,currIndex) = K(currIndex,currIndex)+0.5*det(T)*B'*D*B;
end


%% Assemble body forces (in this case just gravity)

f1 = zeros(nNodes,1); f2 = t*g*rho*ones(nNodes,1);

% triangle vertices
it1 = tri(1,:); it2=tri(2,:); it3=tri(3,:);

% edge vectors
a21 = points(:,it2)-points(:,it1); a31=points(:,it3)-points(:,it1);

% area of triangles
area = abs(a21(1,:).*a31(2,:)-a21(2,:).*a31(1,:))/2;

% assembly
f1h=(f1(it1)+f1(it2)+f1(it3)).*area'/9; 
f2h=(f2(it1)+f2(it2)+f2(it3)).*area'/9; 
ff1=sparse(it1,1,f1h,nNodes,1)+sparse(it2,1,f1h,nNodes,1)+sparse(it3,1,f1h,nNodes,1); 
ff2=sparse(it2,1,f2h,nNodes,1)+sparse(it2,1,f2h,nNodes,1)+sparse(it3,1,f2h,nNodes,1); 

% right-hand side
f=zeros(nnNodes,1);
f(1:2:nnNodes)=full(ff1);
f(2:2:nnNodes)=full(ff2);

FGrav = f;


%% Pressure force vector
FPress=zeros(nnNodes,1);
for i=1:length(edge)
    node1_temp=edge(i,1);
    node2_temp=edge(i,2);
    border=edge(i,3);
    
    x1=points(1,node1_temp);
    x2=points(1,node2_temp);
    
    %reorder the node
    if x2<x1
        node1=node2_temp;
        node2=node1_temp;
    else
        node1=node1_temp;
        node2=node2_temp;
    end
    
    x1=points(1,node1);
    x2=points(1,node2);
    y1=points(2,node1);
    y2=points(2,node2);
    L=((x1-x2)^2+(y1-y2)^2)^(1/2);
   
    if border~=1
    if border==5  
        tx=(y2-y1)/L*P*t;
        ty=-(x2-x1)/L*P*t;
    elseif border==4 
        tx=-(y2-y1)/L*P*t;
        ty=(x2-x1)/L*P*t;
    elseif border==3 
        tx=(y2-y1)/L*P*t;
        ty=-(x2-x1)/L*P*t;
    elseif border==2
        tx=-(y2-y1)/L*P*t;
        ty=(x2-x1)/L*P*t;
    end
    
    FPress(2*(node1-1)+1,1)=FPress(2*(node1-1)+1,1)+L/2*tx;
    FPress(2*(node1-1)+2,1)=FPress(2*(node1-1)+2,1)+L/2*ty;
    FPress(2*(node2-1)+1,1)=FPress(2*(node2-1)+1,1)+L/2*tx;
    FPress(2*(node2-1)+2,1)=FPress(2*(node2-1)+2,1)+L/2*ty;
    end    
    
    
end

%% Reduce system order, put together forces, and solve the linear system

F = FPress+FGrav;

K(ibc,ibc)=K(ibc,ibc)+Penalty*speye(nnb);

F(ibc)=0;

U = K\F;

%% Calculate stresses
vonMises = zeros(nElems,1);
for currElem = 1:nElems
	currNodes = tri(1:3,currElem); 
	currIndex = 2*tri([1,1,2,2,3,3],currElem)-[1;0;1;0;1;0]; 
	
	Be = BData(:,:,currElem);
	
	sigma = D*Be*U(currIndex);
	
	vonMises(currElem) = sqrt( sigma(1)^2-sigma(1)*sigma(2)+sigma(2)^2+3*sigma(3)^2);
end

disp(nElems)

[mavVM, locMaxVM] = max(vonMises);

disp(mavVM);



[Sh,Vms]=elas2dsvms(points,tri,U,E,nu);


%% Visualization

scale = 100;

uu=reshape(U,2,nNodes)';

[minDisp, locMinDisp] = min(uu(:,2));
disp(minDisp);

pu=points(1:2,:)'+scale*uu;
colormap('jet') 
trisurf(tri(1:3,:)',pu(:,1),pu(:,2),zeros(nNodes,1),Vms,'facecolor','interp','EdgeColor','none') 
view(2)
axis equal
hold on

% patch(pu(tri(1:3,locMaxVM),1),pu(tri(1:3,locMaxVM),2),[1 0 0],'LineWidth',2)
plot(pu(locMinDisp,1),pu(locMinDisp,2),'ro','LineWidth',2,'MarkerSize',10);
plot(mean(pu(tri(1:3,locMaxVM),1)),mean(pu(tri(1:3,locMaxVM),2)),'-mo',...
                'LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','r',...
                'MarkerSize',10)
