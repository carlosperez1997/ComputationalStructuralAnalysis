function plotBeams3D_def(x,Tnod,nsub,le,uint,factor,Re)
% PLOTBEAMS3D_DEF - Plot deformed structure of 3D BEAMS
% Inputs:
%   x        Nodal coordinates matrix [Nnodes x Ndim]
%   Tnod     Nodal connectivities matrix [Nelements x NnodesXelement]
%   nsub     Number of elements's subdivisions to evaluate displacements and rotations
%	le		 vector with element lengths [Nelements x 1]
%	uint	 matrix with element displacements and rotations in LOCAL axes [12 x Nelements]
%	factor	 Number to scale/amplify the displacements
%	Re		 Rotation matrix [12 x 12 x Nelements]


fig1 = figure();
set(fig1,'Name','Deformed structure','Color','w');

% Initial plot
sz = 2;
ori = [x(1,1),x(1,2),x(1,3)]-0.1;
hold on;
plot3(ori(1)+[0,sz],ori(2)+[0,0],ori(3)+[0,0],'r','linewidth',1.5); text(ori(1)+1.1*sz,ori(2),ori(3),'X','color','r');
plot3(ori(1)+[0,0],ori(2)+[0,sz],ori(3)+[0,0],'g','linewidth',1.5); text(ori(1),ori(2)+1.1*sz,ori(3),'Y','color','g');
plot3(ori(1)+[0,0],ori(2)+[0,0],ori(3)+[0,sz],'b','linewidth',1.5); text(ori(1),ori(2),ori(3)+1.1*sz,'Z','color','b');
patch('Vertices',x,'Faces',Tnod,'edgecolor',[0.5,0.5,0.5],'linewidth',1);

% Compute deflection polynomial coefficients and plot the deformed beams
% using nsub subdivisions

nel = size(Tnod,1);
puy = zeros(nel,4);
puz = zeros(nel,4);

for e = 1:nel
	temp1 = [
		2,      le(e),   -2,    le(e);
		-3*le(e), -2*le(e)^2, 3*le(e), -le(e)^2;
		0,    le(e)^3,    0,     0;
		le(e)^3,       0,    0,     0];
	temp2 = temp1; temp2(:,[2,4]) = -temp2(:,[2,4]);
	
	matCoeff = zeros(8,12);
	matCoeff(1:4,[2,6,8,12]) = temp1;
	matCoeff(5:8,[3,5,9,11]) = temp2;
	
	coeff = 1/le(e)^3*matCoeff*uint(:,e);
	puy(e,:) = [coeff(1),coeff(2),coeff(3),coeff(4)];
	puz(e,:) = [coeff(5),coeff(6),coeff(7),coeff(8)];
		
	x_coord = linspace(0,le(e),nsub+1);
	uy_coord = polyval(puy(e,:),x_coord);
	uz_coord = polyval(puz(e,:),x_coord);
	ux_coord = linspace(uint(1,e),uint(7,e),nsub+1);
	
	temp_coord = x(Tnod(e,1),:)' +  Re(1:3,1:3,e)' * [x_coord;zeros(2,size(x_coord,2))]  + factor * Re(1:3,1:3,e)' * [ux_coord;uy_coord;uz_coord];
	patch('Vertices',temp_coord','Faces',[1:nsub;2:nsub+1]','FaceColor','none','LineWidth',2);
end

view(30,25);
axis equal;
set(gca,'xcolor','none','ycolor','none','zcolor','none');