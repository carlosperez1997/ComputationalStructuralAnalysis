function plotBeams3D(Tnod,nsub,le,uint,e_plot)
% PLOTBEAMS3D - Plot displacements, rotations for 3D beam in local reference fram
% Inputs:
%   Tnod     Nodal connectivities matrix [Nelements x NnodesXelement]
%   nsub     Number of elements's subdivisions to evaluate displacements and rotations
%	le		 vector with element lengths [Nelements x 1]
%	uint	 matrix with element displacements and rotations in LOCAL axes [12 x Nelements] 
%	e_plot	 Element to plot


% Compute deflection and rotation polynomial coefficients
nel = size(Tnod,1);
puy = zeros(nel,4);
pty = zeros(nel,3);
puz = zeros(nel,4);
ptz = zeros(nel,3);

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
	ptz(e,:) = [3*coeff(1),2*coeff(2),coeff(3)];
	puz(e,:) = [coeff(5),coeff(6),coeff(7),coeff(8)];
	pty(e,:) = [3*coeff(5),2*coeff(6),coeff(7)];
end

x_coord = linspace(0,le(e_plot),nsub);
y_coord = polyval(puy(e_plot,:),x_coord);
z_coord = polyval(puz(e_plot,:),x_coord);
theta_y = polyval(pty(e_plot,:),x_coord);
theta_z = polyval(ptz(e_plot,:),x_coord);

fig = figure();
set(fig,'Name',['Deformation of element ',num2str(e_plot)])

% Plot beam deflections
subplot(2,2,1)
plot(x_coord,y_coord);
set(gca,'XTickLabel',sprintf('%i\n',Tnod(e_plot,:)),'XTick',[0 le(e_plot)]);
xlabel('beam');
ylabel('U_y (m)');
title('Deflection in y');
grid
grid minor
xlim(gca,[0 le(e_plot)]);

subplot(2,2,2)
plot(x_coord,z_coord);
set(gca,'XTickLabel',sprintf('%i\n',Tnod(e_plot,:)),'XTick',[0 le(e_plot)]);
xlabel('beam');
ylabel('U_z (m)');
title('Deflection in z');
grid
grid minor
xlim(gca,[0 le(e_plot)]);

% Plot beam section rotations
subplot(2,2,4)
plot(x_coord,theta_y);
set(gca,'XTickLabel',sprintf('%i\n',Tnod(e_plot,:)),'XTick',[0 le(e_plot)]);
xlabel('beam');
ylabel('\theta_y (rad)');
title('Rotation in y');
grid
grid minor
xlim(gca,[0 le(e_plot)]);

subplot(2,2,3)
plot(x_coord,theta_z);
set(gca,'XTickLabel',sprintf('%i\n',Tnod(e_plot,:)),'XTick',[0 le(e_plot)]);
xlabel('beam');
ylabel('\theta_z (rad)');
title('Rotation in z');
grid
grid minor
xlim(gca,[0 le(e_plot)]);

end