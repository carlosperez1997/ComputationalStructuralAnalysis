clear all
close all
clc

install_UNVARTOP

caso = 3;

%% TEST EXAMPLE
if caso == 0
    F_msg = 'F(n_unkn*find(coord(:,2)==0 & coord(:,1)==nelx),1) = -0.01*nelx;';
    fixed_dofs_msg = 'fixed_dofs = reshape(n_unkn*find(coord(:,1)==0)+(-n_unkn+1:0),1,[]);';
    active_msg = 'active_node = [];';
    passive_msg = 'passive_node = [];';

    %[iter,J, coord, U, connect] = ...
    % UNVARTOP_2D_compliance (nelx,nely,nsteps,Vol0,Vol,k,tau, 
    %F_msg, fixed_dofs_msg, active_msg, passive_msg)

    [iter,J, coord, U, connect] = UNVARTOP_2D_compliance_modified(100, 50, 10, 0, 0.5, 0, 0.5, ...
        F_msg, fixed_dofs_msg, active_msg, passive_msg);

    plot_structure(coord,connect,U, 0.1, true, 'TEST EXAMPLE - Deformed structure with factor = 0.1')
end

%% TEST 1

F_msg = 'F(2*find(coord(:,2)==round(0.2*nely)&coord(:,1)==nelx),1)=-0.01*nelx;';
fixed_dofs_msg = 'fixed_dofs = reshape(2*find(coord(:,1) <= 0.4*nelx & coord(:,2)==nely)+(-1:0),1,[]);';
active_msg = 'active_node = [];';
passive_msg = 'passive_node = find(coord(:,1)>ceil(nelx*0.4) & coord(:,2)>ceil(nely*0.4));';

if caso == 1
    
    %[iter,J, coord, U, connect] = UNVARTOP_2D_compliance_modified(100, 100, 12, 0.36, 0.75, 0, 0.5, ...
    %    F_msg, fixed_dofs_msg, active_msg, passive_msg);
    
    [iter,J, coord, U, connect] = UNVARTOP_2D_compliance_modified(140, 140, 10, 0.36, 0.75, 0, 0.5, ...
        F_msg, fixed_dofs_msg, active_msg, passive_msg);

    plot_structure(coord,connect,U, 0.1, true, 'L-SHAPE STRUCTURE - Deformed structure with factor = 0.1')
end

%% OPTIMAL TEST 1
    
if caso == 1.1
    nels = [10, 30, 50, 60, 70, 80];
    nsteps = [20, 20, 20, 14, 14, 14];
    for i = 1:length(nels)
        [iter,J, coord, U, connect] = ...
            UNVARTOP_2D_compliance_modified(nels(i), nels(i), nsteps(i), 0, 0.36, 0, 0.5, ...
        F_msg, fixed_dofs_msg, active_msg, passive_msg);

        x = coord; Tnod = connect; u = U;
        % Get dimensions
        Ndim = size(x,2);
        Nnodes = size(x,1);
        Nelements = size(Tnod,1);
        NnodesXelement = size(Tnod,2);
        NdofsXnode = 2;

        % X,Y data
        X = x(:,1);
        Y = x(:,2);

        % Displacements
        u = reshape(u,NdofsXnode,Nnodes);
        D = u(1,:).^2 + u(2,:).^2;
        max_Ds(i) = max(D);
        Js(i) = J;
    end
end

%% CASO 2

F_msg = 'F(2*find(coord(:,2)==nely & coord(:,1)==0),1) = -0.01*nelx;';
fixed_dofs_msg = 'fixed_dofs = [reshape(2*find(coord(:,1)==0)-1,1,[]), reshape(2*find(coord(:,1)==nelx & coord(:,2)==0),1,[])];';
active_msg = 'active_node = [];';
passive_msg = 'passive_node = [];';

if caso == 2
    
    %[iter,J, coord, U, connect] = UNVARTOP_2D_compliance_modified(150, 50, 10, 0, 0.6, 0, 1, ...
    %    F_msg, fixed_dofs_msg, active_msg, passive_msg);
    
    [iter,J, coord, U, connect] = UNVARTOP_2D_compliance_modified(190, 76, 10, 0, 0.6, 0, 1, ...
        F_msg, fixed_dofs_msg, active_msg, passive_msg);

    plot_structure(coord,connect,U, 0.1, true, 'MBB BEAM - Deformed structure with factor = 0.1')
end

%% CASO 3

F_msg = 'F(2*find(coord(:,2)>=0.9*nely & coord(:,1)==0) - 1,1) = 0.0001*nelx; F(2*find(coord(:,2)==round(0.9*nely)&coord(:,1)>= 0.9*nelx),2)=0.0001*nelx;';
fixed_dofs_msg = 'fixed_dofs = [reshape(2*find(coord(:,2)==nely),1,[]), reshape(2*find(coord(:,1)==0 & coord(:,2)<=0.1*nely)+(-1:0),1,[])];';
active_msg = 'active_node = [find(coord(:,2)>0.9*nely & coord(:,1)<0.05*nelx); find(coord(:,2)> 0.9*nely & coord(:,2)<=0.95*nely & coord(:,1) >= 0.9*nelx)];';
passive_msg = 'passive_node = [find(coord(:,1)>0.8*nelx & coord(:,1)<0.9*nelx & coord(:,2)>0.8*nely); find(coord(:,1)>=0.9*nelx & coord(:,2)>0.95*nely)];';

if caso == 3

    [iter,J, coord, U, connect] = UNVARTOP_2D_complmechanism_modified(150, 75, 14 , 0, 0.85, -2, 0.5, ...
        F_msg, fixed_dofs_msg, active_msg, passive_msg);
    
    %[iter, J, coord, U, connect] = UNVARTOP_2D_complmechanism_modified(200, 80, 12, 0, 0.7, -2, 0.5, ...
    %    F_msg, fixed_dofs_msg, active_msg, passive_msg);

    plot_structure(coord,connect,U, 0.1, true, 'MBB BEAM - Deformed structure with factor = 0.1')
end



stop = true;

function plot_structure(x,Tnod,u, factor, show_quarters, title_msg)

    % Get dimensions
    Ndim = size(x,2);
    Nnodes = size(x,1);
    Nelements = size(Tnod,1);
    NnodesXelement = size(Tnod,2);
    NdofsXnode = 2;

    % X,Y data
    X = x(:,1);
    Y = x(:,2);

    xmin = min(X);
    xmax = max(X);
    ymin = min(Y);
    ymax = max(Y);

    left = find( x(:,1) == xmin );
    bottom = find( x(:,2) == ymin );
    right  = find( x(:,1) == xmax );
    top = find( x(:,2) == ymax );
    midx = find( x(:,1) == (xmin+xmax)/2 );
    midy = find( x(:,2) == (ymin+ymax)/2 );
    one_fourx = find( x(:,1) == round((xmin+xmax)/4) );
    three_fourx = find( x(:,1) == round(3*(xmin+xmax)/4) );
    one_foury = find( x(:,2) == round((ymin+ymax)/4) );
    three_foury = find( x(:,2) == round(3*(ymin+ymax)/4) );

    % Displacements
    u = reshape(u,NdofsXnode,Nnodes);
    def_X = X + factor * u(1,:)';
    def_Y = Y + factor * u(2,:)';

    figure
    plot(X(bottom),Y(bottom),'k'); hold on; axis equal
    plot(X(left),Y(left),'k'); plot(X(top),Y(top),'k'); 
    plot(X(right),Y(right),'k');
    plot(X(midx),Y(midx)','k');
    plot(X(midy),Y(midy)','k');
    
    plot(def_X(bottom),def_Y(bottom)','r');
    plot(def_X(right),def_Y(right)','r');
    plot(def_X(left),def_Y(left)','r');
    plot(def_X(top),def_Y(top)','r');
    plot(def_X(midx),def_Y(midx)','r');
    plot(def_X(midy),def_Y(midy)','r');
    
    if show_quarters == true
        plot(X(one_fourx),Y(one_fourx)','k');
        plot(X(three_fourx),Y(three_fourx)','k');
        plot(X(one_foury),Y(one_foury)','k');
        plot(X(three_foury),Y(three_foury)','k');
        plot(def_X(one_fourx),def_Y(one_fourx)','r');
        plot(def_X(three_fourx),def_Y(three_fourx)','r');
        plot(def_X(one_foury),def_Y(one_foury)','r');
        plot(def_X(three_foury),def_Y(three_foury)','r');
    end
    title(title_msg,'FontName','latex','Fontsize',12);

    stop =true 

end


function ids = join_ids (id_x, id_y) 
    p = 1;
    ids = [];
    for i = 1:length(id_x)
        for j = 1:length(id_y)
            if id_x(i) == id_y(j)
                ids(p) = i;
            end
        end
    end
end
