function main
    clc; 
    clear all;
    close all;
    
    %% DATA
    g = 9.81;        % Gravity [m/s2]
    M_w = 1450;      % Mass of the wing (kg)
    M_e = 980;       % Mass of the engine (kg)

    %Young Modulus (MPa)
    E_spar = 72.5e9;     
    E_rib = 58.1e9;    

    %Poisson coefficient ()
    mu_spar = 0.3;
    mu_rib = 0.27;

    %Shear Modulus (MPa)
    G_spar = mu_spar*E_spar;
    G_rib = mu_rib*E_rib;

    %Density (kg/m^3)
    rho_spar = 2320;
    rho_rib = 1850;

    %Cross section (m)
    a_spar = 11e-3;      
    b_spar = 87e-3;     
    t_spar = 5e-3;    

    a_rib = 5e-3;      
    b_rib = 12e-3;     
    t_rib = 3e-3; 

    %Lift and Drag parameters in the spars (N/m)
    l1_front = 15e3;
    l2_front = 4.5e3;
    d1_front = 1e3;

    l1_rear = 4.5e3;
    l2_rear = 4.5e3;

    % Fix nodes specification -> nodes 1 and 5
    fixNod = [% Node   DOF  Magnitude
                   1     1      0;
                   1     2      0;
                   1     3      0;
                   1     4      0;
                   1     5      0;
                   1     6      0;
                   5     1      0;
                   5     2      0;
                   5     3      0;
                   5     4      0;
                   5     5      0;
                   5     6      0]';

    % Data Matrix
    mat = [% Young M.   Section    Density      Inertia     G        a           b         t
             E_spar       0        rho_spar      0         G_spar     a_spar    b_spar     t_spar; 
             E_rib        0        rho_rib        0         G_rib      a_rib     b_rib      t_rib]';

    run('input_wing');
    
    %% PREPROCESS
    n_deg = 6; % Degrees of freedom per node
    n_elem = size(Tnod,1); % Number of elements
    n_nel = size(Tnod,2); % Number of nodes in a beam
    n_nod = size(x,1); %Total numbre of nodes
    n_dof = n_nod*n_deg; %Total number of degrees of freedom

    %Degrees of freedom connnectivities table
    T2=zeros(n_deg*n_nel,n_elem); 

    for e=1:n_elem 

        for i=1:n_nel
            for j=0:(n_deg-1)
                T2((i*n_deg-n_deg+1+j),e)=n_deg*Tnod(e,i)-n_deg+1+j;
            end
        end

    end

    %Aerodynamic forces
    L=zeros(n_elem,1);
    D=zeros(n_elem,1);

    for e=1:n_elem

        %Middle position of each beam
        y_beam=(x(Tnod(e,1),2)+x(Tnod(e,2),2))/2;

        if e>=28 && e<=53 %FRONT SPARS

            y_1=x(1,2); %Position Node 1
            y_2=x(29,2); %Position Node 29
            y_3=x(53,2); %Position 53

            %LIFT
            if y_beam<y_2
                L(e)=((l1_front+l2_front+(l1_front-l2_front)*cos(pi*((y_beam-y_1)/(y_2-y_1)))))/2;
            else
                L(e)=l2_front*cos((pi/2)*((y_beam-y_2)/(y_3-y_2)));
            end

            %DRAG
            D(e)=d1_front*(1-((y_beam-y_1)/(y_3-y_1))^2);

        elseif e>=54 && e<=79 %REAR SPARS

            y_1=x(5,2); %Position node 5
            y_2=x(32,2); %Position node 32
            y_3=x(54,2); %Position node 54

            %LIFT
            if y_beam<y_2
                L(e)=0.5*(l1_rear+l2_rear+(l1_rear-l2_rear)*cos(pi*((y_beam-y_1)/(y_2-y_1))));
            else
                L(e)=l2_rear*cos((pi/2)*((y_beam-y_2)/(y_3-y_2)));
            end

         end
    end
    % Tot canviat!!
    
    %% SOLUTION
    
    [u, f_int, u_int, total_mass, total_vol, mass_rib, vol_rib, mass_spar, vol_spar, Lift, Drag , React, React_r, l_elem, R_e] = ...
    solution(g, M_w, M_e, fixNod,  mat, Tmat, x, Tnod, T2, L, D, dat, n_deg, n_elem, n_nel, n_dof)
    
    % tot correcte!
    
    %% PLOTS

    factor=1;% Number to scale/amplify the displacements
    nsub=10; %Number of elements's subdivisions to evaluate displacements and rotations
    plotBeams3D_def(x,Tnod, nsub, l_elem, u_int, factor, R_e)
    
    id_n = [1]; id_qy = [2]; id_qz = [3];
    id_t = [4]; id_my = [5,11]; id_mz = [6,12];
    
    plotWing(x,Tnod, 1, u, u_int, ...
        f_int(id_n,:), f_int(id_qy,:), f_int(id_qz,:),...
        f_int(id_t,:), f_int(id_my,:), f_int(id_mz,:))
    
    %  n	 Axial force of element e in local reference frame [1 x Nelements]
    %  qy	 Shear force in y-direction of element e in local reference frame [1 x Nelements]
    %  qz	 Shear force in z-direction of element e in local reference frame [1 x Nelements]
    %  t	 Torsion moment of element e in local reference frame [1 x Nelements]
    %  my	 Bending moment in y-direction of element e in local reference frame [2 x Nelements]
    %  mz	 Bending moment in z-direction of element e in local reference frame [2 x Nelements]
    
%        fint(1,e) : F  for node 1 of element e in local reference frame
%        fint(2,e) : Qy for node 1 of element e in local reference frame
%        fint(3,e) : Qz for node 1 of element e in local reference frame
%        fint(4,e) : T  for node 1 of element e in local reference frame
%        fint(5,e) : My for node 1 of element e in local reference frame
%        fint(6,e) : Mz for node 1 of element e in local reference frame
%        fint(7,e) : F  for node 2 of element e in local reference frame
%        fint(8,e) : Qy for node 2 of element e in local reference frame
%        fint(9,e) : Qz for node 2 of element e in local reference frame
%        fint(10,e): T  for node 2 of element e in local reference frame
%        fint(11,e): My for node 2 of element e in local reference frame
%        fint(12,e): Mz for node 2 of element e in local reference frame
    
%     plotWing(x, Tnod, 1, u, u_int, n, qy, qz, t, my, mz)
% 
%     plotWing(x, Tnod, u.', u_int, f_int)
% 
%     %% MOST CRITICAL BEAMS
% 
%     % Axial deformation
%     ux = u_int(7,:);
%     % Shear Y
%     Qy = f_int(8,:);
%     % Shear Z
%     Qz = f_int(9,:);
%     % Bending Moment Y
%     My = f_int(11,:);
%     % Bending Moment Z
%     Mz = f_int(12,:);
% 
%     [max_ux, element_ux]=max(ux);
%     [max_Qy, element_Qy]=max(Qy);
%     [max_Qz, element_Qz]=max(Qz);
%     [max_My, element_My]=max(My);
%     [max_Mz, element_Mz]=max(Mz);
% 
%     nsub=10; %Number of elements's subdivisions to evaluate displacements and rotations
%     plotBeams3D(T.',nsub,l_elem,u_int,element_ux) %22
%     plotBeams3D(T.',nsub,l_elem,u_int,element_Mz) %54
%     plotBeams3D(T.',nsub,l_elem,u_int,element_Qz) %28
% 
%     %% CLEAR UNNECESSARY VARIABLES
%     clear a_rib a_spar b_spar b_rib d1_front dat e E_rib E_spar f_int fixNod...
%         g G_rib G_spar i j l1_front l1_rear l2_front l2_rear M_e M_w mat...
%         n_dof n_el n_i n_ne n_nod rho_rib rho_spar T T1 T2 t_rib t_spar Tmat...
%         u_int x y_1 y_2 y_3 y_beam l_elem R_e ux Qy Qz My Mz
    
end

function [l_elem, Lift, Drag, A, I_y, I_z, J, mass_rib, mass_spar, vol_rib, vol_spar, effective_rho, total_mass, total_vol] = density_calculus(n_elem, x, T, Tmat, mat, M_w, L, D, dat)
    
    %Initialization
    Lift=0;
    Drag=0;
    l_elem=zeros(n_elem,1);
    A=zeros(n_elem,1);
    I_y=zeros(n_elem,1);
    I_z=zeros(n_elem,1);
    J=zeros(n_elem,1);
    mass_spar=0;
    vol_spar=0;
    mass_rib=0;
    vol_rib=0;
    
    %Computations
    for e=1:n_elem 
        %LENGTH OF EACH ELEMENT
        %Node 1
        x_1=x(T(e,1),1);
        y_1=x(T(e,1),2);
        z_1=x(T(e,1),3);

        %Node 2
        x_2=x(T(e,2),1);
        y_2=x(T(e,2),2);
        z_2=x(T(e,2),3);

        %Length
        l_elem(e)=sqrt((x_2-x_1)^2+(y_2-y_1)^2+(z_2-z_1)^2);

        %TOTAL AERODYNAMIC FORCES made by the wing
        Lift=Lift+L(e)*l_elem(e); 
        Drag=Drag+D(e)*l_elem(e); 

        %SECTION, SECOND MOMENT OF AREA AND INERTIA OF EACH ELEMENT
        a=mat(6,Tmat(e));
        b=mat(7,Tmat(e));
        t=mat(8,Tmat(e));
        h=dat(e,1);

        A(e)=(h-t)*a+2*b*t;
        I_y(e)=(1/12)*a*(h-t)^3+(1/6)*b*t^3+b*t*(h^2/2);
        I_z(e)=(1/12)*(h-t)*a^3+(1/6)*t*b^3;
        J(e)=(1/3)*(2*b*t^3+h*a^3);

        %COMPUTATION OF THE TOTAL MASS AND VOLUME
        if Tmat(e)==1  %SPAR
            mass_spar = mass_spar + mat(3,Tmat(e)) * A(e)*l_elem(e);
            vol_spar = vol_spar + A(e)*l_elem(e);

        elseif Tmat(e)==2 %RIB
            mass_rib = mass_rib + mat(3,Tmat(e)) * A(e)*l_elem(e); 
            vol_rib=vol_rib+A(e)*l_elem(e);

        end
        
    end
    
    %EFFECTIVE DENSITY
    pseudo_rho=(M_w-mass_spar-mass_rib)/(vol_spar+vol_rib);
    effective_rho=[mat(3,1)+pseudo_rho; mat(3,2)+pseudo_rho];
    
    %TOTAL MASS AND VOLUME
    total_mass=mass_spar+mass_rib;
    total_vol=vol_spar+vol_rib;
end


function [R_e, F_el, K_el] = element(n_nel, n_deg, n_elem, l_elem, A, I_y, I_z, J, Tmat, mat, dat, effective_rho, g, L, D, Drag, M_e)

    %Initialization
    Kel = zeros(n_nel*n_deg,n_nel*n_deg,n_elem);
    Fel = zeros(n_nel*n_deg, n_elem);
    weight_beam=zeros(1,n_elem);
    thrust=zeros(1,n_elem);
    
    %Computations
    for e=1:n_elem
    
   %ROTATION MATRIX
    %Angles for the rotation matrix
    alpha=dat(e,2);
    beta=dat(e,3);
    gamma=dat(e,4);
   
   R_e(:,:,e)=[cos(beta)*cos(gamma)                              cos(beta)*sin(gamma)                           sin(beta) 0  0  0  0  0  0  0  0  0;
       -(sin(alpha)*sin(beta)*cos(gamma))-(cos(alpha)*sin(gamma))   -(sin(alpha)*sin(beta)*sin(gamma))+(cos(alpha)*cos(gamma))  sin(alpha)*cos(beta) 0  0  0  0  0  0  0  0  0;
       -(cos(alpha)*sin(beta)*cos(gamma))+(sin(alpha)*sin(gamma))   -(cos(alpha)*sin(beta)*sin(gamma))-(sin(alpha)*cos(gamma))  cos(alpha)*cos(beta) 0  0  0  0  0  0  0  0  0;
       0  0  0 cos(beta)*cos(gamma)                              cos(beta)*sin(gamma)                           sin(beta)  0  0  0  0  0  0;
       0  0  0 -(sin(alpha)*sin(beta)*cos(gamma))-(cos(alpha)*sin(gamma))   -(sin(alpha)*sin(beta)*sin(gamma))+(cos(alpha)*cos(gamma))  sin(alpha)*cos(beta)  0  0  0  0  0  0;
       0  0  0 -(cos(alpha)*sin(beta)*cos(gamma))+(sin(alpha)*sin(gamma))   -(cos(alpha)*sin(beta)*sin(gamma))-(sin(alpha)*cos(gamma))  cos(alpha)*cos(beta)  0  0  0  0  0  0;
       0  0  0  0  0  0   cos(beta)*cos(gamma)                              cos(beta)*sin(gamma)                           sin(beta)  0  0  0;
       0  0  0  0  0  0   -(sin(alpha)*sin(beta)*cos(gamma))-(cos(alpha)*sin(gamma))   -(sin(alpha)*sin(beta)*sin(gamma))+(cos(alpha)*cos(gamma))  sin(alpha)*cos(beta)  0  0  0;
       0  0  0  0  0  0   -(cos(alpha)*sin(beta)*cos(gamma))+(sin(alpha)*sin(gamma))   -(cos(alpha)*sin(beta)*sin(gamma))-(sin(alpha)*cos(gamma))  cos(alpha)*cos(beta)  0  0  0;
       0  0  0  0  0  0  0  0  0   cos(beta)*cos(gamma)                              cos(beta)*sin(gamma)                           sin(beta);
       0  0  0  0  0  0  0  0  0   -(sin(alpha)*sin(beta)*cos(gamma))-(cos(alpha)*sin(gamma))   -(sin(alpha)*sin(beta)*sin(gamma))+(cos(alpha)*cos(gamma))  sin(alpha)*cos(beta);
       0  0  0  0  0  0  0  0  0   -(cos(alpha)*sin(beta)*cos(gamma))+(sin(alpha)*sin(gamma))   -(cos(alpha)*sin(beta)*sin(gamma))-(sin(alpha)*cos(gamma))  cos(alpha)*cos(beta)]; 

    %FORCES VECTOR
    weight_beam(e)=effective_rho(Tmat(e))*A(e)*g; %Weight of each beam
    
    if e==7 || e==8 %Elements where the engine is placed
        
        weight_beam(e)=weight_beam(e)+((M_e*g)/(2*l_elem(e)));
        thrust(e)=Drag/(2*l_elem(e));
    
    end
    
    %Load applied in each dof
    q_e=[(D(e,1)-thrust(1,e));
            0;
            (L(e,1)-weight_beam(1,e));
            0;
            0;
            0;
            0;
            0;
            0;
            0;
            0;
            0];
        
    q_e_local=R_e(:,:,e)*q_e;
    
    f_e_local=[q_e_local(1)*l_elem(e)/2;
               q_e_local(2)*l_elem(e)/2;
               q_e_local(3)*l_elem(e)/2;
               0;
               -q_e_local(3)*l_elem(e)*l_elem(e)/12;
               q_e_local(2)*l_elem(e)^2/12;
               q_e_local(1)*l_elem(e)/2;  
               q_e_local(2)*l_elem(e)/2;
               q_e_local(3)*l_elem(e)/2;
               0;
               q_e_local(3)*l_elem(e)^2/12;
               -q_e_local(2)*l_elem(e)^2/12];

    F_el(:,e)=R_e(:,:,e).'*f_e_local;
    
    %STIFFNESS MATRIX
    E=mat(1,Tmat(e)); %Yound modulus of the element
    G=mat(5,Tmat(e)); %Shear modulus of the element
     
    K_e_local=[((E*A(e))/l_elem(e)),    0,      0,      0,      0,      0,      (-(E*A(e))/l_elem(e)),    0,      0,      0,      0,      0;
            0, ((12*E*I_z(e))/(l_elem(e)^3)),    0,      0,      0,  ((6*E*I_z(e))/(l_elem(e)^2)),    0, ((-12*E*I_z(e))/(l_elem(e)^3)),    0,      0,      0,  ((6*E*I_z(e))/(l_elem(e)^2));
            0, 0,   ((12*E*I_y(e))/(l_elem(e)^3)),    0,  ((-6*E*I_y(e))/(l_elem(e)^2)),    0, 0, 0, ((-12*E*I_y(e))/(l_elem(e)^3)),    0,  ((-6*E*I_y(e))/(l_elem(e)^2)),    0;
            0,  0,  0, ((G*J(e))/l_elem(e)),   0,      0,      0,      0,      0,  ((-G*J(e))/l_elem(e)),  0,      0;
            0, 0,   ((-6*E*I_y(e))/(l_elem(e)^2)),    0,  ((4*E*I_y(e))/(l_elem(e))),    0, 0, 0, ((6*E*I_y(e))/(l_elem(e)^2)),    0,  ((2*E*I_y(e))/(l_elem(e))),    0;
            0, ((6*E*I_z(e))/(l_elem(e)^2)),    0,      0,      0,  ((4*E*I_z(e))/(l_elem(e))),    0, ((-6*E*I_z(e))/(l_elem(e)^2)),    0,      0,      0,  ((2*E*I_z(e))/(l_elem(e)));
            ((-E*A(e))/l_elem(e)),    0,      0,      0,      0,      0,      ((E*A(e))/l_elem(e)),    0,      0,      0,      0,      0;
            0, ((-12*E*I_z(e))/(l_elem(e)^3)),    0,      0,      0,  ((-6*E*I_z(e))/(l_elem(e)^2)),    0, ((12*E*I_z(e))/(l_elem(e)^3)),    0,      0,      0,  ((-6*E*I_z(e))/(l_elem(e)^2));
            0, 0,   ((-12*E*I_y(e))/(l_elem(e)^3)),    0,  ((6*E*I_y(e))/(l_elem(e)^2)),    0, 0, 0, ((12*E*I_y(e))/(l_elem(e)^3)),    0,  ((6*E*I_y(e))/(l_elem(e)^2)),    0;
            0,  0,  0, ((-G*J(e))/l_elem(e)),   0,      0,      0,      0,      0,  ((G*J(e))/l_elem(e)),  0,      0;
            0, 0,   ((-6*E*I_y(e))/(l_elem(e)^2)),    0,  ((2*E*I_y(e))/(l_elem(e))),    0, 0, 0, ((6*E*I_y(e))/(l_elem(e)^2)),    0,  ((4*E*I_y(e))/(l_elem(e))),    0;
            0, ((6*E*I_z(e))/(l_elem(e)^2)),    0,      0,      0,  ((2*E*I_z(e))/(l_elem(e))),    0, ((-6*E*I_z(e))/(l_elem(e)^2)),    0,      0,      0,  ((4*E*I_z(e))/(l_elem(e)));
            ];

    K_e=R_e(:,:,e).'*K_e_local*R_e(:,:,e);
    
    for r=1:n_nel*n_deg
        for s=1:n_nel*n_deg
            
            K_el(r,s,e)=K_e(r,s);
            
        end
    end
end


end

function [K_G, F_ext] = assembly(n_dof, n_nel, n_deg, n_elem, T2, K_el, F_el)
    
    %Initialization
    K_G=zeros(n_dof,n_dof);
    F_ext=zeros(n_dof,1);
    
    %Computation
    for e=1:n_elem
        for i=1:(n_nel*n_deg)

            F_ext(T2(i,e))=F_ext(T2(i,e))+F_el(i,e);

            for j=1:(n_nel*n_deg)

                K_G(T2(i,e),T2(j,e))=K_G(T2(i,e),T2(j,e))+ K_el(i,j,e);

            end
        end
    end
end


function [u, f_int, u_int, total_mass, total_vol, mass_rib, vol_rib, mass_spar, vol_spar, Lift, Drag , React, React_r, l_elem, R_e] = ...
    solution(g, M_w, M_e, fixNod,  mat, Tmat, x, T, T2, L, D, dat, n_deg, n_elem, n_nel, n_dof)

%PREVIOUS CALCULATIONS
[l_elem, Lift, Drag, A, I_y, I_z, J, ...
    mass_rib, mass_spar, vol_rib, vol_spar, ...
    effective_rho, total_mass, total_vol] = ...
    density_calculus (n_elem, x, T, Tmat, mat, M_w, L, D, dat);

%BIEN!

%COMPUTATION OF THE ELEMENTAL FORCE AND STIFFNESS MATRICES
[R_e, F_el, K_el] = element(n_nel, n_deg, n_elem, l_elem, A, I_y, I_z, J, Tmat, mat, dat, effective_rho, g, L, D, Drag, M_e);

%ASSEMBLY OF THE MATRICES
[K_G, F_ext] = assembly(n_dof, n_nel, n_deg, n_elem, T2, K_el, F_el);

%SOLVE THE SYSTEM
[React, u, React_r] = systemEqs(K_G, F_ext, fixNod, n_deg, n_dof);

%Internal forces and displacements
for e=1:n_elem
  
    for i=1:(n_nel*n_deg)
        u_elem(i,1)=u(1,T2(i,e));
    end
    
    u_int(:,e)=R_e(:,:,e)*u_elem;
    f_int(:,e)=K_el(:,:,e)*u_elem;

end


end

function [React, u, React_r] = systemEqs(K_G, F_ext, fixNod, n_i, n_dof)

    %Initialization
    vr = zeros(1,size(fixNod,2));
    ur = zeros(size(fixNod,2),1);

    for i=1:size(fixNod,2)
        u_r(i,1)= fixNod(3,i);
        v_r(1,i)=n_i*fixNod(1,i)-n_i+fixNod(2,i); 
    end
    
    v_l = setdiff(1:n_dof,v_r);
    
    %Computation
    K_ll = K_G(v_l,v_l);
    K_lr = K_G(v_l,v_r);
    K_rl = K_G(v_r,v_l);
    K_rr = K_G(v_r,v_r);
    f_l = F_ext(v_l,1);
    f_r = F_ext(v_r,1);

    u_l = K_ll\(f_l - K_lr*u_r);
    React_r = K_rr*u_r + K_rl*u_l - f_r;

    u(v_l) = u_l;
    u(v_r) = u_r;

    % Assembly of global reactions
    React(v_l) = 0;
    React(v_r) = React_r;

end



