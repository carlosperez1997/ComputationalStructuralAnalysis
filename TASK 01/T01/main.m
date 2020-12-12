function main
%%%%%%%%% Computational Engineering (220311/205001) %%%%%%%%%
%%%------------------- 30/11/2020 ------------------------%%%
%%%---------------- Code developed by: -------------------%%%
%%%-------- Jaume Puigdelloses and Carlos Pérez ----------%%%
%%%----------------- MSc Aerospace (MUEA)-----------------%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clc; 
    clear all;
    close all;
    
    % Use of struct to simplify workspace and input/output in functions:
    % - d: data
    % - n: node information
    % - m: matrices / middle solutions
    % - s: solution
    
    %% DATA
    d.g = 9.81;        % Gravity [m/s2]
    d.M_w = 1450;      % Mass of the wing (kg)
    d.M_e = 980;       % Mass of the engine (kg)

    %Young Modulus (MPa)
    d.E_spar = 72.5e9;     
    d.E_rib = 58.1e9;    

    %Poisson coefficient ()
    d.mu_spar = 0.3;
    d.mu_rib = 0.27;

    %Shear Modulus (MPa)
    d.G_spar = d.mu_spar*d.E_spar;
    d.G_rib = d.mu_rib*d.E_rib;

    %Density (kg/m^3)
    d.rho_spar = 2320;
    d.rho_rib = 1850;

    %Cross section (m)
    d.a_spar = 11e-3;      
    d.b_spar = 87e-3;     
    d.t_spar = 5e-3;    

    d.a_rib = 5e-3;      
    d.b_rib = 12e-3;     
    d.t_rib = 3e-3; 

    %Lift and Drag parameters in the spars (N/m)
    d.l1_front = 15e3;
    d.l2_front = 4.5e3;
    d.d1_front = 1e3;

    d.l1_rear = 4.5e3;
    d.l2_rear = 4.5e3;

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
    mat = [% Young M.   Section    Density      Inertia         G           a           b         t
             d.E_spar       0        d.rho_spar      0         d.G_spar     d.a_spar    d.b_spar     d.t_spar; 
             d.E_rib        0        d.rho_rib        0         d.G_rib      d.a_rib     d.b_rib      d.t_rib]';

    run('input_wing');
    
    n.Tmat = Tmat; clear Tmat
    n.T = Tnod; clear Tnod
    n.mat = mat; clear mat
    n.x = x; clear x
    
    %% PREPROCESS
    n.n_deg = 6; % Degrees of freedom per node
    n.n_elem = size(n.T,1); % Number of elements
    n.n_nel = size(n.T,2); % Number of nodes in a beam
    n.n_nod = size(n.x,1); %Total numbre of nodes
    n.n_dof = n.n_nod*n.n_deg; %Total number of degrees of freedom

    %Degrees of freedom connnectivities table
    n.T2=zeros(n.n_deg*n.n_nel,n.n_elem); 

    for e=1:n.n_elem 

        for i=1:n.n_nel
            for j=0:(n.n_deg-1)
                n.T2((i*n.n_deg-n.n_deg+1+j),e)=n.n_deg*n.T(e,i)-n.n_deg+1+j;
            end
        end

    end

    %Aerodynamic forces
    L=zeros(n.n_elem,1);
    D=zeros(n.n_elem,1);

    for e=1:n.n_elem

        %Middle position of each beam
        y_beam=(n.x(n.T(e,1),2)+n.x(n.T(e,2),2))/2;

        if e>=28 && e<=53 %FRONT SPARS

            y_1=n.x(1,2); %Position Node 1
            y_2=n.x(29,2); %Position Node 29
            y_3=n.x(53,2); %Position 53

            %LIFT
            if y_beam<y_2
                L(e)=((d.l1_front+d.l2_front+(d.l1_front-d.l2_front)*cos(pi*((y_beam-y_1)/(y_2-y_1)))))/2;
            else
                L(e)=d.l2_front*cos((pi/2)*((y_beam-y_2)/(y_3-y_2)));
            end

            %DRAG
            D(e)=d.d1_front*(1-((y_beam-y_1)/(y_3-y_1))^2);

        elseif e>=54 && e<=79 %REAR SPARS

            y_1=n.x(5,2); %Position node 5
            y_2=n.x(32,2); %Position node 32
            y_3=n.x(54,2); %Position node 54

            %LIFT
            if y_beam<y_2
                L(e)=0.5*(d.l1_rear+d.l2_rear+(d.l1_rear-d.l2_rear)*cos(pi*((y_beam-y_1)/(y_2-y_1))));
            else
                L(e)=d.l2_rear*cos((pi/2)*((y_beam-y_2)/(y_3-y_2)));
            end

         end
    end    
    
    %% SOLUTION
    
    s = solution(d, fixNod, n, L, D, dat);
    
    disp(sprintf('- Lift = %.2f N', s.Lift))
    disp(sprintf('- Drag = %.2f N', s.Drag))
    disp(sprintf('- Total mass = %.2f kg', s.total_mass))
    disp(sprintf('- Spar mass = %.2f kg', s.mass_spar))
    disp(sprintf('- Ribs mass = %.2f kg', s.mass_rib))
    
    disp('REACTIONS IN NODES 1 and 5:')
    disp('Node 1:')
    disp(s.R_r(1:6)')
    disp('Node 5:')
    disp(s.R_r(7:12)')
    disp('Sum of Reactions:')
    disp(s.R_r(1:6)' + s.R_r(7:12)')
    
    %% PLOTS

    factor=1; %Factor to scale the displacements
    nsub=20; %Number of subdivisions
    plotBeams3D_def(n.x, n.T, nsub, s.l_elem, s.u_int, factor, s.R_e)
    
    id_n = [1]; id_qy = [2]; id_qz = [3];
    id_t = [4]; id_my = [5,11]; id_mz = [6,12];
    
    plotWing(n.x, n.T, 1, s.u, s.u_int, ...
        s.f_int(id_n,:), s.f_int(id_qy,:), s.f_int(id_qz,:),...
        s.f_int(id_t,:), s.f_int(id_my,:), s.f_int(id_mz,:))
 
     %% MOST CRITICAL BEAMS
     
     disp('MOST CRITICAL BEAMS:')
     %Beam with the highest axial force
     [max_n, id_beam_max_n] = max(s.f_int(7,:)); 
     disp(sprintf('- Highest Axial Force (%i) = %.1f N',id_beam_max_n, max_n))
     
     %Beam with the highest Shear Force y
     [max_qy, id_beam_max_qy] = max(s.f_int(8,:)); 
     disp(sprintf('- Highest Shear Force in y (%i) = %.1f N',id_beam_max_qy, max_qy))
     
     %Beam with the highest Shear Force z
     [max_qz, id_beam_max_qz] = max(s.f_int(9,:)); 
     disp(sprintf('- Highest Shear Force in z (%i) = %.1f N',id_beam_max_qz, max_qz))
     
     %Beam with the highest Torsional Moment x
     [max_mx, id_beam_max_mx] = max(s.f_int(10,:)); 
     disp(sprintf('- Highest Torsional Moment in x (%i) = %.1f Nm',id_beam_max_mx, max_mx))
     
     %Beam with the highest Bending Moment y
     [max_my, id_beam_max_my] = max(s.f_int(11,:)); 
     disp(sprintf('- Highest Bending Moment in y (%i) = %.1f Nm',id_beam_max_my, max_my))
     
     %Beam with the highest Bending Moment z
     [max_mz, id_beam_max_mz] = max(s.f_int(12,:)); 
     disp(sprintf('- Highest Bending Moment in z (%i) = %.1f Nm',id_beam_max_mz, max_mz))
     
     nsub=20; %Number of elements's subdivisions to evaluate displacements and rotations
     plotBeams3D(n.T, nsub, s.l_elem, s.u_int, id_beam_max_n) %Axial Force
     plotBeams3D(n.T, nsub, s.l_elem, s.u_int, id_beam_max_qy) %Shear Force y
     plotBeams3D(n.T, nsub, s.l_elem, s.u_int, id_beam_max_qz) %Shear Force z
     plotBeams3D(n.T, nsub, s.l_elem, s.u_int, id_beam_max_mx) %Torsional moment x
     plotBeams3D(n.T, nsub, s.l_elem, s.u_int, id_beam_max_my) %Bending Moment y
     plotBeams3D(n.T, nsub, s.l_elem, s.u_int, id_beam_max_mz) %Bending Moment z
    
end

function [l_elem, Area, I_y, I_z, J,...
    mass_rib, mass_spar, vol_rib, vol_spar, rho_effective,...
    total_mass, total_vol, Lift, Drag] = density_calculus(n_elem, x, T, Tmat, mat, M_w, L, D, dat)
    
    %Initialization
    Lift=0; Drag=0;
    l_elem = zeros(n_elem,1);
    Area = zeros(n_elem,1);
    I_y = zeros(n_elem,1);
    I_z=zeros(n_elem,1);
    J=zeros(n_elem,1);
    mass_spar=0;
    vol_spar=0;
    mass_rib=0;
    vol_rib=0;
    
    %Computations
    for e=1:n_elem 
        %LENGTH OF EACH ELEMENT
        x_1=x(T(e,1),1); y_1=x(T(e,1),2); z_1=x(T(e,1),3);
        
        x_2=x(T(e,2),1); y_2=x(T(e,2),2); z_2=x(T(e,2),3);

        l_elem(e)=sqrt((x_2-x_1)^2+(y_2-y_1)^2+(z_2-z_1)^2);

        %AERODYNAMIC FORCES 
        Lift=Lift+L(e)*l_elem(e); 
        Drag=Drag+D(e)*l_elem(e); 

        %SECTION, SECOND MOMENT OF AREA AND INERTIA OF EACH ELEMENT
        a=mat(6,Tmat(e)); b=mat(7,Tmat(e));
        t=mat(8,Tmat(e)); h=dat(e,1);

        Area(e)=(h-t)*a+2*b*t;
        I_y(e)=(1/12)*a*(h-t)^3+(1/6)*b*t^3+b*t*(h^2/2);
        I_z(e)=(1/12)*(h-t)*a^3+(1/6)*t*b^3;
        J(e)=(1/3)*(2*b*t^3+h*a^3);

        %COMPUTATION OF THE TOTAL MASS AND VOLUME
        if Tmat(e)==1  %SPAR
            mass_spar = mass_spar + mat(3,Tmat(e)) * Area(e)*l_elem(e);
            vol_spar = vol_spar + Area(e)*l_elem(e);

        elseif Tmat(e)==2 %RIB
            mass_rib = mass_rib + mat(3,Tmat(e)) * Area(e)*l_elem(e); 
            vol_rib=vol_rib+Area(e)*l_elem(e);

        end
        
    end
    
    %EFFECTIVE DENSITY
    pseudo_rho=(M_w-mass_spar-mass_rib)/(vol_spar+vol_rib);
    rho_effective=[mat(3,1)+pseudo_rho; mat(3,2)+pseudo_rho];
    
    %TOTAL MASS AND VOLUME
    total_mass=mass_spar+mass_rib;
    total_vol=vol_spar+vol_rib;
end


function [R_e, F_el, K_el] = element(n, l_elem, Area, I_y, I_z, J, Tmat, mat, dat, rho_effective, g, L, D, Drag, M_e)

    Kel = zeros(n.n_nel*n.n_deg,n.n_nel*n.n_deg,n.n_elem);
    Fel = zeros(n.n_nel*n.n_deg, n.n_elem);
    weight_beam=zeros(1,n.n_elem);
    thrust=zeros(1,n.n_elem);
    
    %Computations
    for e=1:n.n_elem
    
       %ROTATION MATRIX
        alpha=dat(e,2);
        beta=dat(e,3);
        gamma=dat(e,4);
   
        R = [ cos(beta)*cos(gamma), cos(beta)*sin(gamma), sin(beta);
            -(sin(alpha)*sin(beta)*cos(gamma))-(cos(alpha)*sin(gamma)), -(sin(alpha)*sin(beta)*sin(gamma))+(cos(alpha)*cos(gamma)), sin(alpha)*cos(beta);
            -(cos(alpha)*sin(beta)*cos(gamma))+(sin(alpha)*sin(gamma)), -(cos(alpha)*sin(beta)*sin(gamma))-(sin(alpha)*cos(gamma)), cos(alpha)*cos(beta)];
        
        Re = zeros(12);
        
        Re(1:3,1:3) = R; Re(4:6,4:6) = R;
        Re(7:9,7:9) = R; Re(10:12,10:12) = R;
        
        R_e(:,:,e) = Re;
   
    %FORCES VECTOR
    weight_beam(e)=rho_effective(Tmat(e))*Area(e)*g; %Weight of each beam
    
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
     
    K_e_local=[((E*Area(e))/l_elem(e)),    0,      0,      0,      0,      0,      (-(E*Area(e))/l_elem(e)),    0,      0,      0,      0,      0;
            0, ((12*E*I_z(e))/(l_elem(e)^3)),    0,      0,      0,  ((6*E*I_z(e))/(l_elem(e)^2)),    0, ((-12*E*I_z(e))/(l_elem(e)^3)),    0,      0,      0,  ((6*E*I_z(e))/(l_elem(e)^2));
            0, 0,   ((12*E*I_y(e))/(l_elem(e)^3)),    0,  ((-6*E*I_y(e))/(l_elem(e)^2)),    0, 0, 0, ((-12*E*I_y(e))/(l_elem(e)^3)),    0,  ((-6*E*I_y(e))/(l_elem(e)^2)),    0;
            0,  0,  0, ((G*J(e))/l_elem(e)),   0,      0,      0,      0,      0,  ((-G*J(e))/l_elem(e)),  0,      0;
            0, 0,   ((-6*E*I_y(e))/(l_elem(e)^2)),    0,  ((4*E*I_y(e))/(l_elem(e))),    0, 0, 0, ((6*E*I_y(e))/(l_elem(e)^2)),    0,  ((2*E*I_y(e))/(l_elem(e))),    0;
            0, ((6*E*I_z(e))/(l_elem(e)^2)),    0,      0,      0,  ((4*E*I_z(e))/(l_elem(e))),    0, ((-6*E*I_z(e))/(l_elem(e)^2)),    0,      0,      0,  ((2*E*I_z(e))/(l_elem(e)));
            ((-E*Area(e))/l_elem(e)),    0,      0,      0,      0,      0,      ((E*Area(e))/l_elem(e)),    0,      0,      0,      0,      0;
            0, ((-12*E*I_z(e))/(l_elem(e)^3)),    0,      0,      0,  ((-6*E*I_z(e))/(l_elem(e)^2)),    0, ((12*E*I_z(e))/(l_elem(e)^3)),    0,      0,      0,  ((-6*E*I_z(e))/(l_elem(e)^2));
            0, 0,   ((-12*E*I_y(e))/(l_elem(e)^3)),    0,  ((6*E*I_y(e))/(l_elem(e)^2)),    0, 0, 0, ((12*E*I_y(e))/(l_elem(e)^3)),    0,  ((6*E*I_y(e))/(l_elem(e)^2)),    0;
            0,  0,  0, ((-G*J(e))/l_elem(e)),   0,      0,      0,      0,      0,  ((G*J(e))/l_elem(e)),  0,      0;
            0, 0,   ((-6*E*I_y(e))/(l_elem(e)^2)),    0,  ((2*E*I_y(e))/(l_elem(e))),    0, 0, 0, ((6*E*I_y(e))/(l_elem(e)^2)),    0,  ((4*E*I_y(e))/(l_elem(e))),    0;
            0, ((6*E*I_z(e))/(l_elem(e)^2)),    0,      0,      0,  ((2*E*I_z(e))/(l_elem(e))),    0, ((-6*E*I_z(e))/(l_elem(e)^2)),    0,      0,      0,  ((4*E*I_z(e))/(l_elem(e)));
            ];

    K_e=R_e(:,:,e).'*K_e_local*R_e(:,:,e);
    
    for r=1:n.n_nel*n.n_deg
        for s=1:n.n_nel*n.n_deg
            
            K_el(r,s,e)=K_e(r,s);
            
        end
    end
 end

 % Thrust
   Thrust = F_el(1,7) + F_el(7,7) + F_el(1,8) + F_el(1,8);
   disp(sprintf('- Thrust = %.2f N',Thrust)) 
    
end

function [K_G, F_ext] = assembly(n, T2, K_el, F_el)
    
    %Initialization
    K_G=zeros(n.n_dof,n.n_dof);
    F_ext=zeros(n.n_dof,1);
    
    %Computation
    for e=1:n.n_elem
        for i=1:(n.n_nel*n.n_deg)

            F_ext(T2(i,e))=F_ext(T2(i,e))+F_el(i,e);

            for j=1:(n.n_nel*n.n_deg)

                K_G(T2(i,e),T2(j,e))=K_G(T2(i,e),T2(j,e))+ K_el(i,j,e);

            end
        end
    end
end


function s = solution(d, fixNod,  n, L, D, dat)

%DENSITY AND GEOMETRICAL PARAMETERS
[l_elem, Area, I_y, I_z, J, ...
    mass_rib, mass_spar, vol_rib, vol_spar, ...
    rho_effective, total_mass, total_vol, Lift, Drag] = ...
    density_calculus (n.n_elem, n.x, n.T, n.Tmat, n.mat, d.M_w, L, D, dat);

%COMPUTATION OF THE ELEMENTAL FORCE AND STIFFNESS MATRICES
[R_e, F_el, K_el] = ...
    element(n, l_elem, Area, I_y, I_z, J, n.Tmat, ...
    n.mat, dat, rho_effective, d.g, L, D, Drag, d.M_e);

%MATRICES' ASSEMBLY 
[K_G, F_ext] = assembly(n, n.T2, K_el, F_el);

%SOLVE THE SYSTEM
[R, u, R_r] = systemEqs(n, K_G, F_ext, fixNod);

%INTERNAL FORCES
for e=1:n.n_elem
  
    for i=1:(n.n_nel*n.n_deg)
        u_elem(i,1)=u(1,n.T2(i,e));
    end
    
    u_int(:,e)=R_e(:,:,e)*u_elem;
    f_int(:,e)=K_el(:,:,e)*u_elem;

end

%SOLUTION'S STRUCT
s.l_elem = l_elem;
s.Lift = Lift;
s.Drag = Drag;
s.mass_rib = mass_rib;
s.mass_spar = mass_spar;
s.vol_rib = vol_rib; 
s.vol_spar = vol_spar;
s.rho_effective = rho_effective;
s.total_mass = total_mass;
s.total_vol = total_vol;
s.u = u;
s.u_int = u_int;
s.f_int = f_int; 
s.R = R; 
s.R_r = R_r; 
s.R_e = R_e; 

end

function [R, u, R_r] = systemEqs(n, K_G, F_ext, fixNod)

    vr = zeros(1,size(fixNod,2));
    ur = zeros(size(fixNod,2),1);

    for i=1:size(fixNod,2)
        u_r(i,1)= fixNod(3,i);
        v_r(1,i)=n.n_deg*fixNod(1,i)-n.n_deg+fixNod(2,i); 
    end
    
    v_l = setdiff(1:n.n_dof,v_r);
    
    K_ll = K_G(v_l,v_l);
    K_lr = K_G(v_l,v_r);
    K_rl = K_G(v_r,v_l);
    K_rr = K_G(v_r,v_r);
    f_l = F_ext(v_l,1);
    f_r = F_ext(v_r,1);

    u_l = K_ll\(f_l - K_lr*u_r);
    R_r = K_rr*u_r + K_rl*u_l - f_r;

    % Displacement
    u(v_l) = u_l;
    u(v_r) = u_r;
    
    % Reactions
    R(v_l) = 0;
    R(v_r) = R_r;

end



