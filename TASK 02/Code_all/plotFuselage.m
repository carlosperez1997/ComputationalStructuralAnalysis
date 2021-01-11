classdef plotFuselage < handle
% Call function as:  plotFuselage(xnodes,Tbeams,Tplates,Tmat_beams,Tmat_plates,u,uint_beams,N,Qy,Qz,T,My,Mz,uint_plates,mat_plates)
% With inputs:
%	- xnodes: 		nodal coordinates matrix 												[Nnodes x Ndim]
%	- Tbeams: 		nodal connectivities matrix for beams 									[Nel_beams x NnodesXelement_beams]
%	- Tplates:		nodal connectivities matrix for plates 									[Nel_plates x NnodesXelement_plates]
%	- Tmat_beams:	
%	- Tmat_plates:
%	- u:			global displacements/rotations vector 									[Ndofs x 1]
%	- uint_beams:	matrix with element displacements and rotations in LOCAL axes for beams	[12 x Nel_beams] 
%        uint_beams(1,e) : ux for node 1 of element e in local reference frame
%        uint_beams(2,e) : uy for node 1 of element e in local reference frame
%        uint_beams(3,e) : uz for node 1 of element e in local reference frame
%        uint_beams(4,e) : angle_x for node 1 of element e in local reference frame
%        uint_beams(5,e) : angle_y for node 1 of element e in local reference frame
%        uint_beams(6,e) : angle_z for node 1 of element e in local reference frame
%        uint_beams(7,e) : ux for node 2 of element e in local reference frame
%        uint_beams(8,e) : uy for node 2 of element e in local reference frame
%        uint_beams(9,e) : uz for node 2 of element e in local reference frame
%        uint_beams(10,e): angle_x for node 2 of element e in local reference frame
%        uint_beams(11,e): angle_y for node 2 of element e in local reference frame
%        uint_beams(12,e): angle_z for node 2 of element e in local reference frame
%	- N:			Axial force of element e in local reference frame 						[1 x Nel_beams]
% 	- Qy:			Shear force in y-direction of element e in local reference frame 		[1 x Nel_beams]
%	- Qz:			Shear force in z-direction of element e in local reference frame 		[1 x Nel_beams]
%	- T:			Torsion moment of element e in local reference frame 					[1 x Nel_beams]
%	- My:			Bending moment in y-direction of element e in local reference frame 	[2 x Nel_beams]
%	- Mz:			Bending moment in z-direction of element e in local reference frame 	[2 x Nel_beams]
%	- uint_plates: 	matrix with element displacements and rotations in LOCAL axes for plates[24 x Nel_plates] 
%        uint_plates(1,e) : ux for node 1 of element e in local reference frame
%        uint_plates(2,e) : uy for node 1 of element e in local reference frame
%        uint_plates(3,e) : uz for node 1 of element e in local reference frame
%        uint_plates(4,e) : angle_x for node 1 of element e in local reference frame
%        uint_plates(5,e) : angle_y for node 1 of element e in local reference frame
%        uint_plates(6,e) : angle_z for node 1 of element e in local reference frame
%        uint_plates(7,e) : ux for node 2 of element e in local reference frame
%        uint_plates(8,e) : uy for node 2 of element e in local reference frame
%        uint_plates(9,e) : uz for node 2 of element e in local reference frame
%        uint_plates(10,e): angle_x for node 2 of element e in local reference frame
%        uint_plates(11,e): angle_y for node 2 of element e in local reference frame
%        uint_plates(12,e): angle_z for node 2 of element e in local reference frame
%        uint_plates(13,e) : ux for node 3 of element e in local reference frame
%        uint_plates(14,e) : uy for node 3 of element e in local reference frame
%        uint_plates(15,e) : uz for node 3 of element e in local reference frame
%        uint_plates(16,e) : angle_x for node 3 of element e in local reference frame
%        uint_plates(17,e) : angle_y for node 3 of element e in local reference frame
%        uint_plates(18,e) : angle_z for node 3 of element e in local reference frame
%        uint_plates(19,e) : ux for node 4 of element e in local reference frame
%        uint_plates(20,e) : uy for node 4 of element e in local reference frame
%        uint_plates(21,e) : uz for node 4 of element e in local reference frame
%        uint_plates(22,e): angle_x for node 4 of element e in local reference frame
%        uint_plates(23,e): angle_y for node 4 of element e in local reference frame
%        uint_plates(24,e): angle_z for node 4 of element e in local reference frame
%	- mat_plates:	material properties for plates (columns: rho E nu h, row: skin,floor)	[2 x 4]	properties (Access=public)
	
    properties (Access=public)	
	
		xbeams
		xplates
		Tbeams
		Tplates
		Tmat_beams
		Tmat_plates
		ubeams
		uplates
		
		uint_beams
		N
		Qy
		Qz
		T
		My
		Mz
		
		uint_plates
		mat_plates
		sig_x
		sig_y
		sig_xy
		sig_vm
		p
		
		
	end
	properties (Access=protected)
		
		fg_handle
		ax_handle
		
		obj_stringer
		obj_frame
		obj_floor_spars
		obj_floor_plates
		obj_skin_plates
		
		pnl_selection
		chb_stringer
		chb_frame
		chb_floor_spars
		chb_floor_plates
		chb_skin_plates
		
		chb_sym
		
		txt_scale
		txtedit_scale
		scale			= 1;
		
		dd_color
		
		coord_beams_mapping
		coord_plates_mapping
		connect_beams_mapping
		connect_plates_mapping
		mat_beams_mapping
		mat_plates_mapping
		U_beams_mapping
		U_plates_mapping
		color_beams_mapping
		color_plates_mapping
		
	end
	
	methods
		function obj = plotFuselage(varargin)
			
			obj.Tbeams = varargin{2};
			obj.Tplates = varargin{3};
			
			obj.xbeams = varargin{1}; obj.xbeams = obj.xbeams(reshape(obj.Tbeams',[],1),:);
			obj.xplates = varargin{1}; obj.xplates = obj.xplates(reshape(obj.Tplates',[],1),:);
			
			u_rshp = reshape(varargin{6},6,[]); 
			obj.ubeams = reshape(u_rshp(1:3,reshape(obj.Tbeams',[],1)),[],1);
			obj.uplates = reshape(u_rshp(1:3,reshape(obj.Tplates',[],1)),[],1);
			
			obj.Tbeams = reshape(1:numel(obj.Tbeams),2,[])';
			obj.Tplates = reshape(1:numel(obj.Tplates),4,[])';
			
			obj.Tmat_beams = varargin{4};
			obj.Tmat_plates = varargin{5};
			
			%beams local variables (elemental)
			obj.uint_beams = varargin{7};
			obj.N = repmat(reshape(varargin{8},1,[]),2,1);
			obj.Qy = repmat(reshape(varargin{9},1,[]),2,1);
			obj.Qz = repmat(reshape(varargin{10},1,[]),2,1);
			obj.T = repmat(reshape(varargin{11},1,[]),2,1);
			obj.My = varargin{12};
			obj.Mz = varargin{13};
			
			%plates local variables (elemental)
			obj.uint_plates = varargin{14};
			obj.mat_plates = varargin{15};
			
			E  = obj.mat_plates(obj.Tmat_plates,2);
			nu = obj.mat_plates(obj.Tmat_plates,3);
			a = 2*sqrt(sum((obj.xplates(obj.Tplates(:,2),:)-obj.xplates(obj.Tplates(:,1),:)).^2,2))/2;
			b = 2*sqrt(sum((obj.xplates(obj.Tplates(:,3),:)-obj.xplates(obj.Tplates(:,2),:)).^2,2))/2;
			
			temp1 = E./(1-nu.^2); 
			temp2 = zeros(3,3,size(obj.Tmat_plates,1)); temp2(1,1,:) = 1; temp2(2,2,:) = 1; temp2(3,3,:) = 1; temp2(1,2,:) = nu; temp2(2,1,:) = nu; temp2(3,3,:) = (1-nu)/2;
			temp3_1 = zeros(3,8,size(obj.Tmat_plates,1)); temp3_1(1,1,:) = -b; temp3_1(3,1,:) = -a; temp3_1(2,2,:) = -a; temp3_1(3,2,:) = -b; temp3_1(1,3,:) = b; temp3_1(3,4,:) = b; temp3_1(3,7,:) = a; temp3_1(2,8,:) = a;
			temp3_2 = zeros(3,8,size(obj.Tmat_plates,1)); temp3_2(1,1,:) = -b; temp3_2(3,2,:) = -b; temp3_2(1,3,:) = b; temp3_2(3,3,:) = -a; temp3_2(2,4,:) = -a; temp3_2(3,4,:) = b; temp3_2(3,5,:) = a; temp3_2(2,6,:) = a;
			temp3_3 = zeros(3,8,size(obj.Tmat_plates,1)); temp3_3(3,3,:) = -a; temp3_3(2,4,:) = -a; temp3_3(1,5,:) = b; temp3_3(3,5,:) = a; temp3_3(2,6,:) = a; temp3_3(3,6,:) = b; temp3_3(1,7,:) = -b; temp3_3(3,8,:) = -b;
			temp3_4 = zeros(3,8,size(obj.Tmat_plates,1)); temp3_4(3,1,:) = -a; temp3_4(2,2,:) = -a; temp3_4(1,5,:) = b; temp3_4(3,6,:) = b; temp3_4(1,7,:) = -b; temp3_4(3,7,:) = a; temp3_4(2,8,:) = a; temp3_4(3,8,:) = -b;
			
			sig = zeros(12,size(obj.Tmat_plates,1));
			sig(1:3,:) = permute(reshape(temp1,1,1,[]).*sum( temp2 .*  permute(sum(temp3_1.*permute(obj.uint_plates([1,2,7,8,13,14,19,20],:),[3,1,2]),2),[2,1,3]) ,2  ),[1,3,2]);
			sig(4:6,:) = permute(reshape(temp1,1,1,[]).*sum( temp2 .*  permute(sum(temp3_2.*permute(obj.uint_plates([1,2,7,8,13,14,19,20],:),[3,1,2]),2),[2,1,3]) ,2  ),[1,3,2]);
			sig(7:9,:) = permute(reshape(temp1,1,1,[]).*sum( temp2 .*  permute(sum(temp3_3.*permute(obj.uint_plates([1,2,7,8,13,14,19,20],:),[3,1,2]),2),[2,1,3]) ,2  ),[1,3,2]);
			sig(10:12,:) = permute(reshape(temp1,1,1,[]).*sum( temp2 .*  permute(sum(temp3_4.*permute(obj.uint_plates([1,2,7,8,13,14,19,20],:),[3,1,2]),2),[2,1,3]) ,2  ),[1,3,2]);
			
			obj.sig_x = sig(1:3:end,:);
			obj.sig_y = sig(2:3:end,:);
			obj.sig_xy = sig(3:3:end,:);
			obj.sig_vm = sqrt(obj.sig_x.^2 - obj.sig_x.*obj.sig_y + obj.sig_y.^2 + 3*obj.sig_xy.^2);
			obj.p = 1/3 * (obj.sig_x + obj.sig_y);			
			
			obj.compute_variable_mapping;
			obj.create_figure;
		
		end
		
		function create_figure(obj)
            %% Dimensions
                W0 = 800;
                H0 = 500;
                Wp = 160;
                s = 10;
                Wc = 20;
                Hs = 27;
                Hp = s/2+20+5*Hs;
                Wt = 50; 
			%% Create figure
				obj.fg_handle = figure();
				set(obj.fg_handle,'Name','Fuselage','Color','w','Units','pixels','Position',[50,50,W0,H0],'SizeChangedFcn',@resizeWindow);%,'Resize',@resizeWindow)
			%% Create axes
				obj.ax_handle = axes(obj.fg_handle,'Units','pixels','Position',[Wp+Wc, 3*s, W0-Wp-2*Wc, H0-5*s],'Box','off','XTick', [],'YTick', []);
				axis(obj.ax_handle,'equal','tight');
			%%	Box of selection
				obj.pnl_selection = uipanel('Title','Selection','FontSize',12,...
					 'BackgroundColor','w','Units','pixels',...
					 'Position',[s, H0-Hp-s, Wp-s, Hp]);
				 
				obj.chb_stringer = uicontrol(obj.fg_handle,'Style','checkbox','String','Stringers',...
					'Units','pixels','Position',[2*s, H0-Hp-s/2+4*Hs, Wp-3*s, Hs],...
					'BackgroundColor','w','Value',1,'Callback',@(src,event)obj.plot_stringer);
				
				obj.chb_frame = uicontrol(obj.fg_handle,'Style','checkbox','String','Frames',...
					'Units','pixels','Position',[2*s, H0-Hp-s/2+3*Hs, Wp-3*s, Hs],...
					'BackgroundColor','w','Value',1,'Callback',@(src,event)obj.plot_frame);
				
				obj.chb_floor_spars = uicontrol(obj.fg_handle,'Style','checkbox','String','Floor reinforcements',...
					'Units','pixels','Position',[2*s, H0-Hp-s/2+2*Hs, Wp-3*s, Hs],...
					'BackgroundColor','w','Value',1,'Callback',@(src,event)obj.plot_floor_spars);
				
				obj.chb_floor_plates = uicontrol(obj.fg_handle,'Style','checkbox','String','Floor plates',...
					'Units','pixels','Position',[2*s, H0-Hp-s/2+1*Hs, Wp-3*s, Hs],...
					'BackgroundColor','w','Value',1,'Callback',@(src,event)obj.plot_floor_plates);
				
				obj.chb_skin_plates = uicontrol(obj.fg_handle,'Style','checkbox','String','Skin plates',...
					'Units','pixels','Position',[2*s, H0-Hp-s/2, Wp-3*s, Hs],...
					'BackgroundColor','w','Value',1,'Callback',@(src,event)obj.plot_skin_plates);
			%% Scale
				obj.txt_scale = uicontrol(obj.fg_handle,'Style','text',...
					'Units','pixels','Position',[s, H0-Hp-Hs-2.5*s, Wt, Hs],...
					'String','Scale:','HorizontalAlignment','left',...
					'BackgroundColor','w');

				obj.txtedit_scale = uicontrol(obj.fg_handle,'Style','edit',...
					'Units','pixels','Position',[s+Wt+s, H0-Hp-Hs-2*s, Wp-2*s-Wt, Hs],...
					'String','1',...
					'Callback',@(src,event)obj.update_scale);
			%% Symmetry
				obj.chb_sym = uicontrol(obj.fg_handle,'Style','checkbox','String','Symmetry',...
					'Units','pixels','Position',[s, H0-Hp-2*Hs-3*s, Wp-s, Hs],...
					'BackgroundColor','w','Callback',@(src,event)obj.create_symmetry);
			%% Color
				obj.dd_color = uicontrol(obj.fg_handle,'Style','popupmenu',...
					'Units','pixels','Position',[s, H0-Hp-3*Hs-4*s, Wp-s, Hs],...
					'String',{'None','Displacements','Displacement X','Displacement Y','Displacement Z',...
							'Rotation X','Rotation Y','Rotation Z','Axial force','Shear force Y','Shear force Z',...
							'Torsion','Bending moment Y','Bending moment Z','Stress X','Stress Y','Stress XY','Stress VM','Hydrostatic pressure'},'Value',1,...
					'Callback',@(src,event)obj.change_color);
			%% Print Fuselaje
				% Initial plot
				sz = 1;
				Ax = 0;
				ori = [min(obj.xbeams(:,1)),0,0]-Ax;
				hold on;
				plot3(ori(1)+[0,sz],ori(2)+[0,0],ori(3)+[0,0],'r','linewidth',1.5); text(ori(1)+1.1*sz,ori(2),ori(3),'X','color','r');
				plot3(ori(1)+[0,0],ori(2)+[0,sz],ori(3)+[0,0],'g','linewidth',1.5); text(ori(1),ori(2)+1.1*sz,ori(3),'Y','color','g');
				plot3(ori(1)+[0,0],ori(2)+[0,0],ori(3)+[0,sz],'b','linewidth',1.5); text(ori(1),ori(2),ori(3)+1.1*sz,'Z','color','b');
			
				% Plot beams
				obj.obj_stringer = patch('Vertices',obj.coord_beams_mapping,'Faces',obj.connect_beams_mapping(obj.mat_beams_mapping==2,:),'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',1);
				obj.obj_frame = patch('Vertices',obj.coord_beams_mapping,'Faces',obj.connect_beams_mapping(obj.mat_beams_mapping==1,:),'FaceColor',[0.8 0.8 0.8],'EdgeColor','b','LineWidth',1);
				obj.obj_floor_spars = patch('Vertices',obj.coord_beams_mapping,'Faces',obj.connect_beams_mapping(obj.mat_beams_mapping==3,:),'FaceColor',[0.8 0.8 0.8],'EdgeColor','r','LineWidth',1);
				
				% Plot plates
				obj.obj_floor_plates = patch('Vertices',obj.coord_plates_mapping,'Faces',obj.connect_plates_mapping(obj.mat_plates_mapping==2,:),'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
				obj.obj_skin_plates = patch('Vertices',obj.coord_plates_mapping,'Faces',obj.connect_plates_mapping(obj.mat_plates_mapping==1,:),'FaceColor',[0.6 0.6 0.6],'EdgeColor','none');
				
				view(3)
				axis equal
				set(obj.ax_handle,'xcolor','none','ycolor','none','zcolor','none');

            %% Resize window function
            function resizeWindow(varargin)
                pos = get(obj.fg_handle,'Position');
                W = pos(3);
                H = pos(4);
                set(obj.ax_handle,'Position',[Wp+Wc, 3*s, W-Wp-2*Wc, H-5*s]);
                set(obj.pnl_selection,'Position',[s, H-Hp-s, Wp-s, Hp]);
                set(obj.chb_stringer,'Position',[2*s, H-Hp-s/2+4*Hs, Wp-3*s, Hs]);
                set(obj.chb_frame,'Position',[2*s, H-Hp-s/2+3*Hs, Wp-3*s, Hs]);
                set(obj.chb_floor_spars,'Position',[2*s, H-Hp-s/2+2*Hs, Wp-3*s, Hs]);
                set(obj.chb_floor_plates,'Position',[2*s, H-Hp-s/2+1*Hs, Wp-3*s, Hs]);
                set(obj.chb_skin_plates,'Position',[2*s, H-Hp-s/2, Wp-3*s, Hs]);
                set(obj.txt_scale,'Position',[s, H-Hp-Hs-2.5*s, Wt, Hs]);
                set(obj.txtedit_scale,'Position',[s+Wt+s, H-Hp-Hs-2*s, Wp-2*s-Wt, Hs]);
                set(obj.chb_sym,'Position',[s, H-Hp-2*Hs-3*s, Wp-s, Hs]);
                set(obj.dd_color,'Position',[s, H-Hp-3*Hs-4*s, Wp-s, Hs]);
            end


		end
		
		function plot_stringer(obj)
			
			switch get(obj.chb_stringer,'Value') 
				case 1
					set(obj.obj_stringer,'Visible','on');
				case 0
					set(obj.obj_stringer,'Visible','off');
			end
			
		end
		
		function plot_frame(obj)
			
			switch get(obj.chb_frame,'Value') 
				case 1
					set(obj.obj_frame,'Visible','on');
				case 0
					set(obj.obj_frame,'Visible','off');
			end
			
		end
		
		function plot_floor_spars(obj)
			
			switch get(obj.chb_floor_spars,'Value') 
				case 1
					set(obj.obj_floor_spars,'Visible','on');
				case 0
					set(obj.obj_floor_spars,'Visible','off');
			end
			
		end
		
		function plot_floor_plates(obj)
			
			switch get(obj.chb_floor_plates,'Value') 
				case 1
					set(obj.obj_floor_plates,'Visible','on');
				case 0
					set(obj.obj_floor_plates,'Visible','off');
			end
			
		end
		
		function plot_skin_plates(obj)
			
			switch get(obj.chb_skin_plates,'Value') 
				case 1
					set(obj.obj_skin_plates,'Visible','on');
				case 0
					set(obj.obj_skin_plates,'Visible','off');
			end
			
		end
		
		function compute_variable_mapping(obj)
            % function that updates the mapping coordinates, connectivities and
            % displacement for the corresponding t_ref (t_ref for the initial
            % iteration is equal to 1)
            
            obj.coord_beams_mapping = obj.xbeams;
			obj.coord_plates_mapping = obj.xplates;
			obj.connect_beams_mapping = obj.Tbeams;
			obj.connect_plates_mapping = obj.Tplates;
			obj.mat_beams_mapping = obj.Tmat_beams;
			obj.mat_plates_mapping = obj.Tmat_plates;
			obj.U_beams_mapping = reshape(obj.ubeams,3,[])';
			obj.U_plates_mapping = reshape(obj.uplates,3,[])';
		end
				
		function update_scale(obj)
			
			scl = str2double(get(obj.txtedit_scale,'String'));
            
            if ~ischar(abs(scl))
                obj.scale = scl;
                obj.update_coord(obj.coord_beams_mapping+obj.scale*obj.U_beams_mapping,obj.coord_plates_mapping+obj.scale*obj.U_plates_mapping);
            else
                set(obj.txtedit_scale,'String',num2str(obj.scale));
            end
			
		end
		
		function update_coord(obj,coord_beams_new,coord_plates_new)
			
			set(obj.obj_stringer,'Vertices',coord_beams_new);
			set(obj.obj_frame,'Vertices',coord_beams_new);
			set(obj.obj_floor_spars,'Vertices',coord_beams_new);
			set(obj.obj_floor_plates,'Vertices',coord_plates_new);
			set(obj.obj_skin_plates,'Vertices',coord_plates_new);
			
		end
		
		function update_connect(obj,connect_new_beams,connect_new_plates,mat_new_beams,mat_new_plates)
            
			set(obj.obj_stringer,'Faces',connect_new_beams(mat_new_beams==2,:));
			set(obj.obj_frame,'Faces',connect_new_beams(mat_new_beams==1,:));
			set(obj.obj_floor_spars,'Faces',connect_new_beams(mat_new_beams==3,:));
			set(obj.obj_floor_plates,'Faces',connect_new_plates(mat_new_plates==2,:));
			set(obj.obj_skin_plates,'Faces',connect_new_plates(mat_new_plates==1,:));
			
        end
        
        function update_color(obj,color_beams_new,color_plates_new)
            
			switch get(obj.dd_color,'Value')
				case 1
					set(obj.obj_stringer,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k');
					set(obj.obj_frame,'FaceColor',[0.8 0.8 0.8],'EdgeColor','b');
					set(obj.obj_floor_spars,'FaceColor',[0.8 0.8 0.8],'EdgeColor','r');
					set(obj.obj_floor_plates,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
					set(obj.obj_skin_plates,'FaceColor',[0.6 0.6 0.6],'EdgeColor','none');
				case 2
					set(obj.obj_stringer,'FaceColor','interp','EdgeColor','interp','FaceVertexCData',color_beams_new);
					set(obj.obj_frame,'FaceColor','interp','EdgeColor','interp','FaceVertexCData',color_beams_new);
					set(obj.obj_floor_spars,'FaceColor','interp','EdgeColor','interp','FaceVertexCData',color_beams_new);
					set(obj.obj_floor_plates,'FaceColor','interp','FaceVertexCData',color_plates_new);
					set(obj.obj_skin_plates,'FaceColor','interp','FaceVertexCData',color_plates_new);
				case  {3,4,5,6,7,8,9,10,11,12,13,14}
					set(obj.obj_stringer,'FaceColor','interp','EdgeColor','interp','FaceVertexCData',color_beams_new);
					set(obj.obj_frame,'FaceColor','interp','EdgeColor','interp','FaceVertexCData',color_beams_new);
					set(obj.obj_floor_spars,'FaceColor','interp','EdgeColor','interp','FaceVertexCData',color_beams_new);
					
					set(obj.obj_floor_plates,'FaceColor',[0.8 0.8 0.8],'EdgeColor','none');
					set(obj.obj_skin_plates,'FaceColor',[0.6 0.6 0.6],'EdgeColor','none');
				case {15,16,17,18,19}
					set(obj.obj_stringer,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k');
					set(obj.obj_frame,'FaceColor',[0.8 0.8 0.8],'EdgeColor','b');
					set(obj.obj_floor_spars,'FaceColor',[0.8 0.8 0.8],'EdgeColor','r');
					
					set(obj.obj_floor_plates,'FaceColor','interp','FaceVertexCData',color_plates_new);
					set(obj.obj_skin_plates,'FaceColor','interp','FaceVertexCData',color_plates_new);
			end

        end
		
		function create_symmetry(obj)
			
			if get(obj.chb_sym,'Value')
				n_beams_coord = size(obj.coord_beams_mapping,1); n_plates_coord = size(obj.coord_plates_mapping,1);
				obj.coord_beams_mapping = [obj.coord_beams_mapping; obj.coord_beams_mapping(:,1) 2*min(obj.coord_beams_mapping(:,2))-obj.coord_beams_mapping(:,2) obj.coord_beams_mapping(:,3)];
				obj.coord_plates_mapping = [obj.coord_plates_mapping; obj.coord_plates_mapping(:,1) 2*min(obj.coord_plates_mapping(:,2))-obj.coord_plates_mapping(:,2) obj.coord_plates_mapping(:,3)];
				obj.connect_beams_mapping = [obj.connect_beams_mapping; n_beams_coord+obj.connect_beams_mapping];
				obj.connect_plates_mapping = [obj.connect_plates_mapping; n_plates_coord+obj.connect_plates_mapping];
				obj.mat_beams_mapping = [obj.mat_beams_mapping;obj.mat_beams_mapping];
				obj.mat_plates_mapping = [obj.mat_plates_mapping;obj.mat_plates_mapping];
				
				obj.U_beams_mapping = [obj.U_beams_mapping; obj.U_beams_mapping(:,1) -obj.U_beams_mapping(:,2) obj.U_beams_mapping(:,3)];
				obj.U_plates_mapping = [obj.U_plates_mapping; obj.U_plates_mapping(:,1) -obj.U_plates_mapping(:,2) obj.U_plates_mapping(:,3)];
				
				obj.update_connect(obj.connect_beams_mapping,obj.connect_plates_mapping,obj.mat_beams_mapping,obj.mat_plates_mapping);
				obj.update_scale;
				obj.change_color;
			else
				obj.compute_variable_mapping;
				
				obj.update_connect(obj.connect_beams_mapping,obj.connect_plates_mapping,obj.mat_beams_mapping,obj.mat_plates_mapping);
				obj.update_scale;
				obj.change_color;
			end
			
		end
		
		function change_color(obj)
			
			n_sym = get(obj.chb_sym,'Value');
			
			switch get(obj.dd_color,'Value')
				case 2 %'Displacement'
					u_rshp = reshape(obj.ubeams,3,[])'; Color_beams_mapping = sqrt(sum(u_rshp.^2,2));
					u_rshp = reshape(obj.uplates,3,[])'; Color_plates_mapping = sqrt(sum(u_rshp.^2,2));
					
					Color_beams_mapping = [Color_beams_mapping;repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [Color_plates_mapping;repmat(Color_plates_mapping,n_sym,1)];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 3 %'Displacementx'
					Color_beams_mapping = reshape(obj.uint_beams([1,7],:),[],1);
										
					Color_beams_mapping = [Color_beams_mapping;repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 4 %'Displacementy'
					Color_beams_mapping = reshape(obj.uint_beams([2,8],:),[],1);
										
					Color_beams_mapping = [Color_beams_mapping;-repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 5 %'Displacementz'
					Color_beams_mapping = reshape(obj.uint_beams([3,9],:),[],1);
										
					Color_beams_mapping = [Color_beams_mapping;repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 6 %'Rotx'
					Color_beams_mapping = reshape(obj.uint_beams([4,10],:),[],1);
										
					Color_beams_mapping = [Color_beams_mapping;-repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 7 %'Roty'
					Color_beams_mapping = reshape(obj.uint_beams([5,11],:),[],1);
										
					Color_beams_mapping = [Color_beams_mapping;repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 8 %'Rotz'
					Color_beams_mapping = reshape(obj.uint_beams([6,12],:),[],1);
										
					Color_beams_mapping = [Color_beams_mapping;-repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 9 %'AxialForce'
					Color_beams_mapping = reshape(obj.N,[],1);
										
					Color_beams_mapping = [Color_beams_mapping;repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 10 %'Sheary'
					Color_beams_mapping = reshape(obj.Qy,[],1);
										
					Color_beams_mapping = [Color_beams_mapping;-repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 11 %'Shearz'
					Color_beams_mapping = reshape(obj.Qz,[],1);
										
					Color_beams_mapping = [Color_beams_mapping;repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 12 %'Torsion'
					Color_beams_mapping = reshape(obj.T,[],1);
										
					Color_beams_mapping = [Color_beams_mapping;-repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 13 %'Bendingy'
					Color_beams_mapping = reshape(obj.My,[],1);
										
					Color_beams_mapping = [Color_beams_mapping;repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 14 %'Bendingz'
					Color_beams_mapping = reshape(obj.Mz,[],1);
										
					Color_beams_mapping = [Color_beams_mapping;-repmat(Color_beams_mapping,n_sym,1)];
					Color_plates_mapping = [];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 15 %'sigx'
					Color_plates_mapping = reshape(obj.sig_x,[],1);
					
					Color_beams_mapping = [];
					Color_plates_mapping = [Color_plates_mapping;repmat(Color_plates_mapping,n_sym,1)];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 16 %'sigy'
					Color_plates_mapping = reshape(obj.sig_y,[],1);
					
					Color_beams_mapping = [];
					Color_plates_mapping = [Color_plates_mapping;repmat(Color_plates_mapping,n_sym,1)];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 17 %'sigxy'
					Color_plates_mapping = reshape(obj.sig_xy,[],1);
					
					Color_beams_mapping = [];
					Color_plates_mapping = [Color_plates_mapping;repmat(Color_plates_mapping,n_sym,1)];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 18 %'sigvm'
					Color_plates_mapping = reshape(obj.sig_vm,[],1);
					
					Color_beams_mapping = [];
					Color_plates_mapping = [Color_plates_mapping;repmat(Color_plates_mapping,n_sym,1)];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				case 19 %'p'
					Color_plates_mapping = reshape(obj.p,[],1);
					
					Color_beams_mapping = [];
					Color_plates_mapping = [Color_plates_mapping;repmat(Color_plates_mapping,n_sym,1)];
					
					caxis(obj.ax_handle,[min([Color_beams_mapping(:);Color_plates_mapping(:)]) max([Color_beams_mapping(:);Color_plates_mapping(:)])]); colormap(obj.ax_handle,jet(50)); colorbar;
				otherwise
					Color_beams_mapping = [];
					Color_plates_mapping = [];
					colorbar('off');
			end
			
			obj.update_color(Color_beams_mapping,Color_plates_mapping);
			
		end
		
	end
end
