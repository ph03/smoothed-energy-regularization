%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Smoothed Quadratic Energies on Meshes
%%  J. Martinez Esturo, C. Rössl, and H. Theisel
%%
%%  ACM Transactions on Graphics 2014
%%
%%  Copyright J. Martinez Esturo 2014 (MIT License)
%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef DeformViewControl2D < handle
  %DEFORMVIEWCONTROL2D

  properties
    deform
    deforms = {}
    active_deform_idx = []

    figure_h = []
    surf_h   = []
    texsurf_h = []

    uitext1_h    = []
    uitext2_h    = []
    uibeta_h     = []
    uiconstr_h   = []
    uismoother_h = []
    uideform_h   = []

    drag_lock = false

    show_energy = true
    show_energy_max = 8

    show_detailed_text = 1

    csolver_softweight = 1
  end

  methods
    function obj = DeformViewControl2D(deform)
      %% Constructor
      obj.deform = deform;

      obj.deforms{1} = deform;
      obj.active_deform_idx = 1;
    end

    function addDeform(obj,deform)
      %% call this before show()..
      assert(isempty(obj.figure_h));

      obj.deforms{length(obj.deforms)+1} = deform;
    end

    function show(obj,scale)
      if nargin < 2,
        scale = 2;
      end

      obj.figure_h = figure;

      obj.drawMesh();

      colormap(jet(256));
      caxis([0, obj.show_energy_max]);
      colorbar;

      %% center axis on mesh
      p = obj.deform.mesh.p;
      minx  =  min(p(1,:));  maxx =  max(p(1,:));
      meanx = mean(p(1,:)); meany = mean(p(2,:));

      d = max([scale*(meanx - minx),scale*(maxx - meanx)]);

      axis_arr = [meanx - d meanx + d ...
                  meany - d meany + d];
      axis(axis_arr);
      axis square;

      uicontrol('Style', 'pushbutton', 'String', 'initHandles',...
                'Position', [20 0 45 20],...
                'Callback', @(~,~)initHandles(obj,[]));

      uicontrol('Style', 'pushbutton', 'String', 'resetMesh',...
                'Position', [65 0 45 20],...
                'Callback', @(~,~)resetMesh(obj));

      obj.uitext1_h = text('Interpreter','latex',...
                           'String', '',...
                           'Position', [meanx - 0.9*d meany + 0.9*d], ...
                           'HorizontalAlignment', 'left',...
                           'FontSize', 25);
      obj.uitext2_h = text('Interpreter','latex',...
                           'String', '',...
                           'Position', [meanx - 0.9*d meany + 0.8*d], ...
                           'HorizontalAlignment', 'left',...
                           'FontSize', 25);

      uicontrol('Style','text',...
                'Position',[120 20 100 20],...
                'String', 'beta');

      obj.uibeta_h = uicontrol('Style', 'slider', 'String', 'beta',...
                               'Min',0,'Max',1,'Value',obj.deform.beta,...
                               'Position', [120 0 100 20],...
                               'Callback', @(hObj,event)updateBetaSliderEvent(obj,hObj,event));

      deformstr = obj.deform.name;
      for i=2:length(obj.deforms),
        deformstr = strcat(deformstr,'|',obj.deforms{i}.name);
      end

      obj.uideform_h = uicontrol('Style', 'popup',...
                                 'String', deformstr,...
                                 'Value', obj.active_deform_idx,...
                                 'Position', [20 20 90 20],...
                                 'Callback', @(hObj,event)updateDeformEvent(obj,hObj,event));

      obj.uiconstr_h = uicontrol('Style', 'popup',...
                                 'String', 'CholConstrSolver|SoftConstrSolver|LagrConstrSolver',...
                                 'Value', 1,...
                                 'Position', [230 0 145 20],...
                                 'Callback', @(hObj,event)updateConstrSolverEvent(obj,hObj,event));

      uicontrol('Style','text',...
                'Position',[230 20 60 20],...
                'String', 'soft w');

      uicontrol('Style', 'slider', 'String', 'soft_constr_weight',...
                'Min',1e-4,'Max',1,'SliderStep',[0.1 0.6],'Value',obj.csolver_softweight,...
                'Position', [295 20 80 20],...
                'Callback', @(hObj,event)updateConstrSolverSoftWeightEvent(obj,hObj,event));

      uicontrol('Style', 'checkbox', 'String', 'show_energy',...
                'Value', obj.show_energy,...
                'Position', [380 0 100 20],...
                'Callback', @(hObj,event)updateShowEnergyEvent(obj,hObj,event));

      uicontrol('Style', 'slider', 'String', 'show_energy_max',...
                'Min',0,'Max',10,'Value', obj.show_energy_max,...
                'Position', [380 20 100 20],...
                'Callback', @(hObj,event)updateShowEnergyMaxEvent(obj,hObj,event));

      uicontrol('Style','text',...
                'Position',[490 20 150 20],...
                'String', 'Regularization');

      obj.uismoother_h = uicontrol('Style', 'popup',...
                                   'String', 'SmoothedEnergy|TikhonovGradient',...
                                   'Position', [490 0 150 20],...
                                   'Callback', @(hObj,event)updateSmootherEvent(obj,hObj,event));

      uicontrol('Style', 'checkbox', 'String', 'show_detailed_text',...
                'Value', obj.show_detailed_text,...
                'Position', [650 0 150 20],...
                'Callback', @(hObj,event)updateShowDetailedTextEven(obj,hObj,event));

      set(obj.figure_h, ...
          'KeyPressFcn',           @(src,event)keyEvent(obj,src,event));

      obj.updateText();
    end

    function resetMesh(obj)
      mesh = obj.deform.mesh;
      mesh.p = mesh.p_0;

      obj.resetHandles();

      obj.drawMesh();
    end

    function resetHandles(obj)
      figure(obj.figure_h);

      handles = getappdata(gca,'handles');
      if ~isempty(handles),
        delete(handles);
        setappdata(gca,'handles',[]);
      end
    end

    function updateDeformEvent(obj,hObj,~)
      obj.setActiveDeform(get(hObj,'Value'));
    end

    function updateBetaSliderEvent(obj,hObj,~)
      obj.updateBeta(get(hObj,'Value'));
    end

    function updateConstrSolverEvent(obj,hObj,~)
      obj.updateConstrSolver(get(hObj,'Value'));
    end

    function updateShowEnergyEvent(obj,hObj,~)
      obj.show_energy = get(hObj,'Value');
    end

    function updateShowEnergyMaxEvent(obj,hObj,~)
      obj.show_energy_max = get(hObj,'Value');
      caxis([0, obj.show_energy_max]);
    end

    function updateSmootherEvent(obj,hObj,~)
      obj.updateSmoother(get(hObj,'Value'));
    end

    function updateConstrSolverSoftWeightEvent(obj,hObj,~)
      obj.updateConstrSolverSoftWeight(get(hObj,'Value'));
    end

    function updateShowDetailedTextEven(obj,hObj,~)
      obj.show_detailed_text = get(hObj,'Value');
      obj.updateText();
    end

    function setActiveDeform(obj,deform_idx)
      assert(deform_idx <= length(obj.deforms));

      if deform_idx == obj.active_deform_idx,
        return;
      end

      %% save old deform state
      hidxs     = obj.deform.hidxs;
      beta      = obj.deform.beta;
      csolverf  = obj.deform.csolverf;
      smootherf = obj.deform.smootherf;

      %% switch deform and update deformation
      obj.active_deform_idx = deform_idx;

      obj.deform = obj.deforms{deform_idx};

      obj.deform.beta      = beta;
      obj.deform.csolverf  = csolverf;
      obj.deform.smootherf = smootherf;

      %% update deformation
      obj.deform.init(hidxs);
      obj.deformConstr();

      obj.updateText();
    end

    function beta = updateBeta(obj,beta)
      if nargin < 2,
        beta = obj.deform.beta;
        return;
      end

      obj.deform.beta = beta;

      obj.updateText();

      %% update deformation
      obj.deformConstr();
    end

    function updateConstrSolver(obj,csolver)
      switch csolver,
        case 1,
          obj.deform.csolverf = @(A,cidxs,b)CholConstrSolver(A,cidxs,b);
        case 2,
          obj.deform.csolverf = @(A,cidxs,b)SoftConstrSolver(A,cidxs,b,obj.csolver_softweight);
        case 3,
          obj.deform.csolverf = @(A,cidxs,b)LagrConstrSolver(A,cidxs,b);
      end

      obj.updateText();

      %% update deformation
      obj.deformConstr();
    end

    function updateConstrSolverSoftWeight(obj,csolver_softweight)
      obj.csolver_softweight = csolver_softweight;
      obj.updateConstrSolver(get(obj.uiconstr_h,'Value'));
    end

    function updateSmoother(obj,smoother)
      switch smoother,
        case 1,
          obj.deform.smootherf = @(mesh,beta, en, d,c) ...
            TriSystemSmoothEnergy(mesh,beta, en, d,c);
        case 2,
          obj.deform.smootherf = @(mesh,beta, en, d,c) ...
            TriSystemTikhonovGradientEnergy(mesh,beta, en, d,c);
      end

      obj.updateText();

      %% update deformation
      obj.deformConstr();
    end

    function hidxs = initHandles(obj,hidxs)
      if(ishandle(obj.figure_h))
        figure(obj.figure_h); hold on;
        obj.resetHandles();

        mesh = obj.deform.mesh;

        if nargin < 2 || length(hidxs) < 2,
          %% get user input
          [cx,cy] = getpts();

          hidxs = mesh.closestVertices([cx';cy']);
          hidxs = unique(hidxs);
        end

        vcoords = mesh.p(:,hidxs);

        %% create draggable handles
        handles = zeros(1,length(hidxs));

        for i=1:length(hidxs),
          h = plot(vcoords(1,i),vcoords(2,i),'r');

          draggable(h,@(h)dragEvent(obj,h),'endfcn',@(h)dragEndEvent(obj,h));

          handles(i) = h;
        end

        set(handles,'Marker','o','MarkerSize',10,'MarkerFaceColor','b');

        setappdata(gca,'handles',handles);

        hold off;
      else
        assert(nargin > 1);
      end

      %% init deform implementation
      obj.deform.init(hidxs);
    end
  end

  methods (Access = protected)
    function drawMesh(obj)
      %% (re) draw mesh
      figure(obj.figure_h); hold on;

      if(ishandle(obj.surf_h))
        delete(obj.surf_h);
      end

      obj.surf_h = obj.deform.mesh.draw();

      view(2);
      axis equal;
      axis manual;

      hold off;
    end

    function dragEvent(obj,~)
      if(~obj.drag_lock),
        obj.drag_lock = true;
      else
        return;
      end

      obj.deformConstr();

      obj.drag_lock = false;
    end

    function dragEndEvent(~,~)

    end

    function updateText(obj)
      if obj.show_detailed_text,
        switch get(obj.uismoother_h, 'Value'),
          case 1,
            smoother = 'Smoothed Energy Regularization (our)';
            color = [62 167 16] ./ 255;
          case 2,
            smoother = 'Tikhonov Regularization';
            color = [1 0 1];
        end

        if obj.deform.beta > 0,
          latex_string = sprintf('$\\beta = %.2f$, %s',...
                                 obj.deform.beta,...
                                 smoother);
        else
          latex_string = sprintf('$\\beta = %.2f$, No Regularization', ...
                                 obj.deform.beta);
          color = [0 0 0];
        end

        switch get(obj.uiconstr_h, 'Value'),
          case 1,
            constr = 'Hard Constraints';
          case 2,
            constr = 'Soft Constraints';
          case 3,
            constr = 'Hard Constraints';
        end

        set(obj.uitext1_h, 'String', sprintf('%s, %s', obj.deform.name, constr));
        set(obj.uitext2_h, 'String', latex_string);
        set(obj.uitext2_h, 'Color', color);

      else
        latex_string = sprintf('$\\beta = %.2f$',obj.deform.beta);

        set(obj.uitext1_h, 'String', latex_string);
        set(obj.uitext2_h, 'String', ' ');
      end
    end

    function keyEvent(obj,~,event)
      key = event.Character;

      if (key == 'r'),
        if(~ishandle(obj.surf_h)), return; end;

        %% reset mesh

        % reset handles
        p_0 = obj.deform.mesh.p_0;
        hidxs = obj.deform.hidxs;

        %reset beta
        if(1),
          obj.deform.beta = 0;
          set(obj.uibeta_h,'Value',0);
        end

        if(~isempty(obj.texsurf_h))
          delete(obj.texsurf_h)
          obj.texsurf_h = [];
        end

        set(obj.surf_h,'Visible', 'on');

        obj.updateText();

        obj.deformConstr(p_0(:,hidxs));
      end

      if (key == 't'),
        if(~ishandle(obj.surf_h)), return; end;

        if(~isempty(obj.deform.mesh.texture))
          if(~isempty(obj.texsurf_h))
            delete(obj.texsurf_h);
          end

          obj.deform.mesh.setForceTexture(true);
          obj.texsurf_h = obj.deform.mesh.draw();
          obj.deform.mesh.setForceTexture(false);

          set(obj.surf_h,'Visible', 'off');
        end
      end
    end
  end

  methods (Access = public)
    function deformConstr(obj,hcoords)
      if nargin < 2
        %% get handle coordinates from figure
        figure(obj.figure_h);
        handles = getappdata(gca,'handles');
        xdata = cell2mat(get(handles,'xdata'));
        ydata = cell2mat(get(handles,'ydata'));

        hcoords = [xdata';ydata'];
      else
        %% set handle coordinates
        if(ishandle(obj.figure_h))
          handles = getappdata(gca,'handles');
          for i=1:length(handles),
            set(handles(i),'xdata',hcoords(1,i));
            set(handles(i),'ydata',hcoords(2,i));
          end
        end
      end

      %% call deform implementation
      itn = 1; converged = false;
      while ~converged,
        [converged,u] = obj.deform.deform(hcoords, itn);
        itn = itn + 1;
      end

      if(ishandle(obj.surf_h)),
          set(obj.surf_h,'Vertices',obj.deform.mesh.p');
        if(~isempty(obj.texsurf_h))
          delete(obj.texsurf_h);
          obj.texsurf_h = [];
        end
      end

      %% show energy
      if(obj.show_energy),
        terrEnergy = obj.deform.se.evalTriangleErrors(u);

        set(obj.surf_h,'CData',terrEnergy);
      end
    end
  end
end
