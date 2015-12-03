function CO2B(action)% CO2B This GUI object is designed to be used %        in the Computational Optoelectronics course.%        It stores all the variables which are%        used globally in the second exercise (CO2B).%%        There is one action button:%        RECALCULATE: the values are accepted and the%             intensity distribution on the screen calculated.%             This calculation is done by one of an externally%             GC.%%        On entry, the optional variable ACTION determines%        possible actions to be taken.%--------------------------------------------------------------------%         Test the entry information is in order%--------------------------------------------------------------------	   set(0,'ShowHiddenHandles','on');% Set up default ACTION if it wasn't explicitly passed.if nargin==0,   action = 'OpenReopen'; end;   % ------------------------------------------------------------------% ------------------------------------------------------------------% PART 1:    INITIALIZE THE GRAPHICAL USER INTERFACE% ------------------------------------------------------------------% ------------------------------------------------------------------% This is where the input variables are specified.% You specify the number of variables, their names, their units% and their default values.  The Driver does all the rest.   numInVars = 5;   inputVars = ...	  {'lambda1';'lambda2';'slitWdth';'slitSepn';'numSlits'};   inputUnits = ...	  {   'm'  ;   'm'    ;   'm'    ;    'm'   ; ''};   inputValus = ...	  {'600e-9'; '600e-9' ; '3.0e-7' ; '1.5e-6' ; '20'};  % The input variables are made global in the local workspace % and also in the command workspace.   globalStatement = ['global'];   for i=1:numInVars,      globalStatement = [globalStatement, ' ', char(inputVars(i))];   end;   globalStatement = [globalStatement,';'];   eval(globalStatement);   evalin('base', globalStatement);	   	    % Some extra constants needed,   tealColor  = [0.7,0.7,0.50];   midGray    = [0.50,0.50,0.50];   panelName  = 'CO2B Input panel';   panelTag   = 'CO2B InputPanel';   angleRange = 75;   % limits of the angular detectors % Set up default graphing colours   set(0,'defaultfigurecolor','black');   % Black backgrounds for all graphs   set(0,'defaultaxescolor'  ,'black');   set(0,'defaultaxesXColor' ,'white');   % White printing for all axes   set(0,'defaultaxesYColor' ,'white');if strcmp(action,'OpenReopen')    % ---------------------------------------------------   %        Reopen the figure or make a new one   % ---------------------------------------------------   h_figs   = get(0,'children');   Hf_input = findobj(h_figs,'flat','tag',panelTag);   % Note: the figure used is tagged with the string panelTag.   % If it's not open, the figure will have to be initialized,   % and all the global variables set to their default values.    if length(Hf_input)==0, % If the panel is not not already open perform the following:          % ---------------------------------------------------      %      Give the variables their default values      % ---------------------------------------------------      for i=1:numInVars,         itemStr = char(inputVars(i));      % variable names	     valuStr = char(inputValus(i));     % default values         eval([itemStr '=' valuStr ';']);      end;       					      % ---------------------------------------------------    %         Calculate all display dimensions    % ---------------------------------------------------	% First establish the height and width of the screen	  defaultUnits = get(0,'units');	  set(0,'units','centimeters');     % Our positions are in cm	  scrnSize   = get(0,'screensize'); 	  set(0,'units', defaultUnits);     % Restore root figure	  scrnWidth   = scrnSize(1,3);	  scrnHeight  = scrnSize(1,4);	     % Set dimensions of the input panel	  panelWidth  = 10.0;	  varblHeight = 1.3;	  panelHeight = 4.0 + numInVars*varblHeight;	  panelRect   = [panelWidth panelHeight];					  panelTop    = [scrnWidth scrnHeight] * 0.9;	  panelCorner = panelTop - panelRect;	      % Set dimensions of the frame in the panel	  frameMargin = 0.6;	  frameHeight = (numInVars+0.5) * varblHeight;	  frameWidth  = panelWidth - 2*frameMargin;	  frameLeft   = (panelWidth-frameWidth)/2;	  frameRect   = [frameWidth frameHeight];	  frameTop    = panelHeight - frameMargin;	  frameCorner = [frameLeft frameTop-frameHeight];	      % Set dimensions of the 'recalculate' button	  buttnWidth  = 5.0;	  buttnHeight = frameCorner(1,2) - 2*frameMargin;	  buttnLeft   = (panelWidth-buttnWidth)/2;	  buttnRect   = [buttnWidth buttnHeight];	  buttnCorner = [buttnLeft frameMargin];    % Set dimensions of the frame title	  textHeight  = 0.5;	  textOffset  = 0.8;	  titleWidth  = 3.5;	  titleLeft   = frameLeft + textOffset/2;	  titleBottom = frameTop - textHeight/2;	  titleRect   = [titleWidth textHeight];	      % Set dimensions of the text boxes	  namesWidth   = 2.2;	  namesLeft    = frameLeft + textOffset;	  namesTop     = frameTop;	  namesRect    = [namesWidth textHeight];	  	  unitWidth    = 1.0;	  unitLeft     = frameLeft + frameWidth - textOffset - unitWidth;	  unitRect     = [unitWidth textHeight];	  valueWidth   = frameWidth - namesWidth-unitWidth-2*textOffset;	  valueLeft    = (namesLeft+namesWidth + unitLeft-valueWidth)/2;	  valueRect    = [valueWidth 1.4*textHeight];							      % ---------------------------------------------------    %         Open panel and record handle    % ---------------------------------------------------      Hf_input = figure('color',tealColor,...                          'units','centimeters', ...                          'position',[panelCorner panelRect],...                          'Name',panelName,...						  'HandleVisibility','off',...				          'Tag',panelTag,...				          'NumberTitle','off');      handles(1) = Hf_input;  	  	  			        % ---------------------------------------------------      %    Set up objects to determine/display variables      % ---------------------------------------------------      % ------------------------------------      %         Input variables      % ------------------------------------      % The frame for input variables	  handInpt1 = 2;      handles(handInpt1) = uicontrol(Hf_input,...	                        'style','frame',...                            'backgroundcolor',tealColor,...                            'units','centimeters',...                            'position',[frameCorner frameRect]);							% Note: this handle not kept!      % name for the frame       handles(handInpt1) = uicontrol(Hf_input,...	                        'style','text',...                            'units','centimeters',...                            'backgroundcolor',tealColor,...    	                    'position',[titleLeft titleBottom titleWidth textHeight],...                            'string','Input variables');							% Note: this handle not kept!	  for i=1:numInVars,	     handInpt = handInpt1+i-1;           % number of current handle	     bottom   = namesTop-varblHeight*i;		 grphVar  = char(inputVars(i));		 valuVar  = inputValus(i);		 unit     = char(inputUnits(i));      % name for the variable          handles(handInpt) = uicontrol(Hf_input,...	                        'style','text',...                            'units','centimeters',...                            'backgroundcolor',tealColor,...    	                    'position',[namesLeft bottom namesRect],...                            'string',[grphVar ' =']);							% Note: this handle not kept!       % units for the variable          handles(handInpt) = uicontrol(Hf_input,...	                        'style','text',...                            'units','centimeters',...                            'backgroundcolor',tealColor,...     	                    'position',[unitLeft bottom unitRect],...                            'string',unit);  							% Note: this handle not kept!       % input for the variable         handles(handInpt) = uicontrol(Hf_input,...	                        'Style','edit',...							'Visible','on',...                            'units','centimeters',...                            'position',[valueLeft bottom-0.1 valueRect],...                            'String',valuVar);	  end; 							       % ---------------------------------------------------      %         Set up action buttons      % ---------------------------------------------------						      handActn = handInpt+1;	  % (19) Button for recalculate      handles(handActn) = uicontrol(Hf_input,...	                        'Style','pushbutton',...                            'backgroundcolor',tealColor,...                            'units','centimeters',...                            'position',[buttnCorner buttnRect],...                            'string','Recalculate');	         recalcCall = ...         ['CO2B(''Recalculate'')'];      set(handles(handActn),'callback',recalcCall);	   	         set(Hf_input,'userdata',handles); 	   	  % ------------------------------------------------------ %     If figure pre-existing, retrieve parameters % ------------------------------------------------------      else      handles = get(Hf_input,'userdata');	  figure(Hf_input);	  % bring the panel forward	     end;    % if length(Hf_input)==0   end;     % if strcmp(action,'OpenReopen')% ------------------------------------------------------% ------------------------------------------------------% PART 2:         ACCEPT THE VARIABLES % ------------------------------------------------------% ------------------------------------------------------if strcmp(action,'Recalculate')      % ---------------------------------------------------   %        Recall all the experimental variables   % ---------------------------------------------------   Hf_input   = gcf;   handles    = get(Hf_input,'userdata');   for i=1:numInVars,      itemStr = char(inputVars(i));          % variable names	  valuStr = get(handles(i+1),'String');  % variable values      eval([itemStr '=' char(valuStr) ';']);   end;          % ---------------------------------------------------   %        Reopen the graph or make a new one   % ---------------------------------------------------   h_figs    = get(0,'children');   Hf_screen = findobj(h_figs,'flat','tag','ScreenGraph');   % Note: the figure used is tagged with the name 'ScreenGraph'.   % If it's not open, the figure will have to be initialized,   % and the appropriate axes set up.    if length(Hf_screen)==0,%  If the graph is not not already open perform the following: 	  screenWidth  = 0.90;	  screenHeight = 0.40;	  screenRect   = [screenWidth screenHeight];					  screenCorner = [0.02 0.06];	      % ---------------------------------------------------    %         Open screen and record handle    % ---------------------------------------------------      Hf_screen = figure('color','black',...                         'units','normalized', ...                         'position',[screenCorner screenRect],...                         'Name','Intensity on the screen',...			             'renderer','zbuffer',...			             'Tag','ScreenGraph',...			             'NumberTitle','off');						     % Specify the properties of the intensity graph window      Ha_itot = axes('Color','black',...                     'Position',[0.080 0.400 , 0.850 0.500],...                     'XLim', [-angleRange,+angleRange],...                     'Box','on');					 	  H_xlabel = get(Ha_itot,'xlabel');	  set(H_xlabel,'String','Detector  angle   (degrees)')	  	  	  H_ylabel = get(Ha_itot,'ylabel');	  set(H_ylabel,'String','Intensity   (arbitrary  units)')	  					 	      % Specify the properties of the three intensity graphs       Hl_itot  = line('color',midGray,...					  'xdata',0,'ydata',0);      Hl_itot1 = line('color','yellow',...					  'xdata',0,'ydata',0);      Hl_itot2 = line('color','yellow',...					  'xdata',0,'ydata',0);					      % Specify the pseudocolor contour plot inside the window      Ha_gray = axes('Color','black',...                     'Position',[0.080 0.100 , 0.850 0.150],...                     'XLim', [-angleRange,+angleRange],...                     'Box','on');					        % Record the handles for future plotting	  set(Hf_screen,'userdata',...	                [Ha_itot, Hl_itot, Hl_itot1, Hl_itot2, Ha_gray]);	     else 	  figure(Hf_screen);	  % bring screen graph forward by getting handles for graphs	  handles  = get(Hf_screen,'userdata');	  Ha_itot  = handles(1,1);	  Hl_itot  = handles(1,2);	  Hl_itot1 = handles(1,3);	  Hl_itot2 = handles(1,4);	  Ha_gray  = handles(1,5);	  	  			     end; % if length(Hf_screen)==0,  % ------------------------------------------------------% ------------------------------------------------------% PART 3: CALCULATE AND PLOT THE INTENSITY DISTRIBUTIONS% ------------------------------------------------------% ------------------------------------------------------ % ------------------------------------------ % (1) CALCULATE THE INTERFERENCE PATTERNS % ------------------------------------------   % ------------------------------------------   % Specify the detector angles   % ------------------------------------------   % All the detectors lie around a circle   % They are specified by one coordinate, the angle theta   % which must be chosen so that is it sure to sample the   % intensity distribution at the centre of (possibly very   % narrow) peaks.  This is done by an internal function.   % Each wavelength needs its own vector for theta.     theta1 = DivideAngleAxis(lambda1,angleRange);     theta2 = DivideAngleAxis(lambda2,angleRange);   % Now convert from radians to degrees for plotting     degrees1 = 180/pi*theta1;     degrees2 = 180/pi*theta2;       % ------------------------------------------------------   % Calculate one or two intensity distributions    % ------------------------------------------------------   % Call external function to calculate intensity distribution   % This function must be given the name "GC"   % It requires the array of angles to be passed.     Itot1 = GC(lambda1,theta1);	 Itot2 = GC(lambda2,theta2);	 	 Imax  = max([max(Itot1), max(Itot2)]);	    % -----------------------------------------------------   %  Recalculate vectors for Itot on same scale of theta   % -----------------------------------------------------      % Choose the vector for theta with the most subdivisions   % Recalculate one of the intensities externally     if length(theta1)>length(theta2)	    degrees = degrees1;        Irel1 = Itot1/Imax;        Irel2 = GC(lambda2,theta1)/max(Itot2);     else	    degrees = degrees2;        Irel1 = GC(lambda1,theta2)/max(Itot1);	    Irel2 = Itot2/Imax;     end;       % Treat lambda1=lambda2 as a special case, and that the user   % wants to have only one frequency present.	 if lambda1~=lambda2,	    Itot = Imax * (Irel1 + Irel2);	 else	    Itot = Imax * Irel1;	 end; % ------------------------------------------ % (2) DRAW THE SIMULATED SPECTROGRAM % ------------------------------------------       % open the graphs panel and plot results      figure(Hf_screen);       % Use two Itot vectors with same vector for theta      DrawSpectrogram(Ha_gray,degrees,lambda1,lambda2,Irel1,Irel2);      hold off;         % ------------------------------------------ % (3) PLOT THE THREE INTENSITY LINE GRAPHS % ------------------------------------------   set(Ha_itot, 'YLim',[0,Imax]);    set(Hl_itot,'xdata',degrees,'ydata',Itot);   hold on;      set(Hl_itot1,'color',ColorCode(lambda1),...                'xdata',degrees1,...                'ydata',Itot1);   set(Hl_itot2,'color',ColorCode(lambda2),...                'xdata',degrees2,...                'ydata',Itot2);   zoom on;   hold off;       % -----------------------------------------------------   %   Reset configuration for further calculations   % -----------------------------------------------------   % reopen the input panel   figure(Hf_input);   end;   % if strcmp(action,'Recalculate')% ------------------------------------------------------------------%          ALL INPUT ACTIONS COMPLETE  % ------------------------------------------------------------------   set(0,'ShowHiddenHandles','off');   return;% ------------------------------------------------------------------% ******************************************************************% ******************************************************************% ------------------------------------------------------------------% ------------------------------------------------------------------%          INTERNAL FUNCTION TO DIVIDE UP THE ANGLE AXIS  % ------------------------------------------------------------------   function theta = DivideAngleAxis(lambda, range)  global slitSepn numSlits;     % All the detectors lie around a circle   % The detectors are specified by one coordinate, the angle theta   % However it is important not to miss out on any of the peaks   % which occur when alpha = pi*a*sin(lambda/a)/(lambda*numSlits)   % The number of sample points between each peak will be numSmpls     numSmpls = 20;     dAlpha = pi*slitSepn*sin(lambda/slitSepn)...	                     /(lambda*numSlits*numSmpls);     % The limits of theta are +/- range.  Convert to radians.     range = range*pi/180;     % These correspond to: alpha = +/- pi*slitSepn*sin(range)/lambda.   % This allows us to calculate how many detection points there are     numDtect = 2*fix(pi*slitSepn*sin(range)/lambda/dAlpha) + 1;	    % Theb we calculate row vectors for alpha and theta	      alpha    = dAlpha * [-(numDtect-1)/2 : 1 : +(numDtect-1)/2];     theta    = asin(lambda*alpha/(pi*slitSepn));  return;% ------------------------------------------------------------------%        FUNCTION TO DRAW IN SPECTROSCOPIC REPRESENTATION% ------------------------------------------------------------------   function DrawSpectrogram(Ha_gray,degrees,lambda1,lambda2,Irel1,Irel2);    %%%% FIRST SET UP THE SPECIAL COLOUR MAP %%%% % Select colour for graphs to match the two wavelengths   color1 = ColorCode(lambda1);   color2 = ColorCode(lambda2);    % Interpolate three intermediate colours between them   colorRange = [[0:4]'*color1 + [4:-1:0]'*color2]/4;    % Which of the five colours is "seen" at each detector   % depends on the average of the two intensities.   % If only one intensity is present, the colour  % "seen" is either color1 or color2. % Otherwise it is one of the intermediate colours.   lineColors = fix(5*Irel1./(Irel1+Irel2));    % First construct special color map in the five colours. % This color map consists of five distinct sub-maps % Each has twelve gradations of brightness in a single colour. % There are five different colours: color1 and color2 and % the three intermediate colours.   grayMap = gray(12);   thisColorMap = [grayMap(:,1)*colorRange(1,:);...                   grayMap(:,1)*colorRange(2,:);...                   grayMap(:,1)*colorRange(3,:);...                   grayMap(:,1)*colorRange(4,:);...                   grayMap(:,1)*colorRange(5,:)]; % Prepare to plot intensity as a contourplot with these colours   set(gcf,'CurrentAxes',Ha_gray);    % The brightness within each sub-map, is determined by % the intensity.  It must not be allowed to exceed the % maximim value, else the wrong sub-map  will be chosen.    Irel = Irel1+Irel2;   tooHigh = find(Irel>0.995);   Irel(tooHigh) = 0.995;        % truncate at 0.995   Icol = (Irel + lineColors)/5; % Now do the actual plotting (as a contour) % Use the color map just constructed   colormap(thisColorMap);       %%%% NOW DO THE DRAWING (AS A CONTOUR MAP) %%%% % Set up the three-dimensional grid of points   [xSc, ySc] = meshgrid(degrees,[-1, +1]);   zSc = [Icol;Icol]; % Draw the plot with appropriate plotting parameters   pcolor(xSc,ySc,zSc);   caxis([0,1]);   shading flat;   set(gca,'TickDir','out','YTick',[]);   hold on;        return;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % ------------------------------------------------------% ------------------------------------------------------% PART 4: SUBSIDIARY FUNCTIONS% ------------------------------------------------------% ------------------------------------------------------function thiscolor = ColorCode(lambda)% Return the color appropriate to the supplied wavelength.% Is it assumed the supplied lambda is within the range 380-780 nm% Smaller or higher values are set notionally to the extreme values % All input measurements are in meters.% This approximate conversion from nm to RGB is taken from thiscolor = [0,0,0];lambda    = lambda*1e+9;    % Convert to nm.if lambda<380,    thiscolor = [1,0,1]; end;if (lambda>=380)&(lambda<440),   thiscolor = [(440-lambda)/(440-380),0,1]; end;if (lambda>=440)&(lambda<490),   thiscolor = [0,(lambda-440)/(490-440),1]; end;if (lambda>=490)&(lambda<510),   thiscolor = [0,1,(510-lambda)/(510-490)]; end;if (lambda>=510)&(lambda<580),   thiscolor = [(lambda-510)/(580-510),1,0]; end;if (lambda>=580)&(lambda<645),   thiscolor = [1,(645-lambda)/(645-580),0]; end;if (lambda>=645),   thiscolor = [1,0,0]; end;%  The intensities fall off near limits of visionif lambda>700,   thiscolor = thiscolor * (0.3 + 0.7*(780-lambda)/(780-700)); end;if lambda<420,   thiscolor = thiscolor * (0.3 + 0.7*(lambda-380)/(420-380)); end;return;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 