function FigHandle = figure2(varargin)
MP = get(0, 'MonitorPositions');
if size(MP, 1) == 1  % Single monitor
  FigH = figure(varargin{:});
else                 % Multiple monitors
  % Catch creation of figure with disabled visibility: 
  indexVisible = find(strncmpi(varargin(1:2:end), 'Vis', 3));
  if ~isempty(indexVisible)
    paramVisible = varargin(indexVisible(end) + 1);
  else
    paramVisible = get(0, 'DefaultFigureVisible');
  end
  %
  FigH     = figure(varargin{:}, 'Visible', 'off');
  set(FigH, 'Units', 'pixels');
  pos = get(FigH, 'Position');
  x = MP(1,1) + 0.5*MP(2,1) - 0.5*pos(3) ;
  y = MP(2,2) + 0.5*MP(1,4) - 0.5*pos(4) + 500; 
  set(FigH, 'Position', [x,y, pos(3:4)], ...
            'Visible', paramVisible);
end
if nargout ~= 0
  FigHandle = FigH;
end