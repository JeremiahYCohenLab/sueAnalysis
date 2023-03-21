function FigHandle = figure2Wide(varargin)
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
%   pos = get(FigH, 'Position');
  x = -1000;
  y = 500; 
  w = 900;
  h = 300;
  set(FigH, 'Position', [x,y, w, h], ...
            'Visible', paramVisible);
end
if nargout ~= 0
  FigHandle = FigH;
end